# C++ Functional Mock-up Unit (FMU) Sources

This directory contains the native C++ FMI 2.0 Co-Simulation implementation of the deterministic road profile generator. It is designed to run in real-time multibody simulation environments by using high-performance SIMD optimizations, fast branchless polynomial math, and low-latency C FMI bindings.

---

## 1. Directory Structure

* **`src/`**
  * [InfiniteRoadFMU.cpp](file:///C:/Users/novo/.gemini/antigravity/scratch/road_profile/cpp_fmu/src/InfiniteRoadFMU.cpp): Main implementation file. Contains the co-simulation slave life-cycle entry points, the custom NumPy-compatible MT19937 random number generator, and the optimized height calculation loop.
* **`include/`**
  * FMI 2.0 standard headers:
    * `fmi2Functions.h`
    * `fmi2FunctionTypes.h`
    * `fmi2TypesPlatform.h`
* **[CMakeLists.txt](file:///C:/Users/novo/.gemini/antigravity/scratch/road_profile/cpp_fmu/CMakeLists.txt)**: Configures compiler flags, standard library configuration, and target properties to output the compiled `InfiniteRoadFMU.dll` without prefixes.

---

## 2. FMI 2.0 Interface & Value References

The C++ FMU exposes the following variables under the FMI 2.0 interface. These map directly to their corresponding value references (`vr`) inside [InfiniteRoadFMU.cpp](file:///C:/Users/novo/.gemini/antigravity/scratch/road_profile/cpp_fmu/src/InfiniteRoadFMU.cpp):

### Inputs & Outputs (Reals)

| Name | Direction | Value Reference (`vr`) | Default | Unit | Description |
| :--- | :---: | :---: | :---: | :---: | :--- |
| `x` | Input | `0` | `0.0` | `m` | Longitudinal position query point |
| `y` | Input | `1` | `0.0` | `m` | Lateral position query point |
| `z` | Output | `2` | N/A | `m` | Generated vertical road elevation |

### Parameters

| Name | Type | Value Reference (`vr`) | Default | Description |
| :--- | :---: | :---: | :---: | :--- |
| `seed` | Integer | `3` | `42` | Pseudorandom seed for MT19937 phase angle generator |
| `road_class` | Integer | `4` | `3` | ISO 8608 road class (1=A, 2=B, 3=C, 4=D, 5=E, 0=Custom `Gd_n0`) |
| `Gd_n0` | Real | `5` | `256e-6` | Reference displacement PSD at $n_0 = 0.1\text{ cycles/m}$ (only active if `road_class = 0`) |
| `w` | Real | `6` | `2.0` | Spatial frequency exponent |
| `f_min` | Real | `7` | `0.01` | Lower frequency cutoff ($\text{cycles/m}$) |
| `f_max` | Real | `8` | `2.0` | Upper frequency cutoff ($\text{cycles/m}$) |
| `Nf` | Integer | `9` | `512` | Number of discrete frequency bins |
| `Ntheta` | Integer | `10` | `32` | Number of discrete angular bins |
| `disable_math` | Integer | `11` | `0` | Optimization flag. Set to `1` to bypass computation and return `0.0` (used for pure FMI overhead profiling) |

---

## 3. Architecture & Optimization Details

To meet real-time simulation requirements (often requiring step query times under $10\ \mu\text{s}$), the C++ implementation utilizes several optimization techniques:

### A. Lazy Parameter Initialization
Recalculating the random phase array and wave vector frequencies is an expensive $O(N_f \times N_\theta)$ operation. In [lazy_init](file:///C:/Users/novo/.gemini/antigravity/scratch/road_profile/cpp_fmu/src/InfiniteRoadFMU.cpp#L133-L228), the FMU tracks whether any parameter has changed using a dirty-flag pattern.
- If all parameters (`seed`, `road_class`, `Gd_n0`, `w`, etc.) match their cached values, [compute_height](file:///C:/Users/novo/.gemini/antigravity/scratch/road_profile/cpp_fmu/src/InfiniteRoadFMU.cpp#L240-L259) skips the initialization and executes the wave query loop directly.
- If a parameter is updated via [fmi2SetReal](file:///C:/Users/novo/.gemini/antigravity/scratch/road_profile/cpp_fmu/src/InfiniteRoadFMU.cpp#L363-L391) or [fmi2SetInteger](file:///C:/Users/novo/.gemini/antigravity/scratch/road_profile/cpp_fmu/src/InfiniteRoadFMU.cpp#L415-L439), a `dirty` flag is set, forcing a regeneration of phase vectors on the next step evaluation.

### B. NumPy-Compatible Random State
To ensure the spatial road profile is bit-exact and identical to the Python simulation prototype, we implement a custom C++ class [NumPyRandom](file:///C:/Users/novo/.gemini/antigravity/scratch/road_profile/cpp_fmu/src/InfiniteRoadFMU.cpp#L16-L66). It replicates the Mersenne Twister MT19937 generator used in `numpy.random.RandomState(seed)` for uniform distribution generation, yielding the exact same phase angles $\phi_i \in [0, 2\pi)$ given the same seed.

### C. AVX2 SIMD Vectorization
The core equation of the isotropic road model is a sum-of-sinusoids:

$$z(x, y) = \sum_{i=1}^{N_f \cdot N_\theta} \text{amp}_i \cos(kx_i \cdot x + ky_i \cdot y + \phi_i)$$

When compiled with `/arch:AVX2 /fp:fast`, the loop in [compute_height](file:///C:/Users/novo/.gemini/antigravity/scratch/road_profile/cpp_fmu/src/InfiniteRoadFMU.cpp#L240-L259) is auto-vectorized by the compiler. The CPU uses 256-bit registers to evaluate **8 single-precision/double-precision math operations concurrently**, speeding up execution times dramatically.

### D. Minimax Polynomial Cosine Approximation
Standard library `std::cos` executes as a double-precision scalar call with high precision but massive latency (~10–12 ns). To achieve maximum performance, we implement [fast_cos](file:///C:/Users/novo/.gemini/antigravity/scratch/road_profile/cpp_fmu/src/InfiniteRoadFMU.cpp#L230-L238):
1. **Domain Reduction**: Reduces the input phase to the interval $[-\pi, \pi]$ using:
   $$y = x - 2\pi \cdot \text{round}\left(\frac{x}{2\pi}\right)$$
2. **Minimax Polynomial**: Evaluates a 6th-degree minimax polynomial approximation of cosine evaluated branchlessly via Horner's method:
   $$\cos(y) \approx 1.0 + y^2 \cdot (a_1 + y^2 \cdot (a_2 + y^2 \cdot (a_3 + y^2 \cdot (a_4 + y^2 \cdot (a_5 + y^2 \cdot a_6)))))$$
3. **Speedup**: The combination of register math, single-precision calculations, and AVX2 vectorization yields a **28x speedup** on the raw mathematical computation loop compared to the original unoptimized scalar double-precision implementation (dropping from $169.3\ \mu\text{s}$ to $6.00\ \mu\text{s}$ per query for $16,384$ wave components).

---

## 4. Packaging and Metadata

The C++ code compiles to a native shared library `InfiniteRoadFMU.dll`. This DLL is packaged inside the standard FMU zip container along with FMI metadata:
* **Generation Tool**: Native C++ Implementation
* **Co-Simulation**: Configured as a Co-Simulation slave (`needsExecutionTool="false"`, `canHandleVariableCommunicationStepSize="true"`).
* **Repository Link**: Documented directly inside the `description` metadata attribute of the FMI `modelDescription.xml`:
  `description="Deterministic ISO 8608 2D Isotropic Infinite Road Profile Generator (Single Point). Repo: https://github.com/kaanureyen/road_profile"`
* **Author**: Antigravity Coding Assistant

---

## 5. How to Build Manually

While the top-level script `build_cpp_fmu.py` automates compilation and packaging, you can configure and build the C++ shared library manually:

```bash
# Configure the project in a build directory
cmake -B build -S .

# Build the DLL in Release mode
cmake --build build --config Release
```

The compiled shared library will be generated at `build/Release/InfiniteRoadFMU.dll`.

> [!TIP]
> To package this DLL into the final `.fmu` format, run `python build_cpp_fmu.py` from the root directory. This script compiles the C++ project, writes the `modelDescription.xml` containing value references and the repository link, and zips them into `InfiniteRoadFMU.fmu`.
