# 2D Infinite Isotropic Road Profile Generator FMU (ISO 8608 Compliant)

This repository contains a Functional Mock-up Unit (FMU) that dynamically generates a 2D infinite, isotropic, and deterministic road profile following the **ISO 8608** standard. 

The FMU is compliant with the **FMI 2.0 Co-Simulation standard** and is simplified to a single query point $(x, y) \rightarrow z$. To simulate a full 4-wheel vehicle, the multibody simulator instantiates 4 separate instances of this FMU (one for each tire contact patch) using the same `seed` parameter, ensuring perfect cross-wheel spatial determinism.

---

## 1. ISO 8608 Standard Overview
The **ISO 8608** standard classifies road profiles based on their Power Spectral Density (PSD) of vertical displacement. The 1D spatial frequency PSD $G_d(n)$ is modeled using a power-law relationship:

$$G_d(n) = G_d(n_0) \left( \frac{n}{n_0} \right)^{-w}$$

Where:
- $n$: Spatial frequency (cycles/m).
- $n_0 = 0.1$ cycles/m: Reference spatial frequency.
- $w = 2.0$: Spatial frequency exponent (controlling the slope of the power spectrum).
- $G_d(n_0)$: Reference displacement PSD (roughness coefficient) at $n_0$, defining the road class:

| Road Class | Description | $G_d(n_0) \times 10^{-6}$ ($\text{m}^3$) |
| :--- | :--- | :---: |
| **Class A** | Very good | $16$ |
| **Class B** | Good | $64$ |
| **Class C** | Average | $256$ |
| **Class D** | Poor | $1024$ |
| **Class E** | Very poor | $4096$ |

---

### 2. Mathematical Derivation of the 2D Isotropic Wave Field
To generate a 2D isotropic surface where a 1D slice in *any* direction matches the ISO 8608 power law, we utilize the **Spectral Representation Method (Sum-of-Sinusoids)**.

### Continuous 2D PSD Projection
In 2D isotropic coordinate space, spatial frequency components are defined as:
$$f_x = f_r \cos\theta, \quad f_y = f_r \sin\theta$$
where $f_r$ is the radial spatial frequency and $\theta$ is the spatial angle. 

A 1D slice along any direction (e.g., $x$) corresponds to integrating the 2D PSD $S_{2D}(f_x, f_y)$ over the transverse frequency $f_y$:
$$S_{1D}(f_x) = 2 \int_{-\infty}^{\infty} S_{2D}(f_x, f_y) df_y$$

Assuming isotropy, $S_{2D}(f_x, f_y) = S_{2D}(f_r)$. Setting $S_{2D}(f_r) = C_2 f_r^{-\alpha}$ (where $\alpha = w + 1 = 3.0$), we make the variable substitution $f_y = f_x t$, which yields:
$$S_{1D}(f_x) = 2 \int_{-\infty}^{\infty} C_2 (f_x^2 + f_x^2 t^2)^{-\alpha/2} f_x dt = 2 C_2 f_x^{-(\alpha-1)} \int_{-\infty}^{\infty} (1+t^2)^{-\alpha/2} dt$$

This simplifies to:
$$S_{1D}(f_x) = 2 C_2 f_x^{-w} I(\alpha)$$
where $I(\alpha) = \int_{-\infty}^{\infty} (1+t^2)^{-\alpha/2} dt$ is a standard numerical integral.

### Scaling Factor and Half-Circle Integration
To match the target ISO 8608 single-sided 1D displacement PSD $S_{1D}(f) = C_1 f^{-w}$ (where $C_1 = G_d(n_0) n_0^w$), the physical relation requires $2 C_2 I(\alpha) = C_1 \implies C_2 = \frac{C_1}{2 I(\alpha)}$ when integrating over the entire $2\pi$ circle.

However, since a wave propagating at angle $\theta + \pi$ is spatially collinear with a wave at $\theta$ (just with an independent random phase), the full-circle $[0, 2\pi)$ formulation contains redundant propagation axes. 

By restricting the angular discretization to the half-circle $[0, \pi)$, we eliminate this redundancy. To preserve the total surface variance and target PSD, the continuous power spectral density is doubled, resulting in the continuous scaling coefficient:
$$C_2' = \frac{C_1}{I(\alpha)}$$

Our sum-of-sinusoids model generates independent random phases $\phi \sim \mathcal{U}(0, 2\pi)$ and distributes the spatial angles uniformly over the range $[0, \pi)$. Because opposite directions are omitted, each of the $N_\theta$ angular bins represents a unique wave propagation axis. This doubles the angular resolution for a given $N_\theta$, making the generated 2D surface more isotropic and reducing directional bias.

Using this updated scaling coefficient, the average PSD of the projected slices converges exactly to the target ISO 8608 power law. The codebase incorporates this correct coefficient, achieving minimal statistical error ($<1.5\%$) across all road classes after calibration.

---

## 3. Why 2D is Probabilistic vs. 1D is Deterministic with IFFT
* **1D Generation (Deterministic via IFFT)**: In 1D road generation using the Inverse Fast Fourier Transform (IFFT), you directly populate the discrete frequency bins with the exact target PSD magnitudes and apply random phases. This ensures that the generated road profile matches the target power spectral density **exactly** (with zero magnitude variance).
* **2D Generation (Probabilistic)**: A 2D isotropic road profile is a continuous random surface. When a 1D slice is projected out of this 2D wave field, it represents a linear slice across a multi-directional isotropic grid of wave components. Because the waves are distributed across a continuous spatial grid at random angles and phases, any single 1D slice will experience statistical leakage, wave cancellation, and constructive/destructive interference. 
  Consequently, the PSD of an individual slice is a **random variable**. While any single slice may have minor fluctuations from the target, the **expected value (mean)** of the PSD over multiple realizations converges perfectly to the ISO 8608 curve.

---

### 4. Architectural Decisions & Optimization

1. **Sum-of-Sinusoids (Shinozuka Method) vs. Grid-Based FFT & Caching**:
   * *Grid-Based FFT & Caching*: Pre-generating a discrete 2D mesh of height coordinates requires immense memory (e.g., 12.8 GB for a 1 km x 1 km grid at 0.025 m resolution). Alternatively, dynamically generating and caching 2D grid blocks on the fly to follow the vehicle introduces severe code complexity, boundary-transition interpolation logic, and execution latency spikes. Furthermore, if a vehicle spins, slides, or travels in arbitrary non-linear trajectories, tire contact patches sweep unpredictable paths across the 2D plane, rendering narrow path pre-calculation impossible.
   * *Analytical Procedural Math*: Height $z(x, y)$ is evaluated analytically at the exact tire coordinates on demand. It is truly memoryless ($O(1)$ RAM footprint, storing only $\sim 130$ KB of wave parameter states), spatially infinite (no domain or boundary limits), and requires zero grid-paging or boundary checks. Even if the vehicle spins, drifts, or takes arbitrary 2D trajectories, the height is computed instantly in $\approx 6\ \mu\text{s}$ with zero boundary-crossing or caching overhead.
2. **Native C++ Implementation (FMI 2.0 Compliance)**:
   * To prevent the C-to-Python ctypes wrapper overhead (~170 $\mu$s per time step), the production FMU is implemented in native C++ ([cpp_fmu/src/InfiniteRoadFMU.cpp](cpp_fmu/src/InfiniteRoadFMU.cpp)). This makes it 100% self-contained and Python-independent.
3. **Float-Precision Conversion**:
   * Replaced internal double-precision math with single-precision floating-point (`float`). This halves memory bandwidth and allows the CPU to process twice as many calculations per SIMD register.
4. **Minimax Cosine Polynomial & AVX2 SIMD Autovectorization**:
   * Standard library `std::cos` calls are slow when executed sequentially. Instead, we use a branchless 6th-degree minimax polynomial approximation evaluated via Horner's method. 
   * Compiled with `/arch:AVX2 /fp:fast` flags, the MSVC compiler auto-vectorizes this register-only loop to calculate 8 cosine values in parallel per clock cycle, accelerating computation by **28x** over scalar execution.

---

## 5. Performance Benchmarks

Evaluating a query of **25,000 points** along a road slice (16,384 wave components per query) yields the following performance comparison:

| Importer Wrapper | FMU Implementation | Total Time (s) | Avg Query Time ($\mu\text{s}$/point) | Speedup |
| :--- | :--- | :---: | :---: | :---: |
| **C++ Wrapper** | **C++ FMU (Optimized Minimax)** | **0.158 s** | **6.33 $\mu\text{s}$** | **31.3x** (Best) |
| **Python (FMPy)** | **C++ FMU (Optimized Minimax)** | 0.293 s | 11.73 $\mu\text{s}$ | 16.9x |
| **Python (FMPy)** | **Python FMU** | 4.969 s | 198.76 $\mu\text{s}$ | 1.00x (Baseline) |
| **C++ Wrapper** | **Python FMU** | 4.976 s | 199.04 $\mu\text{s}$ | 1.00x |

*Note: With dynamic memory lookup tables (LUTs) disabled, the C++ FMU's execution time is 99.88% dominated by pure register arithmetic, with FMI wrapper overhead contributing only **7.3 nanoseconds** (0.11%) per query.*

---

## 5. FMI Interface Specification

### Inputs
* `x` (Real): Longitudinal coordinate (m)
* `y` (Real): Lateral coordinate (m)

### Outputs
* `z` (Real): Road profile vertical height (m)

### Parameters
* `seed` (Integer, default = 42): Pseudorandom seed for phase generation.
* `road_class` (Integer, default = 3): ISO 8608 Class (1=A, 2=B, 3=C, 4=D, 5=E, 0=Custom Gd_n0).
* `Gd_n0` (Real, default = 256e-6 $\text{m}^3$): Reference displacement PSD at $n_0=0.1$ cycles/m (active when `road_class=0`).
* `w` (Real, default = 2.0): Spectral exponent.
* `f_min` (Real, default = 0.002 cycles/m): Lower frequency cutoff.
* `f_max` (Real, default = 2000.0 cycles/m): Upper frequency cutoff.

---

## 6. Repository Contents
* [infinite_road_fmu.py](infinite_road_fmu.py): Python source code defining the FMU model.
* [InfiniteRoadFMU.fmu](InfiniteRoadFMU.fmu): Compiled FMI 2.0 Co-Simulation compliant FMU package.
* [run_tests.py](run_tests.py): Verification suite runner that executes all test suites and generates plots and reports.
* [tests/](tests/): Subfolder containing verification suites:
  * [`fmu_validation/`](tests/fmu_validation/): FMI co-simulation compliance, concurrent multi-wheel query safety, and phase seed repeatability.
  * [`distance_homogeneity/`](tests/distance_homogeneity/): Spatial homogeneity validation at distances up to 100 km.
  * [`parameter_fitting/`](tests/parameter_fitting/): Log-log cumulative PSD parameter estimation and calibration mapping across road classes.
  * [`psd_analysis/`](tests/psd_analysis/): Welch-averaged PSD representation and exact isotropic cumulative model validation.
* [.gitignore](.gitignore): Excludes python caches and build/extraction directories.

---

## 7. How to Compile & Run

### Prerequisites
To compile the C++ FMU and execute verification tests, install:
* **C++ Compiler:** Visual Studio 2022 (MSVC) with C++ Desktop Development tools.
* **Build System:** CMake 3.10 or higher.
* **Python Environment:** Install testing libraries:
  ```bash
  pip install pythonfmu fmpy numpy scipy matplotlib
  ```

### Rebuilding the C++ FMU
To compile the C++ shared library (`InfiniteRoadFMU.dll`) and package it into `InfiniteRoadFMU.fmu`, run:
```bash
python build_cpp_fmu.py
```
This script will:
1. Configure and run a CMake build in the `cpp_fmu/build/` directory in Release mode.
2. Compile `InfiniteRoadFMU.cpp` with `/arch:AVX2 /fp:fast` optimizations.
3. Stage the compiled DLL, create an FMI 2.0-compliant `modelDescription.xml`, and package them into a compressed `.fmu` archive at the workspace root.

### Running the Full Verification Suite
To execute all verification scripts (including FMI co-simulation, 100 km homogeneity verification, parameter fitting, and PSD/terrain analysis) and regenerate all plots, run:
```bash
python run_tests.py
```
