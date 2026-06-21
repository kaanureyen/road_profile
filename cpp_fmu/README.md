# C++ Functional Mock-up Unit (FMU) Sources

This directory contains the native C++ FMI 2.0 Co-Simulation implementation of the deterministic road profile generator.

---

## 1. Directory Structure

* **`src/`**
  * [InfiniteRoadFMU.cpp](src/InfiniteRoadFMU.cpp): Main implementation file. Contains the co-simulation slave life-cycle entry points (`fmi2Instantiate`, `fmi2FreeInstance`, `fmi2GetReal`, `fmi2SetReal`, etc.), the custom NumPy-compatible MT19937 random number generator, and the vectorized height calculation loop.
* **`include/`**
  * FMI 2.0 headers standard files:
    * `fmi2Functions.h`
    * `fmi2FunctionTypes.h`
    * `fmi2TypesPlatform.h`
* **[CMakeLists.txt](CMakeLists.txt)**: Configures the MSVC compiler options, sets standard library linking configurations, and sets target properties to output `InfiniteRoadFMU.dll` without prefixes.

---

## 2. Compilation and Optimization Details

To meet real-time simulation requirements, the C++ code is highly optimized:

### A. AVX2 SIMD Vectorization
The core equation of the isotropic road model is a sum-of-sinusoids:
$$z(x, y) = \sum_i \text{amp}_i \cos(kx_i \cdot x + ky_i \cdot y + \phi_i)$$
Under MSVC compilation, the loop is auto-vectorized by enabling `/arch:AVX2`. This enables the compiler to use 256-bit AVX2 SIMD registers to process **8 floating-point wave evaluations concurrently** in a single CPU instruction pipeline.

### B. Fast Math & Minimax Cosine
Standard library `std::cos` executes as a double-precision scalar call with high precision but massive latency overhead (~10–12 ns). To resolve this:
1. **Float Conversion**: Internal math is executed in single-precision floating-point (`float`).
2. **Minimax Polynomial**: We implement a 6th-degree minimax polynomial cosine approximation evaluated branchlessly via Horner's method. This reduces the cosine evaluation to basic multiplications and additions that execute entirely within CPU registers.
3. **Speedup**: The combination of register math, single-precision conversion, and AVX2 vectorization yields a **28x speedup** on the raw mathematical computation loop compared to the original unoptimized scalar double-precision implementation (dropping from $169.3\ \mu\text{s}$ to $6.00\ \mu\text{s}$ per query for 16,384 wave components).

---

## 3. How to Build Manually

While the top-level script `build_cpp_fmu.py` automates compilation and packaging, you can configure and build the C++ shared library manually:

```bash
# Configure the project in a build directory
cmake -B build -S .

# Build the DLL in Release mode
cmake --build build --config Release
```

The compiled shared library will be generated at `build/Release/InfiniteRoadFMU.dll`.
