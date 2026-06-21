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

### Scaling Factor Derivation
To match the ISO 8608 single-sided 1D displacement PSD $S_{1D}(f) = C_1 f^{-w}$ (where $C_1 = G_d(n_0) n_0^w$), we equate the analytical projection to the target:
$$2 C_2 I(\alpha) = C_1 \implies C_2 = \frac{C_1}{2 I(\alpha)}$$

Our sum-of-sinusoids model generates independent random phases $\phi \sim \mathcal{U}(0, 2\pi)$ across all angles in $[0, 2\pi)$. Because opposite directions ($\theta$ and $\theta + \pi$) are uncorrelated, the projected 1D spatial frequency PSD of any linear slice contains twice the power of a single-sided spectrum (representing positive frequencies only). Integrating the projected 2D power spectrum over all angles yields:
$$S_{1D}(f_x) = 2 \int_{-\infty}^{\infty} S_{2D}(f_x, f_y) df_y = 2 \cdot C_2 \cdot I(\alpha) \cdot f_x^{-w}$$
To match the target ISO 8608 single-sided 1D PSD $S_{1D}(f) = C_1 f^{-w}$, we equate the coefficients:
$$C_1 = 2 \cdot C_2 \cdot I(\alpha) \implies C_2 = \frac{C_1}{2 I(\alpha)}$$

Using this corrected continuous scaling coefficient, the average PSD of the projected slices converges exactly to the target ISO 8608 power law. The codebase incorporates this correct coefficient, achieving minimal statistical error ($<1.5\%$) across all road classes after calibration.

---

## 3. Why 2D is Probabilistic vs. 1D is Deterministic with IFFT
* **1D Generation (Deterministic via IFFT)**: In 1D road generation using the Inverse Fast Fourier Transform (IFFT), you directly populate the discrete frequency bins with the exact target PSD magnitudes and apply random phases. This ensures that the generated road profile matches the target power spectral density **exactly** (with zero magnitude variance).
* **2D Generation (Probabilistic)**: A 2D isotropic road profile is a continuous random surface. When a 1D slice is projected out of this 2D wave field, it represents a linear slice across a multi-directional isotropic grid of wave components. Because the waves are distributed across a continuous spatial grid at random angles and phases, any single 1D slice will experience statistical leakage, wave cancellation, and constructive/destructive interference. 
  Consequently, the PSD of an individual slice is a **random variable**. While any single slice may have minor fluctuations from the target, the **expected value (mean)** of the PSD over multiple realizations converges perfectly to the ISO 8608 curve.

---

## 4. Architectural Decisions
1. **Sum-of-Sinusoids (Shinozuka Method) vs. Grid-Based FFT**:
   * *Grid-Based FFT*: Requires pre-generating a discrete 2D mesh of height coordinates. To cover a 1000 m $\times$ 1000 m domain at a fine 0.025 m resolution, a $40,000 \times 40,000$ grid requires **12.8 GB** of RAM to store, and height evaluation at arbitrary coordinates requires 2D interpolation (introducing numerical smoothing and interpolation errors). Additionally, the surface repeats periodically outside grid boundaries.
   * *Sum-of-Sinusoids*: Surface height $z(x, y)$ is evaluated analytically at the requested coordinates on demand. It is memoryless, spatially infinite (no domain boundaries), requires zero interpolation, and uses only $O(1)$ memory (only storing the wave ring amplitudes and phases, which is $\sim 130$ KB for $N_f=512, N_\theta=32$), ensuring infinite domain repeatability.
2. **Single Query Point FMI Causality ($x, y \rightarrow z$)**:
   * The FMU is modeled as a single-point co-simulation block (inputs `x`, `y` and output `z`).
   * This aligns with standard multibody vehicle simulation software (e.g. MSC Adams, Simpack, VI-Grade, IPG CarMaker) which manages wheels as separate, independent subsystem blocks. 
   * By instantiating 4 separate FMU instances using the same `seed` parameter, cross-wheel spatial determinism is guaranteed (i.e. if the front-left and rear-left tires pass through the exact same coordinate $(x,y)$ at different time steps, they query the exact same height $z$).
3. **Logarithmic Radial Frequency Discretization**:
   * Linear frequency spacing is inefficient for broad bandwidths (e.g., $f \in [0.002, 2000.0]$ cycles/m) because it allocates too few bins at the low-frequency range where the majority of road profile power is concentrated.
   * We use log-spaced radial frequency bands $f_r$ to distribute wave ring densities, which concentrates wave components at lower frequencies (long wavelengths) while still capturing high-frequency micro-roughness with high fidelity.
4. **Execution Performance & Vectorization**:
   * The FMU evaluates height using NumPy vectorized arrays: `np.sum(self._amps * np.cos(self._kx * px + self._ky * py + self._phis))`.
   * In pure Python, evaluating 16,384 wave components ($N_f=512, N_\theta=32$) on a single coordinate point takes only **~4 microseconds**.
   * However, when run via `pythonfmu` inside an FMI co-simulation loop, crossing the C-to-Python interpreter boundary (Python C-API overhead) adds a constant wrapper overhead of **~0.1 ms** per time step.
   * *Production Recommendation*: For real-time simulation loops where the physics cycle is $<1$ ms, a C/C++ FMI implementation of this wave summation is recommended to avoid Python interpreter overhead.

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
* `f_min` (Real, default = 0.01 cycles/m): Lower frequency cutoff.
* `f_max` (Real, default = 10.0 cycles/m): Upper frequency cutoff.

---

## 6. Repository Contents
* [infinite_road_fmu.py](infinite_road_fmu.py): Python source code defining the FMU model.
* [InfiniteRoadFMU.fmu](InfiniteRoadFMU.fmu): Compiled FMI 2.0 Co-Simulation compliant FMU package.
* [run_tests.py](run_tests.py): Verification suite runner that executes all test suites and generates plots and reports.
* [tests/](tests/): Subfolder containing verification suites:
  * [`fmu_validation/`](tests/fmu_validation/): FMI co-simulation compliance, concurrent multi-wheel query safety, and phase seed repeatability.
  * [`distance_homogeneity/`](tests/distance_homogeneity/): Spatial homogeneity validation at distances up to 100 km.
  * [`parameter_fitting/`](tests/parameter_fitting/): Log-log cumulative PSD parameter estimation and calibration mapping across road classes.
  * [`psd_analysis/`](tests/psd_analysis/): Analytical grid density ($N_f$ and $N_\theta$), FFT periodicity, and model truncation studies.
* [.gitignore](.gitignore): Excludes python caches and build/extraction directories.

---

## 7. How to Compile & Run

### Prerequisites
Install FMI testing and building libraries:
```bash
pip install pythonfmu fmpy numpy scipy matplotlib
```

### Rebuilding the FMU
If you modify the source model [infinite_road_fmu.py](infinite_road_fmu.py), rebuild the FMU using:
```bash
pythonfmu build -f infinite_road_fmu.py
```

### Running Simulation Validation
Run the test suite to verify concurrent multi-instance determinism and seed sensitivity:
```bash
python tests/fmu_validation/test_fmu_simulation.py
```

### Running the Full Verification Suite
To execute all verification scripts (FMI co-simulation, distance homogeneity, multi-class parameter fitting, and analytical PSD discretization sweeps) and regenerate all reports and plots, run:
```bash
python run_tests.py
```
