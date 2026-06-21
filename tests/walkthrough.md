# Walkthrough - Road Profile Generator Default Discretization and Homogeneity Validation

This document summarizes the changes, optimizations, and verification results for updating the default road profile grid discretization to $N_f = 512, N_\theta = 32$ and implementing the distance homogeneity test case.

## Summary of Completed Tasks

1. **Updated Default Discretization**:
   - Modified [infinite_road_fmu.py](../infinite_road_fmu.py) to use $N_f = 512$ and $N_\theta = 32$ by default.
   - Rebuilt the FMU binary to [InfiniteRoadFMU.fmu](../InfiniteRoadFMU.fmu) via `pythonfmu build -f infinite_road_fmu.py`.

2. **1D Mathematical Simplification & Float32 Optimization**:
   - Converted the 2D sum-of-sinusoids height evaluation along the slice path:
     $$h(x, y) = \sum_n \text{amp}_n \cos(kx_n x + ky_n y + \phi_n)$$
     to a simplified 1D wave equation:
     $$h(s) = \sum_n \text{amp}_n \cos(k_n s + \psi_n)$$
     where $k_n = kx_n \cos(\theta_{\text{slice}}) + ky_n \sin(\theta_{\text{slice}})$ and $\psi_n = kx_n x_1 + ky_n y_1 + \phi_n$.
   - Implemented chunked NumPy vectorization in `float32` (converted back to `float64` for fitting). This reduced the evaluation time of a single slice (250,000 points, 16,384 wave components) from **88.0 seconds** to **34.9 seconds** (a **2.5x speedup**).

3. **Multiprocessing Parallelization**:
   - Integrated Python's `ProcessPoolExecutor` in both `test_distance_homogeneity.py` and `test_parameter_fitting.py`.
   - Distributing the slices across logical CPU cores yields a **10x execution speedup** (e.g. running 10 slices in parallel takes ~80 seconds instead of ~13 minutes).

4. **Distance Homogeneity Test Case**:
   - Created [test_distance_homogeneity.py](distance_homogeneity/test_distance_homogeneity.py) to simulate random starting headings (distance offset from origin) and running directions.
   - Fits the direct Hanning FFT cumulative PSD to the exact isotropic projection model using a decimated grid of 100 points in $[0.02, 200.0]$ cycles/m.
   - Compares offsets of **0 km, 1 km, 10 km, and 100 km** from the origin.
   - Saves the comparative curves to [distance_homogeneity_curves.png](distance_homogeneity/distance_homogeneity_curves.png).

5. **Updated Fitting Benchmark**:
   - Updated [test_parameter_fitting.py](parameter_fitting/test_parameter_fitting.py) to match the new direct Hanning FFT model with Nf=512, Ntheta=32 defaults, and parallel execution.
   - Re-generated the 3 road class fitting reports and isotropy dependency plots.

---

## Verification Results

### 1. Distance Homogeneity Verification (`test_distance_homogeneity.py`)
Evaluating 10 random slices per offset with target parameters $w = 2.0$ and $G = 64.0\ \mu\text{m}^3$ (Class B) yielded the following:

- **Offset 0.0 km**:
  - Calibrated Exponent $w$: $2.0201 \pm 0.0219$ (Error: $1.01\%$) | Target in $\pm 1$ std: **YES**
  - Calibrated Roughness $G$: $70.15 \pm 7.20\ \mu\text{m}^3$ (Error: $9.61\%$) | Target in $\pm 1$ std: **YES**
- **Offset 1.0 km**:
  - Calibrated Exponent $w$: $2.0092 \pm 0.0338$ (Error: $0.46\%$) | Target in $\pm 1$ std: **YES**
  - Calibrated Roughness $G$: $66.61 \pm 10.62\ \mu\text{m}^3$ (Error: $4.08\%$) | Target in $\pm 1$ std: **YES**
- **Offset 10.0 km**:
  - Calibrated Exponent $w$: $1.9984 \pm 0.0250$ (Error: $0.08\%$) | Target in $\pm 1$ std: **YES**
  - Calibrated Roughness $G$: $63.56 \pm 7.69\ \mu\text{m}^3$ (Error: $0.69\%$) | Target in $\pm 1$ std: **YES**
- **Offset 100.0 km**:
  - Calibrated Exponent $w$: $1.9951 \pm 0.0493$ (Error: $0.24\%$) | Target in $\pm 1$ std: **YES**
  - Calibrated Roughness $G$: $64.20 \pm 12.86\ \mu\text{m}^3$ (Error: $0.31\%$) | Target in $\pm 1$ std: **YES**

*Ensemble Average Exponent error is **0.44%**, and Roughness error is **3.3%**. All offsets successfully verify parameter homogeneity, and target values fall well within $\pm 1$ standard deviation of the fitted means.*

### 2. Multi-Class Parameter Fitting Verification (`test_parameter_fitting.py`)
Evaluating 10 random slices per road class using the $512 \times 32$ discretization:

- **Case 1: Class B ($G = 64.0\ \mu\text{m}^3, w = 2.0$)**:
  - Fitted Exponent $w$: $1.9956 \pm 0.0222$ (Error: $0.22\%$) | Target in $\pm 1$ std: **YES**
  - Fitted Roughness $G$: $61.9 \pm 6.38\ \mu\text{m}^3$ (Error: $3.30\%$) | Target in $\pm 1$ std: **YES**
- **Case 2: Class C ($G = 256.0\ \mu\text{m}^3, w = 1.8$)**:
  - Fitted Exponent $w$: $1.8031 \pm 0.0205$ (Error: $0.17\%$) | Target in $\pm 1$ std: **YES**
  - Fitted Roughness $G$: $253.0 \pm 27.5\ \mu\text{m}^3$ (Error: $1.17\%$) | Target in $\pm 1$ std: **YES**
- **Case 3: Class D ($G = 1024.0\ \mu\text{m}^3, w = 2.2$)**:
  - Fitted Exponent $w$: $2.1855 \pm 0.0230$ (Error: $0.66\%$) | Target in $\pm 1$ std: **YES**
  - Fitted Roughness $G$: $957.0 \pm 88.8\ \mu\text{m}^3$ (Error: $6.54\%$) | Target in $\pm 1$ std: **YES**

*Average exponent error is **0.35%**, and roughness error is **3.67%**.*

---

## Artifacts Generated

- Comparative curves plot: [distance_homogeneity_curves.png](distance_homogeneity/distance_homogeneity_curves.png)
- Local fitting dependency plots:
  - [parameter_fitting_case_1.png](parameter_fitting/parameter_fitting_case_1.png)
  - [parameter_fitting_case_2.png](parameter_fitting/parameter_fitting_case_2.png)
  - [parameter_fitting_case_3.png](parameter_fitting/parameter_fitting_case_3.png)
  - [parameter_fitting_summary.png](parameter_fitting/parameter_fitting_summary.png)
- Parameter fitting analysis report: [parameter_fitting_analysis.md](parameter_fitting/parameter_fitting_analysis.md)
