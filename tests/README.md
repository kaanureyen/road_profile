# Road Profile Generator Verification Tests

This directory contains verification tests, diagnostic reports, and plots validating the 2D isotropic road profile generator using the compiled FMI Co-Simulation binary `InfiniteRoadFMU.fmu`.

Every test suite is organized into its own subdirectory containing a `README.md` report that explains the test, presents the results, and displays the corresponding validation plots. This allows you to inspect all tests directly on GitHub.

---

## Test Directory Structure

### 1. [FMI Compliance & Concurrent Simulation](fmu_validation/)
- **Test:** [`test_fmu_simulation.py`](fmu_validation/test_fmu_simulation.py)
- **FMU Used:** **Yes** (interacts directly with the compiled `InfiniteRoadFMU.fmu` binary)
- **Report:** [FMI Validation Report](fmu_validation/README.md)
- **Validates:** Instantiation of multiple independent FMU instances, spatial coordinate querying safety, cross-instance determinism (exact height match), and seed sensitivity (different realization seeds).
- **Result Plot:** `tests/fmu_validation/fmu_simulation_results.png`

### 2. [Spatial Homogeneity Verification](distance_homogeneity/)
- **Test:** [`test_distance_homogeneity.py`](distance_homogeneity/test_distance_homogeneity.py)
- **FMU Used:** **Yes** (queries the `InfiniteRoadFMU.fmu` in parallel via `fmu_helper.py`)
- **Report:** [Distance Homogeneity Report](distance_homogeneity/README.md)
- **Validates:** Spatial homogeneity of the FMU isotropic wave field at large offsets from the origin (0 km, 1 km, 10 km, and 100 km) using parallel FMU query workers.
- **Result Plot:** `tests/distance_homogeneity/distance_homogeneity_curves.png`

### 3. [ISO 8608 Parameter Fitting & Calibration](parameter_fitting/)
- **Test:** [`test_parameter_fitting.py`](parameter_fitting/test_parameter_fitting.py)
- **FMU Used:** **Yes** (queries the `InfiniteRoadFMU.fmu` in parallel via `fmu_helper.py`)
- **Report:** [Parameter Fitting Report](parameter_fitting/README.md)
- **Validates:** Exponent $w$ and roughness $G$ estimation accuracy across Class B, Class C, and Class D road types. Evaluates isotropy and homogeneity using log-log cumulative PSD projection curve-fitting.
- **Result Plots:** 
  - `tests/parameter_fitting/parameter_fitting_case_*.png` (Individual cases)
  - `tests/parameter_fitting/parameter_fitting_summary.png` (Summary curve match)

### 4. [Power Spectral Density & Sensitivity sweeps](psd_analysis/)
- **Tests:**
  - [`plot_direct_fft.py`](psd_analysis/plot_direct_fft.py) (**FMU Used: Yes**, via `fmu_helper.py`)
  - [`plot_infinite_comparison.py`](psd_analysis/plot_infinite_comparison.py) (**FMU Used: Yes**, via `fmu_helper.py`)
  - [`plot_psd_comparison.py`](psd_analysis/plot_psd_comparison.py) (**FMU Used: Yes**, via `fmu_helper.py`)
  - [`run_advanced_analysis.py`](psd_analysis/run_advanced_analysis.py) (**FMU Used: Yes**, via `fmu_helper.py`)
  - [`run_advanced_sensitivity.py`](psd_analysis/run_advanced_sensitivity.py) (**FMU Used: Yes**, via `fmu_helper.py`)
- **Report:** [PSD Analysis Master Report](psd_analysis/README.md)
- **Validates:** Frequency grid rings ($N_f$) and angular division ($N_\theta$) discretization sensitivities on the actual FMU, periodogram discretization valleys, and infinite vs. band-limited model truncation.
- **Result Plots:** Multiple sensitivity curves (raw vs. cumulative PSD) saved in `tests/psd_analysis/`.

---

## Summary of Used FMU Parameters

The table below summarizes the exact FMI variables set on the FMU for each of the test scenarios:

| Test Suite / Scenario | `seed` | `road_class` | `Gd_n0` ($m^3$) | `w` | `f_min` (c/m) | `f_max` (c/m) | `Nf` | `Ntheta` |
|---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| **FMI Simulation (Wheel 1-4)** | `42` | `2` (Class B) | *Default (Inactive)* | `2.0` | `0.002` | `2000.0` | `512` | `32` |
| **FMI Simulation (Sensitivity)** | `99` | `2` (Class B) | *Default (Inactive)* | `2.0` | `0.002` | `2000.0` | `512` | `32` |
| **Distance Homogeneity** | `dist + slice_idx + 2026` | `0` (Custom) | `64e-6` | `2.0` | `0.002` | `2000.0` | `512` | `32` |
| **Parameter Fitting (Case 1)** | `200 + slice_idx` | `0` (Custom) | `64e-6` | `2.0` | `0.002` | `2000.0` | `512` | `32` |
| **Parameter Fitting (Case 2)** | `201 + slice_idx` | `0` (Custom) | `256e-6` | `1.8` | `0.002` | `2000.0` | `512` | `32` |
| **Parameter Fitting (Case 3)** | `202 + slice_idx` | `0` (Custom) | `1024e-6` | `2.2` | `0.002` | `2000.0` | `512` | `32` |
| **Direct FFT Comparison** | `42` | `0` (Custom) | `64e-6` | `2.0` | `0.002` | `20.0` | `64` | `16` |
| **Infinite Model Comparison** | `42` | `0` (Custom) | `64e-6` | `2.0` | `0.002` | `20.0` | `64` | `16` |
| **PSD Comparison** | `42` | `0` (Custom) | `64e-6` | `2.0` | `0.002` | `20.0` | `64` | `16` |
| **Parameter Sweep Case 1** | `42` | `0` (Custom) | `4e-6` | `1.5` | `0.002` | `2000.0` | `64` | `16` |
| **Parameter Sweep Case 2** | `42` | `0` (Custom) | `64e-6` | `2.0` | `0.002` | `2000.0` | `64` | `16` |
| **Parameter Sweep Case 3** | `42` | `0` (Custom) | `256e-6` | `3.0` | `0.002` | `2000.0` | `64` | `16` |
| **Parameter Sweep Case 4** | `42` | `0` (Custom) | `16384e-6` | `4.5` | `0.002` | `2000.0` | `64` | `16` |
| **Grid Sensitivity ($N_f$ Sweep)** | `42` | `0` (Custom) | `64e-6` | `2.0` | `0.002` | `2000.0` | `16` to `1024` | `16` |
| **Grid Sensitivity ($N_\theta$ Sweep)** | `42` | `0` (Custom) | `64e-6` | `2.0` | `0.002` | `2000.0` | `64` | `4` to `64` |

---

## Run All Verification Tests

Execute the runner script in the repository root to automatically run all the verification scripts in sequence and regenerate all output plots and markdown reports:

```bash
python run_tests.py
```
