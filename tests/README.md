# Road Profile Generator Verification Tests

This directory contains verification tests, diagnostic reports, and plots validating the 2D isotropic road profile generator. 

The codebase separates verification into two distinct categories: **Analytical Python Simulations** and **Compiled FMU Binary Tests**.

---

## 1. Analytical Python Simulations
These tests evaluate the mathematical correctness and statistical isotropy/homogeneity of the 2D sum-of-sinusoids algorithms directly in Python (without compiling or loading the FMU binary).

*   **`test_distance_homogeneity.py`**:
    *   **Evaluates**: Spatial homogeneity at offsets of 0 km, 1 km, 10 km, and 100 km from the origin.
    *   **Method**: Generates random start angles and slice headings using analytical sum-of-sinusoids waves, fits the cumulative PSD to the exact isotropic projection model, and verifies that the fitted parameters fall within $\pm 1$ standard deviation of the targets.
    *   **Output**: Plots saved to `tests/distance_homogeneity_curves.png`.
*   **`test_parameter_fitting.py`**:
    *   **Evaluates**: Enexponent $w$ and roughness $G$ fitting errors across three standard road classes (Class B, Class C, Class D).
    *   **Method**: Generates 10 random slices per class, applies direct Hanning FFT cumulative projection fitting, corrects windowing/discretization bias using calibrated mappings, and outputs statistical confidence limits.
    *   **Output**: Saved to `tests/parameter_fitting_case_*.png`, `tests/parameter_fitting_summary.png`, and a markdown report `tests/parameter_fitting_analysis.md`.
*   **`test_fmu_isotropy.py`**:
    *   **Evaluates**: Isotropic direction independence. It queries slices in different angles in Python to verify that direction does not bias fitted exponent and roughness values.
*   **`test_perfect_fmu.py`**:
    *   **Evaluates**: Fits target cumulative PSD curves to coordinate profiles.

---

## 2. Compiled FMU Binary Tests (FMI Compliance)
These tests interact directly with the compiled FMI Co-Simulation binary `InfiniteRoadFMU.fmu` to verify integration correctness and FMI standard compliance.

*   **`test_fmu_simulation.py`**:
    *   **Evaluates**: Core FMI compliance, multi-wheel concurrent query safety, and determinism.
    *   **Method**: Uses the `fmpy` library to unzip, extract, and instantiate **4 concurrent instances** of the compiled FMU binary (`InfiniteRoadFMU.fmu`).
    *   **Checks**:
        1.  *Concurrent queries*: Asserts independent coordinate querying on separate FMU instances representing different wheels.
        2.  *determinism & Repeatability*: Queries the same coordinate across all 4 FMU instances and asserts that they return identical height values (max diff = `0.00e+00` m).
        3.  *Seed sensitivity*: Instantiates FMUs with different seed parameters (`seed=42` and `seed=99`) and asserts that they generate independent, distinct random realizations of the road surface.

---

## 3. Automation Scripts
*   **`run_tests.py`** (located in the repository root):
    *   Runs the primary analytical validation suites (`test_distance_homogeneity.py` and `test_parameter_fitting.py`) in parallel using a multiprocessing executor.
    *   Saves the resulting diagnostic plots and markdown reports directly into this `tests/` directory for immediate viewing on GitHub.
