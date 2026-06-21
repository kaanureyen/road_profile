import os
import sys
import subprocess
import time

def run_script(script_path):
    print(f"\n========================================================")
    print(f"Running {script_path}...")
    print(f"========================================================")
    t0 = time.time()
    
    script_dir = os.path.dirname(script_path)
    script_name = os.path.basename(script_path)
    cwd = os.path.abspath(script_dir) if script_dir else None
    
    # Run with unbuffered output so logs print live, specifying cwd
    res = subprocess.run([sys.executable, "-u", script_name], cwd=cwd, capture_output=False)
    
    elapsed = time.time() - t0
    if res.returncode == 0:
        print(f"SUCCESS: {script_path} finished in {elapsed:.1f}s")
        return True
    else:
        print(f"FAILED: {script_path} exited with code {res.returncode} in {elapsed:.1f}s")
        return False

def main():
    print("=== ROAD PROFILE GENERATOR VERIFICATION RUNNER ===")
    
    scripts = [
        "tests/fmu_validation/test_fmu_simulation.py",
        "tests/distance_homogeneity/test_distance_homogeneity.py",
        "tests/parameter_fitting/test_parameter_fitting.py",
        "tests/psd_analysis/plot_direct_fft.py",
        "tests/psd_analysis/plot_infinite_comparison.py",
        "tests/psd_analysis/plot_psd_comparison.py",
        "tests/psd_analysis/run_advanced_sensitivity.py",
        "tests/psd_analysis/run_advanced_analysis.py"
    ]
    
    success = True
    for script in scripts:
        if not run_script(script):
            success = False
            
    print("\n========================================================")
    if success:
        print("ALL TESTS COMPLETED SUCCESSFULLY!")
        print("Regenerated artifacts are located in their respective subdirectories:")
        print("  - tests/fmu_validation/ (verified concurrent FMI simulation & determinism)")
        print("  - tests/distance_homogeneity/distance_homogeneity_curves.png")
        print("  - tests/parameter_fitting/parameter_fitting_case_1.png")
        print("  - tests/parameter_fitting/parameter_fitting_case_2.png")
        print("  - tests/parameter_fitting/parameter_fitting_case_3.png")
        print("  - tests/parameter_fitting/parameter_fitting_summary.png")
        print("  - tests/parameter_fitting/parameter_fitting_analysis.md")
        print("  - tests/psd_analysis/psd_direct_fft_comparison.png")
        print("  - tests/psd_analysis/psd_infinite_comparison.png")
        print("  - tests/psd_analysis/psd_comparison_curves.png")
        print("  - tests/psd_analysis/psd_sensitivity_Nf_raw_cum.png")
        print("  - tests/psd_analysis/psd_sensitivity_Ntheta_raw_cum.png")
        print("  - tests/psd_analysis/psd_multi_params.png")
        print("  - tests/psd_analysis/psd_sensitivity_Nf.png")
        print("  - tests/psd_analysis/advanced_psd_analysis.md")
        sys.exit(0)
    else:
        print("SOME TESTS ENCOUNTERED ERRORS. Please check console output.")
        sys.exit(1)

if __name__ == "__main__":
    main()
