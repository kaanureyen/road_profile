import os
import sys
import subprocess
import time

def run_script(script_name):
    print(f"\n========================================================")
    print(f"Running {script_name}...")
    print(f"========================================================")
    t0 = time.time()
    
    # Run with unbuffered output so logs print live
    res = subprocess.run([sys.executable, "-u", script_name], capture_output=False)
    
    elapsed = time.time() - t0
    if res.returncode == 0:
        print(f"SUCCESS: {script_name} finished in {elapsed:.1f}s")
        return True
    else:
        print(f"FAILED: {script_name} exited with code {res.returncode} in {elapsed:.1f}s")
        return False

def main():
    print("=== ROAD PROFILE GENERATOR VERIFICATION RUNNER ===")
    
    # Ensure tests directory exists
    os.makedirs("tests", exist_ok=True)
    
    scripts = [
        "test_distance_homogeneity.py",
        "test_parameter_fitting.py"
    ]
    
    success = True
    for script in scripts:
        if not run_script(script):
            success = False
            
    print("\n========================================================")
    if success:
        print("ALL TESTS COMPLETED SUCCESSFULLY!")
        print("Regenerated artifacts are located in the 'tests/' directory:")
        print("  - tests/distance_homogeneity_curves.png")
        print("  - tests/parameter_fitting_case_1.png")
        print("  - tests/parameter_fitting_case_2.png")
        print("  - tests/parameter_fitting_case_3.png")
        print("  - tests/parameter_fitting_summary.png")
        print("  - tests/parameter_fitting_analysis.md")
        sys.exit(0)
    else:
        print("SOME TESTS ENCOUNTERED ERRORS. Please check console output.")
        sys.exit(1)

if __name__ == "__main__":
    main()
