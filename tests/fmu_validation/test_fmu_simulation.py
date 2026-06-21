import os
import shutil
import numpy as np
from fmpy import read_model_description, extract
from fmpy.fmi2 import FMU2Slave

def main():
    fmu_filename = "InfiniteRoadFMU.fmu"
    if not os.path.exists(fmu_filename):
        candidate = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", fmu_filename))
        if os.path.exists(candidate):
            fmu_filename = candidate
        else:
            print(f"Error: {fmu_filename} not found!")
            return

    print("=== READING FMU MODEL DESCRIPTION ===")
    model_description = read_model_description(fmu_filename)
    print(f"Model Name:        {model_description.modelName}")
    print(f"FMI Version:       {model_description.fmiVersion}")
    print(f"GUID:              {model_description.guid}")
    print(f"Description:       {model_description.description}")
    
    print("\nModel Variables:")
    for v in model_description.modelVariables:
        print(f"  {v.name:<12} | causality: {str(v.causality):<10} | type: {v.type:<8} | description: {v.description}")
        
    print("\n=== EXTRACTING AND INSTANTIATING 4 CONCURRENT FMUs ===")
    # Extract the FMU archive to a temporary directory
    unzipdir = extract(fmu_filename)
    
    slaves = []
    try:
        # We will instantiate 4 slaves concurrently to simulate 4 wheel contacts
        for i in range(4):
            slave = FMU2Slave(
                guid=model_description.guid,
                unzipDirectory=unzipdir,
                modelIdentifier=model_description.coSimulation.modelIdentifier,
                instanceName=f"road_wheel_{i+1}"
            )
            slave.instantiate()
            slave.setupExperiment(startTime=0.0)
            slave.enterInitializationMode()
            
            # Set the seed for all slaves to 42 (matching seed)
            var_refs = {v.name: v.valueReference for v in model_description.modelVariables}
            slave.setInteger([var_refs['seed']], [42])
            
            slave.exitInitializationMode()
            slaves.append(slave)
            
        print("Successfully instantiated and initialized 4 independent FMU instances!")
        
        # Helper dictionary for value references
        var_refs = {v.name: v.valueReference for v in model_description.modelVariables}
        x_ref = var_refs['x']
        y_ref = var_refs['y']
        z_ref = var_refs['z']
        
        # TEST 1: Simultaneous different coordinate queries (4 wheels of a car)
        print("\n--- TEST 1: Querying 4 different wheel coordinates simultaneously ---")
        # Define 4 wheel positions:
        # Wheel 1 (Front-Left):  (10.0, 5.0)
        # Wheel 2 (Front-Right): (10.0, 3.5)
        # Wheel 3 (Rear-Left):   (7.0, 5.0)
        # Wheel 4 (Rear-Right):  (7.0, 3.5)
        wheel_coords = [
            (10.0, 5.0),
            (10.0, 3.5),
            (7.0, 5.0),
            (7.0, 3.5)
        ]
        
        heights = []
        for i, (wx, wy) in enumerate(wheel_coords):
            slave = slaves[i]
            slave.setReal([x_ref, y_ref], [wx, wy])
            # Retrieve height
            wz = slave.getReal([z_ref])[0]
            heights.append(wz)
            print(f"Wheel {i+1} at ({wx:5.1f}, {wy:5.1f}) returned height z = {wz:+.6f} m")
            
        # TEST 2: Determinism check across separate instances
        print("\n--- TEST 2: Determinism & Repeatability (Cross-Instance Query) ---")
        # Let's query Wheel 1's coordinate (10.0, 5.0) on all 4 independent instances
        print("Querying coordinate (10.0, 5.0) on all 4 instances:")
        test_heights = []
        for i, slave in enumerate(slaves):
            slave.setReal([x_ref, y_ref], [10.0, 5.0])
            wz = slave.getReal([z_ref])[0]
            test_heights.append(wz)
            print(f"  Instance {i+1} returned z = {wz:+.12f} m")
            
        # Verify that all instances returned the exact same value
        first_val = test_heights[0]
        for val in test_heights:
            assert np.isclose(first_val, val, atol=1e-15), "Determinism failed! Different instances returned different heights for the same coordinate."
        print(f"SUCCESS: Cross-instance repeatability verified! Max diff = {np.max(np.abs(np.array(test_heights) - first_val)):.2e} m")
        
        # TEST 3: Seed sensitivity
        print("\n--- TEST 3: Seed Sensitivity (Realization Check) ---")
        # Instantiate another slave but with seed = 99
        slave_seed99 = FMU2Slave(
            guid=model_description.guid,
            unzipDirectory=unzipdir,
            modelIdentifier=model_description.coSimulation.modelIdentifier,
            instanceName="road_wheel_seed99"
        )
        slave_seed99.instantiate()
        slave_seed99.setupExperiment(startTime=0.0)
        slave_seed99.enterInitializationMode()
        slave_seed99.setInteger([var_refs['seed']], [99])
        slave_seed99.exitInitializationMode()
        
        # Query (10.0, 5.0) on seed 99 instance
        slave_seed99.setReal([x_ref, y_ref], [10.0, 5.0])
        z_seed99 = slave_seed99.getReal([z_ref])[0]
        print(f"Coordinate (10.0, 5.0) on Seed 42 instance: z = {first_val:+.12f} m")
        print(f"Coordinate (10.0, 5.0) on Seed 99 instance: z = {z_seed99:+.12f} m")
        
        diff = first_val - z_seed99
        print(f"Difference (Seed 42 vs 99): {diff:+.6f} m")
        assert abs(diff) > 1e-4, "Seed parameter did not change the road realization!"
        print("SUCCESS: Different seeds generate independent random realizations!")
        
        # Generate plot for validation report
        import matplotlib.pyplot as plt
        print("\n--- Generating validation plots for the report ---")
        s_eval = np.linspace(0.0, 20.0, 500)
        
        # We will reuse the first wheel instance (seed 42) and the seed 99 instance
        z_wheel1 = []
        z_wheel2 = []
        z_seed99_arr = []
        
        for sx in s_eval:
            slaves[0].setReal([x_ref, y_ref], [sx, 5.0])
            z_wheel1.append(slaves[0].getReal([z_ref])[0])
            
            slaves[1].setReal([x_ref, y_ref], [sx, 3.5])
            z_wheel2.append(slaves[1].getReal([z_ref])[0])
            
            slave_seed99.setReal([x_ref, y_ref], [sx, 5.0])
            z_seed99_arr.append(slave_seed99.getReal([z_ref])[0])
            
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
        
        ax1.plot(s_eval, z_wheel1, label='Wheel 1 (y = 5.0 m)', color='#1f77b4')
        ax1.plot(s_eval, z_wheel2, label='Wheel 2 (y = 3.5 m)', color='#ff7f0e', linestyle='--')
        ax1.set_title('Concurrent Querying on Distinct Spatial Paths (Seed 42)')
        ax1.set_xlabel('Longitudinal Position x (m)')
        ax1.set_ylabel('Road Height z (m)')
        ax1.grid(True, linestyle='--', alpha=0.6)
        ax1.legend()
        
        ax2.plot(s_eval, z_wheel1, label='Seed 42 (Realization 1)', color='#1f77b4')
        ax2.plot(s_eval, z_seed99_arr, label='Seed 99 (Realization 2)', color='#d62728', linestyle=':')
        ax2.set_title('Seed Sensitivity: Independent Realizations at y = 5.0 m')
        ax2.set_xlabel('Longitudinal Position x (m)')
        ax2.set_ylabel('Road Height z (m)')
        ax2.grid(True, linestyle='--', alpha=0.6)
        ax2.legend()
        
        plt.tight_layout()
        
        script_dir = os.path.dirname(os.path.abspath(__file__))
        plot_path = os.path.join(script_dir, "fmu_simulation_results.png")
        plt.savefig(plot_path, dpi=150)
        print(f"Validation plot saved to {plot_path}")
        
        # Copy to artifact folder if available
        artifact_dir = os.environ.get("ANTIGRAVITY_ARTIFACT_DIR")
        if not artifact_dir:
            artifact_dir = r"C:\Users\novo\.gemini\antigravity\brain\eb2516c5-ab48-42c6-b6d8-0b90cc4ca6ca"
        try:
            os.makedirs(artifact_dir, exist_ok=True)
            plt.savefig(os.path.join(artifact_dir, "fmu_simulation_results.png"), dpi=150)
        except Exception as e:
            print(f"Could not save copy to artifact: {e}")
        plt.close()
        
        slave_seed99.terminate()
        slave_seed99.freeInstance()

        
    finally:
        # Clean up slaves
        for slave in slaves:
            try:
                slave.terminate()
                slave.freeInstance()
            except Exception:
                pass
        # Clean up unzip directory
        shutil.rmtree(unzipdir, ignore_errors=True)

if __name__ == "__main__":
    main()
