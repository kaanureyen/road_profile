import os
import shutil
import numpy as np
from fmpy import read_model_description, extract
from fmpy.fmi2 import FMU2Slave

def main():
    fmu_filename = "InfiniteRoadFMU.fmu"
    
    if not os.path.exists(fmu_filename):
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
