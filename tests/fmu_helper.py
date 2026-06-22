import os
import shutil
import atexit
import numpy as np
from fmpy import read_model_description, extract
from fmpy.fmi2 import FMU2Slave

class FMURoadQuery:
    def __init__(self, fmu_path=None):
        if fmu_path is None:
            # find fmu_path relative to this helper file
            fmu_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "InfiniteRoadFMU.fmu"))
        self.fmu_path = fmu_path
        self.model_description = read_model_description(self.fmu_path)
        self.unzipdir = extract(self.fmu_path)
        self.guid = self.model_description.guid
        self.model_identifier = self.model_description.coSimulation.modelIdentifier
        self.var_refs = {v.name: v.valueReference for v in self.model_description.modelVariables}
        atexit.register(self.cleanup)

    def cleanup(self):
        if hasattr(self, 'unzipdir') and os.path.exists(self.unzipdir):
            shutil.rmtree(self.unzipdir, ignore_errors=True)

    def get_slave(self, seed=42, Gd_n0=256e-6, w=2.0, f_min=0.005, f_max=100.0, Nf=512, Ntheta=32, road_class=0, instance_name="fmu_instance"):
        slave = FMU2Slave(
            guid=self.guid,
            unzipDirectory=self.unzipdir,
            modelIdentifier=self.model_identifier,
            instanceName=instance_name
        )
        slave.instantiate()
        slave.setupExperiment(startTime=0.0)
        slave.enterInitializationMode()
        
        # Set parameters
        slave.setInteger([self.var_refs['seed']], [int(seed)])
        slave.setInteger([self.var_refs['road_class']], [int(road_class)])
        slave.setReal([self.var_refs['Gd_n0']], [float(Gd_n0)])
        slave.setReal([self.var_refs['w']], [float(w)])
        slave.setReal([self.var_refs['f_min']], [float(f_min)])
        slave.setReal([self.var_refs['f_max']], [float(f_max)])
        
        # Set Nf and Ntheta if exposed as variables
        if 'Nf' in self.var_refs:
            slave.setInteger([self.var_refs['Nf']], [int(Nf)])
        if 'Ntheta' in self.var_refs:
            slave.setInteger([self.var_refs['Ntheta']], [int(Ntheta)])
            
        slave.exitInitializationMode()
        return slave

    def query_profile(self, slave, px, py):
        """Query coordinates px, py (can be float or numpy array)."""
        x_ref = self.var_refs['x']
        y_ref = self.var_refs['y']
        z_ref = self.var_refs['z']
        
        if isinstance(px, (int, float)) or np.isscalar(px):
            slave.setReal([x_ref, y_ref], [float(px), float(py)])
            return slave.getReal([z_ref])[0]
        else:
            # array input
            z = np.zeros_like(px)
            for i in range(len(px)):
                slave.setReal([x_ref, y_ref], [float(px[i]), float(py[i])])
                z[i] = slave.getReal([z_ref])[0]
            return z

    def query_profile_parallel(self, px, py, num_threads=8, seed=42, Gd_n0=256e-6, w=2.0, f_min=0.005, f_max=100.0, Nf=512, Ntheta=32, road_class=0):
        """Query coordinates px, py in parallel using multiple FMI slave instances."""
        from concurrent.futures import ThreadPoolExecutor
        
        n_points = len(px)
        chunk_size = int(np.ceil(n_points / num_threads))
        z = np.zeros(n_points)
        
        def worker(thread_idx):
            start_idx = thread_idx * chunk_size
            end_idx = min(start_idx + chunk_size, n_points)
            if start_idx >= end_idx:
                return
                
            slave = self.get_slave(
                seed=seed, Gd_n0=Gd_n0, w=w, f_min=f_min, f_max=f_max,
                Nf=Nf, Ntheta=Ntheta, road_class=road_class,
                instance_name=f"parallel_slave_{thread_idx}"
            )
            
            x_ref = self.var_refs['x']
            y_ref = self.var_refs['y']
            z_ref = self.var_refs['z']
            
            x_vals = px[start_idx:end_idx]
            y_vals = py[start_idx:end_idx]
            z_vals = np.zeros(len(x_vals))
            
            for i in range(len(x_vals)):
                slave.setReal([x_ref, y_ref], [float(x_vals[i]), float(y_vals[i])])
                z_vals[i] = slave.getReal([z_ref])[0]
                
            z[start_idx:end_idx] = z_vals
            
            slave.terminate()
            slave.freeInstance()
            
        with ThreadPoolExecutor(max_workers=num_threads) as executor:
            list(executor.map(worker, range(num_threads)))
            
        return z
