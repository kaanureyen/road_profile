import numpy as np
from pythonfmu import Fmi2Causality, Fmi2Slave, Fmi2Variability, Real, Integer

class InfiniteRoadFMU(Fmi2Slave):
    author = "Antigravity Coding Assistant"
    description = "Deterministic ISO 8608 2D Isotropic Infinite Road Profile Generator (Single Point)"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        
        # Expose inputs as standard class attributes
        self.x = 0.0
        self.y = 0.0
        
        # Expose parameters as standard class attributes
        self.seed = 42
        self.road_class = 3  # 1=A, 2=B, 3=C, 4=D, 5=E, 0=Custom Gd_n0
        self.Gd_n0 = 256e-6
        self.w = 2.0
        self.f_min = 0.01
        self.f_max = 10.0
        self.Nf = 512
        self.Ntheta = 32
        
        # Register variables for FMI interface
        self.register_variable(Real("x", causality=Fmi2Causality.input))
        self.register_variable(Real("y", causality=Fmi2Causality.input))
        
        # Output defined as property to dynamically query current input coords
        self.register_variable(Real("z", causality=Fmi2Causality.output))
        
        self.register_variable(Integer("seed", causality=Fmi2Causality.parameter, variability=Fmi2Variability.fixed))
        self.register_variable(Integer("road_class", causality=Fmi2Causality.parameter, variability=Fmi2Variability.fixed))
        self.register_variable(Real("Gd_n0", causality=Fmi2Causality.parameter, variability=Fmi2Variability.fixed))
        self.register_variable(Real("w", causality=Fmi2Causality.parameter, variability=Fmi2Variability.fixed))
        self.register_variable(Real("f_min", causality=Fmi2Causality.parameter, variability=Fmi2Variability.fixed))
        self.register_variable(Real("f_max", causality=Fmi2Causality.parameter, variability=Fmi2Variability.fixed))
        self.register_variable(Integer("Nf", causality=Fmi2Causality.parameter, variability=Fmi2Variability.fixed))
        self.register_variable(Integer("Ntheta", causality=Fmi2Causality.parameter, variability=Fmi2Variability.fixed))
        
        # Internal wave representation cache
        self._last_params = None
        self._amps = None
        self._kx = None
        self._ky = None
        self._phis = None

    def _get_I(self, alpha):
        """Numerical integration of (1 + t^2)^(-alpha/2)."""
        t = np.linspace(-2000, 2000, 200000)
        dt = t[1] - t[0]
        return np.sum((1.0 + t**2)**(-alpha/2.0)) * dt

    def _lazy_init(self):
        current_params = (self.seed, self.road_class, self.Gd_n0, self.w, self.f_min, self.f_max, self.Nf, self.Ntheta)
        if self._amps is not None and self._last_params == current_params:
            return
            
        rng = np.random.RandomState(self.seed)
        
        # Map road class to Gd(n0)
        class_map = {
            1: 16e-6,   # Class A
            2: 64e-6,   # Class B
            3: 256e-6,  # Class C
            4: 1024e-6, # Class D
            5: 4096e-6  # Class E
        }
        Gd_n0_val = class_map.get(self.road_class, self.Gd_n0)
        n0 = 0.1
        C1 = Gd_n0_val * (n0**self.w)
        
        alpha = self.w + 1.0
        I_val = self._get_I(alpha)
        
        # Corrected continuous scaling coefficient C2:
        # C2 = C1 / (2.0 * I_val)
        # We use 2.0 in the denominator because the 1D slice projection of the 2D wave field
        # yields S_1D(f) = 2.0 * C2 * I_val * f^-w.
        C2 = C1 / (2.0 * I_val)
        
        Nf = self.Nf
        Ntheta = self.Ntheta

        
        # Logarithmic frequency spacing to capture low-frequency energy accurately
        f_r = np.logspace(np.log10(self.f_min), np.log10(self.f_max), Nf + 1)
        df_r = np.diff(f_r)
        f_centers = 0.5 * (f_r[:-1] + f_r[1:])
        
        # Angular grid
        theta = np.linspace(0, 2*np.pi, Ntheta, endpoint=False)
        dtheta = 2*np.pi / Ntheta
        
        amps = []
        kx = []
        ky = []
        phis = []
        
        for i in range(Nf):
            fc = f_centers[i]
            dfc = df_r[i]
            S_2D_val = C2 * (fc**(-alpha))
            power_per_angle = S_2D_val * fc * dfc * dtheta
            amp = np.sqrt(2.0 * power_per_angle)
            
            for j in range(Ntheta):
                th = theta[j]
                phi = rng.uniform(0, 2*np.pi)
                
                amps.append(amp)
                kx.append(2.0 * np.pi * fc * np.cos(th))
                ky.append(2.0 * np.pi * fc * np.sin(th))
                phis.append(phi)
                
        self._amps = np.array(amps)
        self._kx = np.array(kx)
        self._ky = np.array(ky)
        self._phis = np.array(phis)
        self._last_params = current_params

    def _compute_height(self, px, py):
        self._lazy_init()
        # Fast vectorized sum of cosine wave components
        phases = self._kx * px + self._ky * py + self._phis
        return float(np.sum(self._amps * np.cos(phases)))

    # Property for output to dynamically query coordinates
    @property
    def z(self):
        return self._compute_height(self.x, self.y)
    
    @z.setter
    def z(self, val):
        pass

    def do_step(self, current_time, step_size):
        # Outputs are dynamically computed via properties when queried,
        # but returning True is required by FMI standard to indicate success.
        return True
