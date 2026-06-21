import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate

def get_I(alpha):
    t = np.linspace(-2000, 2000, 200000)
    dt = t[1] - t[0]
    return np.sum((1.0 + t**2)**(-alpha/2.0)) * dt

def custom_welch(y, fs, nperseg):
    win = np.hanning(nperseg)
    win_norm = np.sum(win**2)
    step = nperseg // 2
    num_segments = (len(y) - nperseg) // step + 1
    freqs = np.fft.rfftfreq(nperseg, d=1.0/fs)
    psd_accum = np.zeros(len(freqs))
    for i in range(num_segments):
        start = i * step
        end = start + nperseg
        seg = y[start:end]
        seg_detrended = seg - np.mean(seg)
        seg_windowed = seg_detrended * win
        fft_seg = np.fft.rfft(seg_windowed)
        psd_accum += np.abs(fft_seg)**2
    psd = (2.0 / (fs * win_norm)) * (psd_accum / num_segments)
    return freqs, psd

class SumOfSinusoidsRoad:
    def __init__(self, Gd_n0=64e-6, w=2.0, f_min=0.002, f_max=20.0, Nf=64, Ntheta=16):
        self.Gd_n0 = Gd_n0
        self.w = w
        self.f_min = f_min
        self.f_max = f_max
        self.Nf = Nf
        self.Ntheta = Ntheta
        self._init_waves()

    def _init_waves(self):
        rng = np.random.RandomState(42)
        n0 = 0.1
        C1 = self.Gd_n0 * (n0**self.w)
        alpha = self.w + 1.0
        I_val = get_I(alpha)
        self.C2 = C1 / (2.0 * I_val)
        
        f_r = np.logspace(np.log10(self.f_min), np.log10(self.f_max), self.Nf + 1)
        df_r = np.diff(f_r)
        f_centers = 0.5 * (f_r[:-1] + f_r[1:])
        
        theta = np.linspace(0, 2*np.pi, self.Ntheta, endpoint=False)
        dtheta = 2*np.pi / self.Ntheta
        
        self.amps = []
        self.kx = []
        self.ky = []
        self.phis = []
        
        for i in range(self.Nf):
            fc = f_centers[i]
            dfc = df_r[i]
            S_2D_val = self.C2 * (fc**(-alpha))
            power_per_angle = S_2D_val * fc * dfc * dtheta
            amp = np.sqrt(2.0 * power_per_angle)
            for j in range(self.Ntheta):
                th = theta[j]
                phi = rng.uniform(0, 2*np.pi)
                self.amps.append(amp)
                self.kx.append(2.0 * np.pi * fc * np.cos(th))
                self.ky.append(2.0 * np.pi * fc * np.sin(th))
                self.phis.append(phi)
                
        self.amps = np.array(self.amps)
        self.kx = np.array(self.kx)
        self.ky = np.array(self.ky)
        self.phis = np.array(self.phis)

    def height(self, x, y):
        h = np.zeros_like(x)
        for amp, kx, ky, phi in zip(self.amps, self.kx, self.ky, self.phis):
            h += amp * np.cos(kx * x + ky * y + phi)
        return h

# Exact isotropic cumulative projection model
def exact_isotropic_cum_model(f_array, C1, w):
    alpha = w + 1.0
    I_val = get_I(alpha)
    results = []
    for f_val in f_array:
        val, _ = integrate.quad(lambda f_2D: (f_2D**(-w)) * np.arccos(f_val / f_2D), f_val, 20.0)
        results.append(C1 * (2.0 / I_val) * val)
    return np.array(results)

def main():
    G_target = 64e-6
    w_target = 2.0
    
    road = SumOfSinusoidsRoad(Gd_n0=G_target, w=w_target)
    
    # 10,000 meters slice to resolve frequencies down to 0.0001 cycles/m
    slice_length = 10000.0
    dx = 0.025
    N_slice = int(slice_length / dx)
    fs = 1.0 / dx
    
    print("Generating 10,000m road slice...", flush=True)
    np.random.seed(42)
    x1 = np.random.uniform(-5000.0, 5000.0)
    y1 = np.random.uniform(-5000.0, 5000.0)
    theta = np.random.uniform(0, 2*np.pi)
    
    s = np.linspace(0, slice_length, N_slice, endpoint=False)
    px = x1 + s * np.cos(theta)
    py = y1 + s * np.sin(theta)
    z = road.height(px, py)
    
    print("Computing Welch PSD with nperseg=32768...", flush=True)
    nperseg = 32768 # High resolution window
    freqs, psd = custom_welch(z, fs=fs, nperseg=nperseg)
    freqs = freqs[1:]
    psd = psd[1:]
    df = freqs[1] - freqs[0]
    
    cum_psd = np.cumsum(psd[::-1])[::-1] * df
    
    # Theoretical infinite cumulative curve
    C1_target = G_target * (0.1**w_target)
    infinite_cum = (C1_target / (w_target - 1.0)) * (freqs**(-(w_target - 1.0)))
    
    # Exact cumulative PSD model (incorporating 20.0 cycles/m truncation)
    print("Computing exact cumulative model...", flush=True)
    exact_cum = exact_isotropic_cum_model(freqs, C1_target, w_target)
    
    # Plotting
    plt.figure(figsize=(10, 7))
    
    # Welch Cumulative PSD
    plt.loglog(freqs, cum_psd, color='blue', linewidth=2.5, label='Welch Cumulative PSD (10,000m slice)')
    
    # Infinite Theoretical Cumulative curve
    plt.loglog(freqs, infinite_cum, color='red', linestyle='--', linewidth=2, label='Infinite Theoretical Model (No Cutoff)')
    
    # Exact Cumulative Model (Truncated at 20.0 cycles/m)
    plt.loglog(freqs, exact_cum, color='green', linestyle=':', linewidth=2, label='Exact Model (Truncated at 20.0 cycles/m)')
    
    # Plot limits and markers
    plt.axvline(0.002, color='orange', linestyle='--', alpha=0.7, label='Lower Cutoff (f_min = 0.002)')
    plt.axvline(0.0002, color='purple', linestyle=':', alpha=0.7, label='1/10th of Lower Cutoff (0.0002)')
    
    plt.xlim(0.0001, 20.0)
    plt.title("Cumulative PSD: Infinite vs. Truncated Isotropic Model")
    plt.xlabel("Spatial Frequency (cycles/m)")
    plt.ylabel("Cumulative Power (m²)")
    plt.grid(True, which="both", linestyle='--', alpha=0.5)
    plt.legend()
    plt.tight_layout()
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_path = os.path.join(script_dir, "psd_infinite_comparison.png")
    plt.savefig(output_path, dpi=150)
    plt.close()
    print(f"Plot saved successfully to: {output_path}", flush=True)
    
    artifact_dir = os.environ.get("ANTIGRAVITY_ARTIFACT_DIR")
    if not artifact_dir:
        artifact_dir = r"C:\Users\novo\.gemini\antigravity\brain\6bd8f97d-a5dd-4779-9267-df20885b87f3"
    try:
        os.makedirs(artifact_dir, exist_ok=True)
        plt.savefig(os.path.join(artifact_dir, "psd_infinite_comparison.png"), dpi=150)
    except Exception as e:
        print(f"Could not save to artifact directory: {e}")

if __name__ == "__main__":
    main()
