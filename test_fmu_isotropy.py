import os
import sys
import time
import numpy as np
import matplotlib.pyplot as plt

def get_ram_usage():
    try:
        import psutil
        process = psutil.Process(os.getpid())
        return f"{process.memory_info().rss / (1024**2):.1f} MB"
    except ImportError:
        return "N/A"

def get_I(alpha):
    t = np.linspace(-2000, 2000, 200000)
    dt = t[1] - t[0]
    return np.sum((1.0 + t**2)**(-alpha/2.0)) * dt

def custom_welch(y, fs, nperseg):
    """
    Custom implementation of Welch's periodogram method to estimate PSD,
    correcting for window energy loss.
    """
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
    def __init__(self, seed=42, road_class='C', w=2.0, f_min=0.002, f_max=20.0):
        self.seed = seed
        self.road_class = road_class
        self.w = w
        self.f_min = f_min
        self.f_max = f_max
        self._init_waves()

    def _init_waves(self):
        rng = np.random.RandomState(self.seed)
        
        class_map = {
            'A': 16e-6,
            'B': 64e-6,
            'C': 256e-6,
            'D': 1024e-6,
            'E': 4096e-6
        }
        Gd_n0 = class_map.get(self.road_class.upper(), 256e-6)
        n0 = 0.1
        C1 = Gd_n0 * (n0**self.w)
        
        alpha = self.w + 1.0
        I_val = get_I(alpha)
        
        # Corrected continuous scaling coefficient C2:
        # C2 = C1 / (4.0 * I_val)
        # We use 4.0 in the denominator because our sum-of-sinusoids model generates independent random
        # phases for all angles in [0, 2pi), meaning opposite directions (theta and theta + pi) are
        # uncorrelated, which doubles the projected variance compared to a conjugate symmetric FFT.
        C2 = C1 / (4.0 * I_val)
        
        # Grid sizes for wave components
        Nf = 250
        Ntheta = 16
        
        f_r = np.linspace(self.f_min, self.f_max, Nf + 1)
        df_r = np.diff(f_r)
        
        theta = np.linspace(0, 2*np.pi, Ntheta, endpoint=False)
        dtheta = 2*np.pi / Ntheta
        
        amps = []
        kx = []
        ky = []
        phis = []
        
        for i in range(Nf):
            # Randomize frequencies within linear bins (Shinozuka Method) to smooth out spectral lines
            fc = rng.uniform(f_r[i], f_r[i+1])
            dfc = df_r[i]
            
            S_2D_val = C2 * (fc**(-alpha))
            power_per_angle = S_2D_val * fc * dfc * dtheta
            amp = np.sqrt(2.0 * power_per_angle)
            
            for j in range(Ntheta):
                # Randomize angles within sectors to ensure continuous isotropic direction distribution
                th = theta[j] + rng.uniform(-dtheta/2, dtheta/2)
                phi = rng.uniform(0, 2*np.pi)
                
                amps.append(amp)
                kx.append(2.0 * np.pi * fc * np.cos(th))
                ky.append(2.0 * np.pi * fc * np.sin(th))
                phis.append(phi)
                
        self.amps = np.array(amps)
        self.kx = np.array(kx)
        self.ky = np.array(ky)
        self.phis = np.array(phis)

    def height(self, x, y):
        # Pure 1D loop: The most cache-friendly and fastest way to evaluate in NumPy
        x = np.atleast_1d(x)
        y = np.atleast_1d(y)
        h = np.zeros_like(x)
        for amp, kx, ky, phi in zip(self.amps, self.kx, self.ky, self.phis):
            h += amp * np.cos(kx * x + ky * y + phi)
        return h

def main():
    w_target = 2.0
    road_class = 'C'
    Gd_n0_target = 256e-6
    
    # Extended frequency cutoffs
    f_min = 0.002
    f_max = 20.0
    
    print("=== ISO 8608 FMU ISOTROPY & HOMOGENEITY VALIDATION ===", flush=True)
    print(f"Road Class: {road_class} (Gd(n0) = 2.56e-4 m^3)", flush=True)
    print(f"Cutoff Frequencies: f_min = {f_min} Hz, f_max = {f_max} Hz", flush=True)
    
    # Initialize road surface model
    road = SumOfSinusoidsRoad(seed=42, road_class=road_class, w=w_target, f_min=f_min, f_max=f_max)
    print(f"Initialized {len(road.amps)} wave components. (RAM: {get_ram_usage()})", flush=True)
    
    # Slice parameters
    num_slices = 100
    slice_length = 1000.0  # 1 km long slices
    dx_slice = 0.025       # spacing 0.025m (fs = 40 Hz, Nyquist = 20 Hz)
    N_slice = int(slice_length / dx_slice)
    fs = 1.0 / dx_slice
    
    # Welch settings: segment length of 4096 points (~100m) to smooth out wave spikes
    nperseg = 4096
    
    np.random.seed(54321)
    slice_paths = []
    all_psds = []
    
    # Welch frequencies (exclude DC)
    freqs_welch, _ = custom_welch(np.zeros(N_slice), fs, nperseg)
    freqs_welch = freqs_welch[1:]
    
    print(f"\nGenerating {num_slices} random slices of length {slice_length}m ({N_slice} points each)...", flush=True)
    
    start_time = time.time()
    for i in range(num_slices):
        x1 = np.random.uniform(-5000.0, 5000.0)
        y1 = np.random.uniform(-5000.0, 5000.0)
        theta = np.random.uniform(0, 2*np.pi)
        
        x2 = x1 + slice_length * np.cos(theta)
        y2 = y1 + slice_length * np.sin(theta)
        
        slice_paths.append(((x1, y1), (x2, y2)))
        
        s = np.linspace(0, slice_length, N_slice, endpoint=False)
        px = x1 + s * np.cos(theta)
        py = y1 + s * np.sin(theta)
        
        # Compute heights (continuous and 1D loop)
        z_slice = road.height(px, py)
        
        # Estimate PSD using Welch's method
        f_w, psd_w = custom_welch(z_slice, fs=fs, nperseg=nperseg)
        all_psds.append(psd_w[1:])
        
        # Progress logging & ETA
        elapsed_total = time.time() - start_time
        avg_time_per_slice = elapsed_total / (i + 1)
        remaining_slices = num_slices - (i + 1)
        eta = avg_time_per_slice * remaining_slices
        
        if (i + 1) % 10 == 0 or i == 0 or (i + 1) == num_slices:
            print(f"Slice {i+1:3d}/{num_slices:3d} | "
                  f"Avg Time: {avg_time_per_slice:.2f}s | "
                  f"Elapsed: {elapsed_total:.1f}s | "
                  f"ETA: {eta:.1f}s | "
                  f"RAM: {get_ram_usage()}", flush=True)

    all_psds = np.array(all_psds)
    avg_psd = np.mean(all_psds, axis=0)
    
    # Target continuous PSD
    C1 = Gd_n0_target * (0.1**w_target)
    psd_target = C1 * (freqs_welch**(-w_target))
    
    # FITTING INDIVIDUAL SLICES
    # Fit from f_fit_min to f_max.
    # We set f_fit_min = 0.01 cycles/m because the Welch segment length (nperseg * dx = 102.4m)
    # cannot resolve any frequencies below 1/102.4 = 0.0098 cycles/m. Fitting below this limit
    # introduces statistical leakage and bias from the Welch window.
    f_fit_min = 0.01
    fit_idx = (freqs_welch >= f_fit_min) & (freqs_welch <= f_max)
    
    B_fits = []
    Gd_n0_fits = []
    
    for i in range(num_slices):
        p = np.polyfit(np.log(freqs_welch[fit_idx]), np.log(all_psds[i, fit_idx]), 1)
        w_fit = -p[0]
        C1_fit = np.exp(p[1])
        Gd_n0_fit = C1_fit / (0.1**w_fit)
        
        B_fits.append(w_fit)
        Gd_n0_fits.append(Gd_n0_fit)
        
    B_fits = np.array(B_fits)
    Gd_n0_fits = np.array(Gd_n0_fits)
    
    # Fit of the average PSD
    p_avg = np.polyfit(np.log(freqs_welch[fit_idx]), np.log(avg_psd[fit_idx]), 1)
    w_avg = -p_avg[0]
    C1_avg = np.exp(p_avg[1])
    Gd_n0_avg = C1_avg / (0.1**w_avg)
    
    print("\n--- STATISTICAL RESULTS OVER 100 INFINITE SLICES ---", flush=True)
    print(f"Target Exponent (w):        {w_target:.4f}", flush=True)
    print(f"Fitted Exponent (w_fit):    {np.mean(B_fits):.4f} \u00b1 {np.std(B_fits):.4f}", flush=True)
    print(f"Target Gd(n0):             {Gd_n0_target:.4e} m^3", flush=True)
    print(f"Fitted Gd(n0_fit):         {np.mean(Gd_n0_fits):.4e} \u00b1 {np.std(Gd_n0_fits):.4e} m^3", flush=True)
    print("-" * 60, flush=True)
    print("Fit of the Average PSD (global estimation):")
    print(f"  Global Exponent w:        {w_avg:.4f} (Error: {np.abs(w_avg - w_target)/w_target*100:.2f}%)", flush=True)
    print(f"  Global Roughness Gd(n0):  {Gd_n0_avg:.4e} m^3 (Error: {np.abs(Gd_n0_avg - Gd_n0_target)/Gd_n0_target*100:.2f}%)", flush=True)
    
    # Save the statistics to a text summary
    summary_filename = "fmu_isotropy_results.txt"
    with open(summary_filename, 'w') as f_out:
        f_out.write("=== FMU CONTINUOUS ROAD MULTI-SLICE SPECTRAL VALIDATION ===\n")
        f_out.write(f"Target Parameters: Gd(n0) = {Gd_n0_target:.4e} m^3, w = {w_target:.3f}\n")
        f_out.write(f"Extended Cutoffs: f_min = {f_min} Hz, f_max = {f_max} Hz\n")
        f_out.write(f"Welch Resolution Limit (Fit Min): {f_fit_min} Hz\n")
        f_out.write(f"Number of Slices: {num_slices}, Length: {slice_length}m, dx_slice: {dx_slice}m\n\n")
        f_out.write("Spectral Fitting Results of Individual Slices:\n")
        f_out.write(f"  Mean w: {np.mean(B_fits):.4f} (std: {np.std(B_fits):.4f})\n")
        f_out.write(f"  Mean Gd(n0): {np.mean(Gd_n0_fits):.4e} (std: {np.std(Gd_n0_fits):.4e} m^3)\n\n")
        f_out.write("Fit of the Average PSD (global estimation):\n")
        f_out.write(f"  Global Exponent w: {w_avg:.4f}\n")
        f_out.write(f"  Global Roughness Gd(n0): {Gd_n0_avg:.4e} m^3\n")
    print(f"\nSaved statistics summary to {summary_filename}", flush=True)
    
    # Plotting
    print("Generating validation plots...", flush=True)
    fig = plt.figure(figsize=(15, 10))
    
    # Plot 1: Slices on continuous coordinate space
    ax1 = fig.add_subplot(2, 2, 1)
    for idx, path in enumerate(slice_paths[:15]):
        p1, p2 = path
        ax1.plot([p1[0], p2[0]], [p1[1], p2[1]], label=f'Slice {idx+1}' if idx < 5 else '')
    ax1.set_title("Arbitrary Slice Paths in Infinite 2D Plane")
    ax1.set_xlabel("X (m)")
    ax1.set_ylabel("Y (m)")
    ax1.grid(True, linestyle='--', alpha=0.5)
    ax1.legend(loc='upper right')
    
    # Plot 2: PSD Comparison (Log-Log) over the extended band
    ax2 = fig.add_subplot(2, 2, 2)
    for psd in all_psds[:10]:
        ax2.loglog(freqs_welch, psd, color='gray', alpha=0.3)
    ax2.loglog(freqs_welch, avg_psd, 'blue', linewidth=2.0, label="Estimated Slice PSD (Average)")
    ax2.loglog(freqs_welch, psd_target, 'k--', linewidth=2.0, label="ISO 8608 Target PSD")
    ax2.axvline(f_fit_min, color='green', linestyle=':', label='Fitting Limits (Welch Resolution)')
    ax2.axvline(f_max, color='green', linestyle=':')
    ax2.set_title("Extended PSD Comparison (Log-Log)")
    ax2.set_xlabel("Spatial Frequency (cycles/m)")
    ax2.set_ylabel("PSD ($m^3$)")
    ax2.legend()
    ax2.grid(True, which="both", linestyle='--', alpha=0.5)
    
    # Plot 3: Fitted w vs Angle (Isotropy)
    ax3 = fig.add_subplot(2, 2, 3)
    angles = []
    for path in slice_paths:
        p1, p2 = path
        angle = np.arctan2(p2[1] - p1[1], p2[0] - p1[0]) * 180.0 / np.pi
        if angle < 0:
            angle += 180.0
        angles.append(angle)
        
    ax3.scatter(angles, B_fits, color='blue', alpha=0.7)
    ax3.axhline(w_target, color='black', linestyle='--', label=f'Target w ({w_target})')
    ax3.set_title("Exponent w vs. Slice Angle (Isotropy)")
    ax3.set_xlabel("Slice Angle (degrees)")
    ax3.set_ylabel("Fitted Exponent w")
    ax3.set_ylim(w_target - 0.4, w_target + 0.4)
    ax3.legend()
    ax3.grid(True, linestyle='--', alpha=0.5)
    
    # Plot 4: Fitted w vs Spatial Offset (Homogeneity)
    ax4 = fig.add_subplot(2, 2, 4)
    offsets = [np.sqrt(p[0][0]**2 + p[0][1]**2) for p in slice_paths]
    ax4.scatter(offsets, B_fits, color='blue', alpha=0.7)
    ax4.axhline(w_target, color='black', linestyle='--')
    ax4.set_title("Exponent w vs. Distance from Origin (Homogeneity)")
    ax4.set_xlabel("Starting Point Distance (m)")
    ax4.set_ylabel("Fitted Exponent w")
    ax4.set_ylim(w_target - 0.4, w_target + 0.4)
    ax4.grid(True, linestyle='--', alpha=0.5)
    
    plt.tight_layout()
    output_filename = "fmu_isotropy_validation.png"
    plt.savefig(output_filename, dpi=150)
    print(f"\nPlots saved to {output_filename}", flush=True)
    
    # Copy to artifact folder
    artifact_dir = r"C:\Users\novo\.gemini\antigravity\brain\5a4d8382-93f9-4de6-ba07-708b80ce6613"
    if os.path.exists(artifact_dir):
        dest_path = os.path.join(artifact_dir, "fmu_isotropy_validation.png")
        plt.savefig(dest_path, dpi=150)
        print(f"Plots copied to artifact dir: {dest_path}", flush=True)
        
        # Copy summary text to artifact folder too
        with open(os.path.join(artifact_dir, "fmu_isotropy_results.txt"), 'w') as f_out:
            f_out.write("=== FMU CONTINUOUS ROAD MULTI-SLICE SPECTRAL VALIDATION ===\n")
            f_out.write(f"Target Parameters: Gd(n0) = {Gd_n0_target:.4e} m^3, w = {w_target:.3f}\n")
            f_out.write(f"Extended Cutoffs: f_min = {f_min} Hz, f_max = {f_max} Hz\n")
            f_out.write(f"Welch Resolution Limit (Fit Min): {f_fit_min} Hz\n")
            f_out.write(f"Number of Slices: {num_slices}, Length: {slice_length}m, dx_slice: {dx_slice}m\n\n")
            f_out.write("Spectral Fitting Results of Individual Slices:\n")
            f_out.write(f"  Mean w: {np.mean(B_fits):.4f} (std: {np.std(B_fits):.4f})\n")
            f_out.write(f"  Mean Gd(n0): {np.mean(Gd_n0_fits):.4e} (std: {np.std(Gd_n0_fits):.4e} m^3)\n\n")
            f_out.write("Fit of the Average PSD (global estimation):\n")
            f_out.write(f"  Global Exponent w: {w_avg:.4f}\n")
            f_out.write(f"  Global Roughness Gd(n0): {Gd_n0_avg:.4e} m^3\n")

if __name__ == "__main__":
    main()
