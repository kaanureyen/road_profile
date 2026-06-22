import os
import sys
import warnings
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from scipy.integrate import IntegrationWarning

warnings.filterwarnings("ignore", category=IntegrationWarning)

# Add tests/ to path to import fmu_helper
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from fmu_helper import FMURoadQuery

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

# Exact isotropic cumulative projection model
def exact_isotropic_cum_model(f_array, C1, w, f_min=0.01, f_max=2.0):
    alpha = w + 1.0
    I_val = get_I(alpha)
    results = []
    for f_val in f_array:
        if f_val >= f_max:
            results.append(0.0)
        else:
            f_start = max(f_val, f_min)
            u_start = np.arccos(np.clip(f_val / f_start, -1.0, 1.0))
            u_end = np.arccos(np.clip(f_val / f_max, -1.0, 1.0))
            val, _ = integrate.quad(
                lambda u: u * np.sin(u) * (np.cos(u)**(w - 2.0)),
                u_start,
                u_end
            )
            results.append(C1 * (2.0 / I_val) * (f_val**(1.0 - w)) * val)
    return np.array(results)

def main():
    G_target = 64e-6
    w_target = 2.0
    
    # Initialize FMU query helper
    fmu_query = FMURoadQuery()
    slave = fmu_query.get_slave(Gd_n0=G_target, w=w_target, f_min=0.01, f_max=2.0, Nf=512, Ntheta=32)
    
    slice_length = 5000.0
    dx = 0.25
    N_slice = int(slice_length / dx)
    fs = 1.0 / dx
    nperseg = 400
    df = fs / nperseg
    
    # Generate 1 slice
    np.random.seed(42)
    x1 = np.random.uniform(-5000.0, 5000.0)
    y1 = np.random.uniform(-5000.0, 5000.0)
    theta = np.random.uniform(0, 2*np.pi)
    
    s = np.linspace(0, slice_length, N_slice, endpoint=False)
    px = x1 + s * np.cos(theta)
    py = y1 + s * np.sin(theta)
    
    # Query using FMU
    z = fmu_query.query_profile(slave, px, py)
    slave.terminate()
    slave.freeInstance()
    
    # 1. Welch PSD for raw plot (smooth representation)
    freqs, psd = custom_welch(z, fs=fs, nperseg=nperseg)
    freqs = freqs[1:]
    psd = psd[1:]
    
    # 2. Direct FFT PSD for exact cumulative PSD representation
    win_full = np.hanning(len(z))
    win_full_norm = np.sum(win_full**2)
    z_det = z - np.mean(z)
    z_win = z_det * win_full
    fft_full = np.fft.rfft(z_win)
    freqs_fft = np.fft.rfftfreq(len(z), d=dx)
    psd_fft = (2.0 / (fs * win_full_norm)) * (np.abs(fft_full)**2)
    
    freqs_fft = freqs_fft[1:]
    psd_fft = psd_fft[1:]
    df_fft = freqs_fft[1] - freqs_fft[0]
    
    # Cumulative PSD (Direct FFT)
    cum_psd = np.cumsum(psd_fft[::-1])[::-1] * df_fft
    
    # Target analytical PSD
    C1_target = G_target * (0.1**w_target)
    target_psd = C1_target * freqs**(-w_target)
    
    # Exact cumulative PSD model (using true target parameters and Direct FFT frequency bins)
    exact_cum = exact_isotropic_cum_model(freqs_fft, C1_target, w_target)
    
    # Plotting
    fig, axes = plt.subplots(1, 2, figsize=(15, 6))
    
    # Left: PSD comparison (Target vs Welch)
    ax1 = axes[0]
    ax1.loglog(freqs, psd, color='gray', alpha=0.6, linewidth=0.8, label='Welch PSD (Sliced Profile)')
    ax1.loglog(freqs, target_psd, color='red', linestyle='--', linewidth=2, label=f'Target 1D PSD (w={w_target:.1f}, G={G_target*1e6:.1f} \u03bcm\u00b3)')
    ax1.set_title("Power Spectral Density (PSD) Comparison")
    ax1.set_xlabel("Spatial Frequency (cycles/m)")
    ax1.set_ylabel("PSD (m³)")
    ax1.grid(True, which="both", linestyle='--', alpha=0.5)
    ax1.legend()
    
    # Right: Cumulative PSD comparison
    ax2 = axes[1]
    ax2.loglog(freqs_fft, cum_psd, color='blue', linewidth=2, label='Cumulative PSD (Direct FFT)')
    ax2.loglog(freqs_fft, exact_cum, color='red', linestyle='--', linewidth=2, label='Exact Isotropic Cumulative Model')
    ax2.set_title("Cumulative PSD (Residual Height Variance)")
    ax2.set_xlabel("Spatial Frequency (cycles/m)")
    ax2.set_ylabel("Cumulative Power (m²)")
    ax2.grid(True, which="both", linestyle='--', alpha=0.5)
    ax2.legend()
    
    plt.tight_layout()
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_path = os.path.join(script_dir, "psd_comparison_curves.png")
    plt.savefig(output_path, dpi=150)
    plt.close()
    print(f"Plot saved successfully to: {output_path}", flush=True)
    
    artifact_dir = os.environ.get("ANTIGRAVITY_ARTIFACT_DIR")
    if not artifact_dir:
        artifact_dir = r"C:\Users\novo\.gemini\antigravity\brain\fd4ff96c-fd17-4c02-94be-eb8b0fc6fd62"
    try:
        os.makedirs(artifact_dir, exist_ok=True)
        plt.savefig(os.path.join(artifact_dir, "psd_comparison_curves.png"), dpi=150)
    except Exception as e:
        print(f"Could not save to artifact directory: {e}")

if __name__ == "__main__":
    main()
