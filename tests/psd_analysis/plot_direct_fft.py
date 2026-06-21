import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# Add tests/ to path to import fmu_helper
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from fmu_helper import FMURoadQuery

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

def main():
    G_target = 64e-6
    w_target = 2.0
    
    # Initialize FMU query helper
    fmu_query = FMURoadQuery()
    slave = fmu_query.get_slave(Gd_n0=G_target, w=w_target, f_min=0.002, f_max=2000.0, Nf=64, Ntheta=16)
    
    slice_length = 1000.0
    dx = 0.025
    N_slice = int(slice_length / dx)
    fs = 1.0 / dx
    
    # Generate 1 slice
    np.random.seed(42)
    x1 = np.random.uniform(-5000.0, 5000.0)
    y1 = np.random.uniform(-5000.0, 5000.0)
    theta = np.random.uniform(0, 2*np.pi)
    
    s = np.linspace(0, slice_length, N_slice, endpoint=False)
    px = x1 + s * np.cos(theta)
    py = y1 + s * np.sin(theta)
    
    # Query road profile using FMU
    z = fmu_query.query_profile(slave, px, py)
    slave.terminate()
    slave.freeInstance()
    
    # 1. Direct FFT (Periodogram) on entire signal
    win_full = np.hanning(len(z))
    win_full_norm = np.sum(win_full**2)
    
    z_detrended = z - np.mean(z)
    z_windowed = z_detrended * win_full
    
    fft_full = np.fft.rfft(z_windowed)
    freqs_fft = np.fft.rfftfreq(len(z), d=dx)
    
    psd_fft = (2.0 / (fs * win_full_norm)) * (np.abs(fft_full)**2)
    
    freqs_fft = freqs_fft[1:]
    psd_fft = psd_fft[1:]
    df_fft = freqs_fft[1] - freqs_fft[0]
    
    cum_psd_fft = np.cumsum(psd_fft[::-1])[::-1] * df_fft
    
    # 2. Welch PSD for comparison
    nperseg = 4096
    freqs_welch, psd_welch = custom_welch(z, fs=fs, nperseg=nperseg)
    freqs_welch = freqs_welch[1:]
    psd_welch = psd_welch[1:]
    df_welch = freqs_welch[1] - freqs_welch[0]
    cum_psd_welch = np.cumsum(psd_welch[::-1])[::-1] * df_welch
    
    # Plotting
    fig, axes = plt.subplots(1, 2, figsize=(15, 6))
    
    ax1 = axes[0]
    ax1.loglog(freqs_fft, psd_fft, color='gray', alpha=0.3, linewidth=0.5, label='Direct FFT PSD (All 40,000 pts)')
    ax1.loglog(freqs_welch, psd_welch, color='blue', alpha=0.8, linewidth=1.2, label='Welch PSD (nperseg=4096)')
    
    C1_target = G_target * (0.1**w_target)
    target_psd = C1_target * freqs_welch**(-w_target)
    ax1.loglog(freqs_welch, target_psd, color='red', linestyle='--', linewidth=2, label='Target PSD')
    
    ax1.set_title("Direct FFT PSD vs Welch PSD")
    ax1.set_xlabel("Spatial Frequency (cycles/m)")
    ax1.set_ylabel("PSD (m³)")
    ax1.grid(True, which="both", linestyle='--', alpha=0.5)
    ax1.legend()
    
    ax2 = axes[1]
    ax2.loglog(freqs_fft, cum_psd_fft, color='green', linewidth=1.5, label='Cumulative Direct FFT PSD')
    ax2.loglog(freqs_welch, cum_psd_welch, color='blue', linestyle='--', linewidth=1.5, label='Cumulative Welch PSD')
    
    ax2.set_title("Cumulative PSD Comparison")
    ax2.set_xlabel("Spatial Frequency (cycles/m)")
    ax2.set_ylabel("Cumulative Power (m²)")
    ax2.grid(True, which="both", linestyle='--', alpha=0.5)
    ax2.legend()
    
    plt.tight_layout()
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_path = os.path.join(script_dir, "psd_direct_fft_comparison.png")
    plt.savefig(output_path, dpi=150)
    plt.close()
    print(f"Plot saved successfully to: {output_path}", flush=True)
    
    artifact_dir = os.environ.get("ANTIGRAVITY_ARTIFACT_DIR")
    if not artifact_dir:
        artifact_dir = r"C:\Users\novo\.gemini\antigravity\brain\eb2516c5-ab48-42c6-b6d8-0b90cc4ca6ca"
    try:
        os.makedirs(artifact_dir, exist_ok=True)
        plt.savefig(os.path.join(artifact_dir, "psd_direct_fft_comparison.png"), dpi=150)
    except Exception as e:
        print(f"Could not save to artifact directory: {e}")

if __name__ == "__main__":
    main()
