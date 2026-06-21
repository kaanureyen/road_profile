import os
import sys
import time
import numpy as np
import matplotlib.pyplot as plt

# Add tests/ to path to import fmu_helper
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from fmu_helper import FMURoadQuery

def main():
    G_target = 64e-6
    w_target = 2.0
    
    slice_length = 500.0
    dx = 0.00025  # Nyquist frequency of 1.0 / (2 * 0.00025) = 2000.0 Hz/cycles/m
    N_slice = int(slice_length / dx)
    
    print(f"Generating road slice of length {slice_length}m with dx={dx}m ({N_slice:,} points)...", flush=True)
    np.random.seed(42)
    x1 = np.random.uniform(-1000.0, 1000.0)
    y1 = np.random.uniform(-1000.0, 1000.0)
    theta = np.random.uniform(0, 2*np.pi)
    
    s = np.linspace(0, slice_length, N_slice, endpoint=False)
    px = x1 + s * np.cos(theta)
    py = y1 + s * np.sin(theta)
    
    print("Evaluating road height using parallel FMU queries...", flush=True)
    t0 = time.time()
    
    fmu_query = FMURoadQuery()
    z = fmu_query.query_profile_parallel(
        px, py, num_threads=8, seed=42, Gd_n0=G_target, w=w_target,
        f_min=0.002, f_max=2000.0, Nf=512, Ntheta=32, road_class=0
    )
    
    elapsed = time.time() - t0
    print(f"Evaluation completed in {elapsed:.2f} seconds.", flush=True)
    
    # Remove mean for FFT
    z_detrended = z - np.mean(z)
    
    # Compute FFT and get single-sided amplitude spectrum
    print("Computing FFT...", flush=True)
    yf = np.fft.rfft(z_detrended)
    freqs = np.fft.rfftfreq(N_slice, d=dx)
    
    # Amplitude normalized by N (factor of 2 for single-sided)
    amplitudes = np.abs(yf) / N_slice * 2.0
    
    # Slice to range 0.002 - 2000 Hz
    mask = (freqs >= 0.002) & (freqs <= 2000.0)
    freqs_plot = freqs[mask]
    amplitudes_plot = amplitudes[mask]
    
    # Plotting both linear and log-log plots
    print("Generating plots...", flush=True)
    fig, axes = plt.subplots(1, 2, figsize=(15, 6))
    
    # Subplot 1: Linear plot
    axes[0].plot(freqs_plot, amplitudes_plot, color='#1f77b4', linewidth=0.5, alpha=0.8)
    axes[0].set_title("Amplitude Spectrum (Linear Scale)")
    axes[0].set_xlabel("Spatial Frequency f (cycles/m or Hz)")
    axes[0].set_ylabel("Amplitude (m)")
    axes[0].grid(True, linestyle='--', alpha=0.5)
    
    # Subplot 2: Log-Log plot
    axes[1].loglog(freqs_plot, amplitudes_plot, color='#9467bd', linewidth=0.5, alpha=0.6)
    # Add a reference line showing the theoretical PSD slope (A ~ f^(-w/2) = f^(-1.0))
    ref_f = np.logspace(np.log10(0.002), np.log10(2000.0), 100)
    ref_amp = 0.0001 * (ref_f**(-1.0))
    axes[1].loglog(ref_f, ref_amp, color='black', linestyle='--', linewidth=1.5, label='Theoretical Amplitude Decay Slope ($f^{-1}$)')
    
    axes[1].set_title("Amplitude Spectrum (Log-Log Scale)")
    axes[1].set_xlabel("Spatial Frequency f (cycles/m or Hz)")
    axes[1].set_ylabel("Amplitude (m)")
    axes[1].grid(True, which="both", linestyle='--', alpha=0.5)
    axes[1].legend()
    
    fig.suptitle(f"FMU Road Profile FFT Amplitude Spectrum (0.002 - 2000 Hz)\nFMU Settings: seed=42, Gd_n0=64e-6, w=2.0, Nf=512, Ntheta=32", fontsize=13, fontweight='bold')
    plt.tight_layout()
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_path = os.path.join(script_dir, "psd_amplitude_frequency.png")
    plt.savefig(output_path, dpi=150)
    plt.close()
    print(f"Plot saved successfully to: {output_path}", flush=True)
    
    # Copy to artifacts directory
    artifact_dir = os.environ.get("ANTIGRAVITY_ARTIFACT_DIR")
    if not artifact_dir:
        artifact_dir = r"C:\Users\novo\.gemini\antigravity\brain\eb2516c5-ab48-42c6-b6d8-0b90cc4ca6ca"
    try:
        os.makedirs(artifact_dir, exist_ok=True)
        plt.savefig(os.path.join(artifact_dir, "psd_amplitude_frequency.png"), dpi=150)
    except Exception as e:
        print(f"Could not save to artifact directory: {e}")

if __name__ == "__main__":
    main()
