import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate

# Add tests/ to path to import fmu_helper
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from fmu_helper import FMURoadQuery

def get_I(alpha):
    t = np.linspace(-2000, 2000, 200000)
    dt = t[1] - t[0]
    return np.sum((1.0 + t**2)**(-alpha/2.0)) * dt

# Exact cumulative PSD model
def exact_isotropic_cum_model(f_array, C1, w, f_max=2000.0):
    alpha = w + 1.0
    I_val = get_I(alpha)
    results = []
    for f_val in f_array:
        val, _ = integrate.quad(lambda f_2D: (f_2D**(-w)) * np.arccos(f_val / f_2D), f_val, f_max)
        results.append(C1 * (2.0 / I_val) * val)
    return np.array(results)

def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.makedirs(script_dir, exist_ok=True)
    
    artifact_dir = os.environ.get("ANTIGRAVITY_ARTIFACT_DIR")
    if not artifact_dir:
        artifact_dir = r"C:\Users\novo\.gemini\antigravity\brain\eb2516c5-ab48-42c6-b6d8-0b90cc4ca6ca"
    try:
        os.makedirs(artifact_dir, exist_ok=True)
    except Exception:
        pass
    
    # Path configuration
    slice_length = 500.0
    dx = 0.002
    N_slice = int(slice_length / dx) # 250,000 points
    fs = 1.0 / dx
    
    np.random.seed(42)
    x1 = np.random.uniform(-1000.0, 1000.0)
    y1 = np.random.uniform(-1000.0, 1000.0)
    theta = np.random.uniform(0, 2*np.pi)
    
    s = np.linspace(0, slice_length, N_slice, endpoint=False)
    px = x1 + s * np.cos(theta)
    py = y1 + s * np.sin(theta)
    
    G_target = 64e-6
    w_target = 2.0
    C1_target = G_target * (0.1**w_target)
    
    # Initialize FMU query helper
    fmu_query = FMURoadQuery()
    
    # ----------------------------------------------------
    # Study 1: Nf Sensitivity (Frequency Bins) - Raw & Cumulative
    # ----------------------------------------------------
    print("Running Nf study...", flush=True)
    Nf_values = [16, 64, 256, 1024]
    colors = ['purple', 'cyan', 'magenta', 'blue']
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))
    
    # Pre-calculate target PSD
    freqs_ref = np.fft.rfftfreq(N_slice, d=dx)[1:]
    target_psd = C1_target * freqs_ref**(-w_target)
    ax1.loglog(freqs_ref, target_psd, color='red', linestyle='--', linewidth=2, label='Target 1D PSD')
    
    # Pre-calculate exact cumulative model
    plot_idx = np.round(np.logspace(0, np.log10(len(freqs_ref)-1), 300)).astype(int)
    plot_idx = np.unique(np.clip(plot_idx, 0, len(freqs_ref)-1))
    exact_cum = exact_isotropic_cum_model(freqs_ref[plot_idx], C1_target, w_target, f_max=2000.0)
    ax2.loglog(freqs_ref[plot_idx], exact_cum, color='red', linestyle='--', linewidth=2, label='Exact Model')
    
    for idx, Nf in enumerate(Nf_values):
        print(f"  Nf = {Nf}...", flush=True)
        slave = fmu_query.get_slave(Gd_n0=G_target, w=w_target, f_min=0.002, f_max=2000.0, Nf=Nf, Ntheta=16)
        z = fmu_query.query_profile(slave, px, py)
        slave.terminate()
        slave.freeInstance()
        
        # Periodogram PSD
        win = np.hanning(N_slice)
        win_norm = np.sum(win**2)
        z_detrended = z - np.mean(z)
        z_windowed = z_detrended * win
        
        fft_z = np.fft.rfft(z_windowed)
        psd = (2.0 / (fs * win_norm)) * (np.abs(fft_z[1:])**2)
        df = freqs_ref[1] - freqs_ref[0]
        cum_psd = np.cumsum(psd[::-1])[::-1] * df
        
        # Plot Raw PSD (decimated for visual clarity, taking every 5th point)
        ax1.loglog(freqs_ref[::5], psd[::5], color=colors[idx], alpha=0.6, linewidth=0.7, label=f'Nf = {Nf}')
        # Plot Cumulative PSD
        ax2.loglog(freqs_ref, cum_psd, color=colors[idx], linewidth=1.5, label=f'Nf = {Nf}')
        
    ax1.set_xlim(0.001, 2000.0)
    ax1.set_ylim(1e-18, 1e2)
    ax1.set_title("Raw PSD: Frequency Bin Study (Nf)")
    ax1.set_xlabel("Spatial Frequency (cycles/m)")
    ax1.set_ylabel("PSD (m³)")
    ax1.grid(True, which="both", linestyle='--', alpha=0.5)
    ax1.legend()
    
    ax2.set_xlim(0.001, 2000.0)
    ax2.set_ylim(1e-10, 1e-2)
    ax2.set_title("Cumulative PSD: Frequency Bin Study (Nf)")
    ax2.set_xlabel("Spatial Frequency (cycles/m)")
    ax2.set_ylabel("Cumulative Power (m²)")
    ax2.grid(True, which="both", linestyle='--', alpha=0.5)
    ax2.legend()
    
    plt.tight_layout()
    output_Nf_plot = os.path.join(script_dir, "psd_sensitivity_Nf_raw_cum.png")
    plt.savefig(output_Nf_plot, dpi=150)
    try:
        plt.savefig(os.path.join(artifact_dir, "psd_sensitivity_Nf_raw_cum.png"), dpi=150)
    except Exception:
        pass
    plt.close()
    print(f"Saved Nf plot to: {output_Nf_plot}", flush=True)
    
    # ----------------------------------------------------
    # Study 2: Ntheta Sensitivity (Directional/Angular Bins) - Raw & Cumulative
    # ----------------------------------------------------
    print("Running Ntheta study...", flush=True)
    Ntheta_values = [4, 8, 16, 64]
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))
    
    # Target and Exact lines
    ax1.loglog(freqs_ref, target_psd, color='red', linestyle='--', linewidth=2, label='Target 1D PSD')
    ax2.loglog(freqs_ref[plot_idx], exact_cum, color='red', linestyle='--', linewidth=2, label='Exact Model')
    
    for idx, Ntheta in enumerate(Ntheta_values):
        print(f"  Ntheta = {Ntheta}...", flush=True)
        slave = fmu_query.get_slave(Gd_n0=G_target, w=w_target, f_min=0.002, f_max=2000.0, Nf=64, Ntheta=Ntheta)
        z = fmu_query.query_profile(slave, px, py)
        slave.terminate()
        slave.freeInstance()
        
        # Periodogram PSD
        win = np.hanning(N_slice)
        win_norm = np.sum(win**2)
        z_detrended = z - np.mean(z)
        z_windowed = z_detrended * win
        
        fft_z = np.fft.rfft(z_windowed)
        psd = (2.0 / (fs * win_norm)) * (np.abs(fft_z[1:])**2)
        df = freqs_ref[1] - freqs_ref[0]
        cum_psd = np.cumsum(psd[::-1])[::-1] * df
        
        # Plot Raw PSD (decimated)
        ax1.loglog(freqs_ref[::5], psd[::5], color=colors[idx], alpha=0.6, linewidth=0.7, label=f'Ntheta = {Ntheta}')
        # Plot Cumulative PSD
        ax2.loglog(freqs_ref, cum_psd, color=colors[idx], linewidth=1.5, label=f'Ntheta = {Ntheta}')
        
    ax1.set_xlim(0.001, 2000.0)
    ax1.set_ylim(1e-18, 1e2)
    ax1.set_title("Raw PSD: Angular Directional Bin Study (Ntheta)")
    ax1.set_xlabel("Spatial Frequency (cycles/m)")
    ax1.set_ylabel("PSD (m³)")
    ax1.grid(True, which="both", linestyle='--', alpha=0.5)
    ax1.legend()
    
    ax2.set_xlim(0.001, 2000.0)
    ax2.set_ylim(1e-10, 1e-2)
    ax2.set_title("Cumulative PSD: Angular Directional Bin Study (Ntheta)")
    ax2.set_xlabel("Spatial Frequency (cycles/m)")
    ax2.set_ylabel("Cumulative Power (m²)")
    ax2.grid(True, which="both", linestyle='--', alpha=0.5)
    ax2.legend()
    
    plt.tight_layout()
    output_Ntheta_plot = os.path.join(script_dir, "psd_sensitivity_Ntheta_raw_cum.png")
    plt.savefig(output_Ntheta_plot, dpi=150)
    try:
        plt.savefig(os.path.join(artifact_dir, "psd_sensitivity_Ntheta_raw_cum.png"), dpi=150)
    except Exception:
        pass
    plt.close()
    print(f"Saved Ntheta plot to: {output_Ntheta_plot}", flush=True)

if __name__ == "__main__":
    main()
