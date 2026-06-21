import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from scipy.optimize import curve_fit

# Add tests/ to path to import fmu_helper
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from fmu_helper import FMURoadQuery

def get_I(alpha):
    t = np.linspace(-2000, 2000, 200000)
    dt = t[1] - t[0]
    return np.sum((1.0 + t**2)**(-alpha/2.0)) * dt

# Exact cumulative PSD model with numerical integration of the 2D projection factor
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
    
    # Path configuration matching wavelengths 0.005m and 50m
    slice_length = 500.0
    dx = 0.002
    N_slice = int(slice_length / dx) # 250,000 points
    fs = 1.0 / dx
    
    print(f"Path Configuration: Length = {slice_length}m, dx = {dx}m, N = {N_slice} points", flush=True)
    
    # Slicing parameters
    np.random.seed(42)
    x1 = np.random.uniform(-1000.0, 1000.0)
    y1 = np.random.uniform(-1000.0, 1000.0)
    theta = np.random.uniform(0, 2*np.pi)
    
    s = np.linspace(0, slice_length, N_slice, endpoint=False)
    px = x1 + s * np.cos(theta)
    py = y1 + s * np.sin(theta)
    
    # Initialize FMU query helper
    fmu_query = FMURoadQuery()
    
    # ----------------------------------------------------
    # Task 1: Generate Plots for Different w & G Parameters
    # ----------------------------------------------------
    cases = [
        {'G': 4e-6,    'w': 1.5,  'label': 'Smoother than A (w=1.5)'},
        {'G': 64e-6,   'w': 2.0,  'label': 'Class B (w=2.0)'},
        {'G': 256e-6,  'w': 3.0,  'label': 'Class C (w=3.0)'},
        {'G': 16384e-6,'w': 4.5,  'label': 'Harsher than E (w=4.5)'}
    ]
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    for idx, case in enumerate(cases):
        print(f"Processing Case {idx+1}: G = {case['G']:.2e}, w = {case['w']:.1f}...", flush=True)
        slave = fmu_query.get_slave(Gd_n0=case['G'], w=case['w'], f_min=0.002, f_max=2000.0, Nf=64, Ntheta=16)
        z = fmu_query.query_profile(slave, px, py)
        slave.terminate()
        slave.freeInstance()
        
        # Periodogram PSD with full Hanning window
        win = np.hanning(N_slice)
        win_norm = np.sum(win**2)
        z_detrended = z - np.mean(z)
        z_windowed = z_detrended * win
        
        fft_z = np.fft.rfft(z_windowed)
        freqs = np.fft.rfftfreq(N_slice, d=dx)
        psd = (2.0 / (fs * win_norm)) * (np.abs(fft_z)**2)
        
        freqs = freqs[1:]
        psd = psd[1:]
        df = freqs[1] - freqs[0]
        
        cum_psd = np.cumsum(psd[::-1])[::-1] * df
        
        # Theoretical infinite cumulative PSD
        C1_target = case['G'] * (0.1**case['w'])
        infinite_cum = (C1_target / (case['w'] - 1.0)) * (freqs**(-(case['w'] - 1.0)))
        
        # Exact cumulative model truncated at 2000.0 cycles/m (decimated for speed in plotting)
        plot_idx = np.round(np.logspace(0, np.log10(len(freqs)-1), 300)).astype(int)
        plot_idx = np.unique(np.clip(plot_idx, 0, len(freqs)-1))
        
        exact_cum_dec = exact_isotropic_cum_model(freqs[plot_idx], C1_target, case['w'], f_max=2000.0)
        
        # Plotting in subplots
        ax = axes[idx // 2, idx % 2]
        
        # Cumulative curves
        ax.loglog(freqs, cum_psd, color='blue', linewidth=2, label='Cumulative Direct FFT PSD')
        ax.loglog(freqs, infinite_cum, color='red', linestyle='--', linewidth=1.5, label='Infinite Theoretical Model')
        ax.loglog(freqs[plot_idx], exact_cum_dec, color='green', linestyle=':', linewidth=2, label='Exact Model (Truncated 2000Hz)')
        
        # Shading the target validation range [0.02, 200.0] cycles/m (wavelengths 50m to 0.005m)
        ax.axvspan(0.02, 200.0, color='orange', alpha=0.1, label='Fit Window [0.02, 200] c/m')
        
        ax.set_xlim(0.001, 2000.0)
        ax.set_ylim(1e-15, 1e2)
        ax.set_title(f"Case {idx+1}: {case['label']} | G={case['G']*1e6:.1f} \u03bcm\u00b3")
        ax.set_xlabel("Spatial Frequency (cycles/m)")
        ax.set_ylabel("Cumulative Power (m²)")
        ax.grid(True, which="both", linestyle='--', alpha=0.5)
        ax.legend()
        
    plt.tight_layout()
    output_cases_plot = os.path.join(script_dir, "psd_multi_params.png")
    plt.savefig(output_cases_plot, dpi=150)
    try:
        plt.savefig(os.path.join(artifact_dir, "psd_multi_params.png"), dpi=150)
    except Exception:
        pass
    plt.close()
    print(f"Saved parameters plot to: {output_cases_plot}", flush=True)
    
    # ----------------------------------------------------
    # Task 2: Sensitivity to Number of Bins (Nf)
    # ----------------------------------------------------
    print("Running Nf sensitivity analysis for Case 2...", flush=True)
    Nf_values = [16, 64, 256, 1024]
    case_sens = cases[1] # Case 2: Class B, w=2.0
    
    plt.figure(figsize=(10, 7))
    
    # Pre-calculate exact model for reference
    C1_ref = case_sens['G'] * (0.1**case_sens['w'])
    
    # Direct FFT frequencies
    freqs_ref = np.fft.rfftfreq(N_slice, d=dx)[1:]
    plot_idx = np.round(np.logspace(0, np.log10(len(freqs_ref)-1), 300)).astype(int)
    plot_idx = np.unique(np.clip(plot_idx, 0, len(freqs_ref)-1))
    
    exact_cum_ref = exact_isotropic_cum_model(freqs_ref[plot_idx], C1_ref, case_sens['w'], f_max=2000.0)
    plt.loglog(freqs_ref[plot_idx], exact_cum_ref, color='red', linestyle='--', linewidth=2.5, label='Exact Model (Reference)')
    
    colors = ['purple', 'cyan', 'magenta', 'blue']
    for idx_nf, Nf in enumerate(Nf_values):
        print(f"  Generating surface for Nf = {Nf}...", flush=True)
        slave = fmu_query.get_slave(Gd_n0=case_sens['G'], w=case_sens['w'], f_min=0.002, f_max=2000.0, Nf=Nf, Ntheta=16)
        z = fmu_query.query_profile(slave, px, py)
        slave.terminate()
        slave.freeInstance()
        
        # Periodogram PSD
        win = np.hanning(N_slice)
        win_norm = np.sum(win**2)
        z_detrended = z - np.mean(z)
        z_windowed = z_detrended * win
        
        fft_z = np.fft.rfft(z_windowed)
        freqs_z = np.fft.rfftfreq(N_slice, d=dx)[1:]
        psd_z = (2.0 / (fs * win_norm)) * (np.abs(fft_z[1:])**2)
        df_z = freqs_z[1] - freqs_z[0]
        
        cum_psd_z = np.cumsum(psd_z[::-1])[::-1] * df_z
        
        plt.loglog(freqs_z, cum_psd_z, color=colors[idx_nf], linewidth=1.5, alpha=0.8, label=f'Nf = {Nf} bins')
        
    plt.xlim(0.001, 2000.0)
    plt.ylim(1e-10, 1e-2)
    plt.title(f"Cumulative PSD Sensitivity to Number of Freq Bins (Nf)\nTarget: G = {case_sens['G']*1e6:.1f} \u03bcm\u00b3, w = {case_sens['w']:.1f}")
    plt.xlabel("Spatial Frequency (cycles/m)")
    plt.ylabel("Cumulative Power (m²)")
    plt.grid(True, which="both", linestyle='--', alpha=0.5)
    plt.legend()
    plt.tight_layout()
    
    output_sens_plot = os.path.join(script_dir, "psd_sensitivity_Nf.png")
    plt.savefig(output_sens_plot, dpi=150)
    try:
        plt.savefig(os.path.join(artifact_dir, "psd_sensitivity_Nf.png"), dpi=150)
    except Exception:
        pass
    plt.close()
    print(f"Saved sensitivity plot to: {output_sens_plot}", flush=True)
    print("Done!", flush=True)

if __name__ == "__main__":
    main()
