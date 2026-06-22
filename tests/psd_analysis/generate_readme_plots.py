import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import shutil

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
    print("=== GENERATING README ILLUSTRATION PLOTS ===", flush=True)
    fmu_query = FMURoadQuery()
    
    seed = 42
    road_class = 3  # Class C
    Nf = 512
    Ntheta = 32
    f_min = 0.005
    f_max = 100.0
    
    # 1. 1D C-Type Profile Plot (500m length, dx=0.002m, generated via IFFT with random phases)
    print("Generating 1D Class C profile data via IFFT...", flush=True)
    slice_length = 5000.0
    dx = 0.25
    N_slice = int(slice_length / dx)
    fs = 1.0 / dx
    nperseg = 400
    
    Gd_n0 = 256e-6
    w = 2.0
    n0 = 0.1
    
    freqs_fft = np.fft.rfftfreq(N_slice, d=dx)
    df_fft = freqs_fft[1] - freqs_fft[0]
    
    psd_target = np.zeros_like(freqs_fft)
    psd_target[1:] = Gd_n0 * (freqs_fft[1:] / n0)**(-w)
    psd_target[freqs_fft < f_min] = 0.0
    
    magnitudes = N_slice * np.sqrt(psd_target * df_fft / 2.0)
    
    rng = np.random.RandomState(seed)
    phases = rng.uniform(0, 2*np.pi, len(freqs_fft))
    
    X = magnitudes * np.exp(1j * phases)
    X[0] = 0.0
    if N_slice % 2 == 0:
        X[-1] = np.abs(X[-1])
        
    z_1d = np.fft.irfft(X, n=N_slice)
    s_1d = np.linspace(0, slice_length, N_slice, endpoint=False)
    
    freqs_1d, psd_1d = custom_welch(z_1d, fs=fs, nperseg=nperseg)
    freqs_1d = freqs_1d[1:]
    psd_1d = psd_1d[1:]
    
    # Target Class C PSD (Gd_n0 = 256e-6, w = 2.0)
    C1 = Gd_n0 * (0.1**w)
    target_psd_1d = C1 * (freqs_1d**(-w))
    
    # Plot 1D road profile
    fig1, axes1 = plt.subplots(1, 2, figsize=(15, 6))
    
    # Left Panel: Road Profile Height
    axes1[0].plot(s_1d, z_1d, color='#1f77b4', linewidth=1.0)
    axes1[0].set_title("1D Class C Road Elevation Profile", fontsize=13, fontweight='bold')
    axes1[0].set_xlabel("Longitudinal Distance (m)")
    axes1[0].set_ylabel("Elevation Z (m)")
    axes1[0].grid(True, linestyle='--', alpha=0.5)
    
    # Right Panel: PSD comparison
    axes1[1].loglog(freqs_1d, psd_1d, color='#7f7f7f', alpha=0.8, linewidth=1.2, label='Welch PSD (Generated Profile)')
    axes1[1].loglog(freqs_1d, target_psd_1d, color='#d62728', linestyle='--', linewidth=2.0, label='ISO 8608 Class C Target')
    axes1[1].set_title("Power Spectral Density (PSD) Comparison", fontsize=13, fontweight='bold')
    axes1[1].set_xlabel("Spatial Frequency (cycles/m)")
    axes1[1].set_ylabel("PSD (m³)")
    axes1[1].grid(True, which="both", linestyle='--', alpha=0.5)
    axes1[1].legend()
    
    fig1.suptitle("ISO 8608 Class C Road Profile & Welch PSD Verification", fontsize=15, fontweight='bold', y=0.98)
    plt.tight_layout()
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_1d = os.path.join(script_dir, "readme_1d_road.png")
    plt.savefig(output_1d, dpi=150)
    plt.close(fig1)
    print(f"Saved 1D plot to: {output_1d}", flush=True)
    
    # 2. 2D Surface Map & Offset Slices Isotropy Verification
    # Grid: 500x500m with 2.0m spacing (251x251 points)
    print("Generating 2D Class C elevation map grid...", flush=True)
    x_grid = np.linspace(0.0, 500.0, 251)
    y_grid = np.linspace(0.0, 500.0, 251)
    X_grid, Y_grid = np.meshgrid(x_grid, y_grid)
    
    Z_grid = fmu_query.query_profile_parallel(
        X_grid.ravel(), Y_grid.ravel(), num_threads=8, seed=seed,
        Nf=Nf, Ntheta=Ntheta, f_min=f_min, f_max=f_max, road_class=road_class
    ).reshape(X_grid.shape)
    
    # Define 4 offset slices (not passing through origin)
    # Format: (label, start_x, start_y, angle_deg, color)
    slices_def = [
        ("Slice A (0° at Y=100m)", 0.0, 100.0, 0.0, '#2ca02c'),
        ("Slice B (45° from (50,50)m)", 50.0, 50.0, 45.0, '#ff7f0e'),
        ("Slice C (90° at X=350m)", 350.0, 0.0, 90.0, '#9467bd'),
        ("Slice D (135° from (450,50)m)", 450.0, 50.0, 135.0, '#e377c2')
    ]
    
    slice_len_2d = 5000.0
    dx_2d = 0.25
    N_slice_2d = int(slice_len_2d / dx_2d)
    fs_2d = 1.0 / dx_2d
    nperseg_2d = 400
    
    s_slice = np.linspace(0.0, slice_len_2d, N_slice_2d, endpoint=False) # 500m long slices
    slice_psds = []
    
    fig2 = plt.figure(figsize=(18, 9.5))
    gs = fig2.add_gridspec(2, 3, width_ratios=[1.3, 1.0, 1.0])
    
    # Left Panel: 2D Elevation contour map
    ax_map = fig2.add_subplot(gs[:, 0])
    im = ax_map.imshow(
        Z_grid, 
        extent=[0, 500.0, 0, 500.0], 
        origin='lower',
        cmap='terrain', 
        aspect='equal'
    )
    ax_map.set_title("2D Elevation Map with Slice Trajectories", fontsize=13, fontweight='bold')
    ax_map.set_xlabel("Longitudinal Position X (m)")
    ax_map.set_ylabel("Lateral Position Y (m)")
    cbar = fig2.colorbar(im, ax=ax_map, fraction=0.046, pad=0.04)
    cbar.set_label("Elevation Z (m)")
    
    for label, sx, sy, deg, color in slices_def:
        rad = np.deg2rad(deg)
        # Compute coordinates for plotting trajectory
        px = sx + s_slice * np.cos(rad)
        py = sy + s_slice * np.sin(rad)
        
        # Draw on map (decimated for visual clarity and performance)
        ax_map.plot(px[::100], py[::100], color=color, linewidth=2.5, label=label)
        
        # Query slice profile using FMU
        print(f"Querying {label}...", flush=True)
        z_slice = fmu_query.query_profile_parallel(
            px, py, num_threads=8, seed=seed,
            Nf=Nf, Ntheta=Ntheta, f_min=f_min, f_max=f_max, road_class=road_class
        )
        
        # Compute PSD
        freqs_s, psd_s = custom_welch(z_slice, fs=fs_2d, nperseg=nperseg_2d)
        slice_psds.append((label, freqs_s[1:], psd_s[1:], color))
        
    ax_map.legend(loc='upper right', fontsize=8.5)
    ax_map.grid(True, linestyle='--', alpha=0.4)
    ax_map.set_xlim(0, 500.0)
    ax_map.set_ylim(0, 500.0)
    
    # Right Panels: 4 separate PSD subplots
    ax_psd_list = [
        fig2.add_subplot(gs[0, 1]), # top left of grid
        fig2.add_subplot(gs[0, 2]), # top right of grid
        fig2.add_subplot(gs[1, 1]), # bottom left of grid
        fig2.add_subplot(gs[1, 2])  # bottom right of grid
    ]
    
    for idx, (label, fs_s, ps_s, color) in enumerate(slice_psds):
        ax_psd = ax_psd_list[idx]
        ax_psd.loglog(fs_s, ps_s, color=color, alpha=0.8, linewidth=1.5, label='Welch PSD')
        
        # Theoretical target line
        target_psd_slice = C1 * (fs_s**(-w))
        ax_psd.loglog(fs_s, target_psd_slice, color='black', linestyle='--', linewidth=2.0, label='ISO 8608 Target')
        
        ax_psd.set_title(f"PSD: {label.split(' (')[0]}", fontsize=11, fontweight='bold')
        ax_psd.set_xlabel("Spatial Frequency (cycles/m)", fontsize=9)
        ax_psd.set_ylabel("PSD (m³)", fontsize=9)
        ax_psd.grid(True, which="both", linestyle='--', alpha=0.5)
        ax_psd.legend(loc='lower left', fontsize=8)
    
    fig2.suptitle("2D Spatial Isotropy Verification (Non-Origin Offset Slices)", fontsize=15, fontweight='bold', y=0.98)
    plt.tight_layout()
    
    output_2d = os.path.join(script_dir, "readme_2d_slices.png")
    plt.savefig(output_2d, dpi=150)
    plt.close(fig2)
    print(f"Saved 2D plot to: {output_2d}", flush=True)
    
    # Copy both plots to artifact folder for user display
    artifact_dir = os.environ.get("ANTIGRAVITY_ARTIFACT_DIR")
    if not artifact_dir:
        artifact_dir = r"C:\Users\novo\.gemini\antigravity\brain\fd4ff96c-fd17-4c02-94be-eb8b0fc6fd62"
    try:
        os.makedirs(artifact_dir, exist_ok=True)
        shutil.copy(output_1d, os.path.join(artifact_dir, "readme_1d_road.png"))
        shutil.copy(output_2d, os.path.join(artifact_dir, "readme_2d_slices.png"))
        print(f"Successfully copied readme plots to artifact directory: {artifact_dir}", flush=True)
    except Exception as e:
        print(f"Could not copy readme plots to artifact directory: {e}")

if __name__ == "__main__":
    main()
