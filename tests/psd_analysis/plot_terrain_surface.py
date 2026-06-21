import os
import sys
import time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

# Add tests/ to path to import fmu_helper
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from fmu_helper import FMURoadQuery

def build_plots(X_macro, Y_macro, Z_macro, X_sub, Y_sub, Z_sub, seed, Nf, Ntheta, f_min, f_max):
    """Generates the combined 2D and 3D terrain plots."""
    fig = plt.figure(figsize=(18, 8))
    
    # Subplot 1: 2D Topographic contour map (Entire 500x500m surface)
    ax1 = fig.add_subplot(1, 2, 1)
    im = ax1.imshow(
        Z_macro, 
        extent=[0, 500.0, 0, 500.0], 
        origin='lower',
        cmap='terrain', 
        aspect='equal'
    )
    ax1.set_title("2D Terrain Elevation Map (500m x 500m)", fontsize=14, fontweight='bold', pad=10)
    ax1.set_xlabel("Longitudinal Distance X (m)", fontsize=12)
    ax1.set_ylabel("Lateral Distance Y (m)", fontsize=12)
    cbar = fig.colorbar(im, ax=ax1, fraction=0.046, pad=0.04)
    cbar.set_label("Elevation Z (m)", fontsize=12)
    ax1.grid(True, linestyle='--', alpha=0.4)
    
    # Subplot 2: 3D Surface Mesh of a high-resolution 100x100m subset
    ax2 = fig.add_subplot(1, 2, 2, projection='3d')
    surf = ax2.plot_surface(
        X_sub, Y_sub, Z_sub, 
        cmap='terrain',
        edgecolor='none', 
        alpha=0.9,
        antialiased=True
    )
    
    ax2.set_title("3D Detailed Terrain View (100m x 100m Subset)", fontsize=14, fontweight='bold', pad=10)
    ax2.set_xlabel("X (m)", fontsize=12)
    ax2.set_ylabel("Y (m)", fontsize=12)
    ax2.set_zlabel("Z (m)", fontsize=12)
    ax2.view_init(elev=35, azim=-45)
    
    fig.suptitle("ISO 8608 Isotropic 2D Road Terrain Profile Simulation\n"
                 f"FMU Class C Configuration (seed={seed}, Nf={Nf}, Ntheta={Ntheta}, f_min={f_min}, f_max={f_max} cycles/m)",
                 fontsize=16, fontweight='bold', y=0.98)
    
    plt.tight_layout()
    return fig

def main():
    print("=== SAMPLING 500x500m TERRAIN SURFACE USING FMU ===")
    
    fmu_query = FMURoadQuery()
    
    seed = 42
    road_class = 3  # Class C
    Nf = 512
    Ntheta = 32
    f_min = 0.002
    f_max = 2000.0
    
    t0 = time.time()
    
    # 1. Macro Grid: 500x500m with 2.0m spacing (251x251 points = 63,001 coords)
    print("Evaluating 500x500m macro grid (2.0m spacing)...", flush=True)
    x_macro = np.linspace(0.0, 500.0, 251)
    y_macro = np.linspace(0.0, 500.0, 251)
    X_macro, Y_macro = np.meshgrid(x_macro, y_macro)
    
    Z_macro = fmu_query.query_profile_parallel(
        X_macro.ravel(), Y_macro.ravel(), num_threads=8, seed=seed, 
        Nf=Nf, Ntheta=Ntheta, f_min=f_min, f_max=f_max, road_class=road_class
    ).reshape(X_macro.shape)
    
    # 2. Detailed Grid: 100x100m with 0.5m spacing (201x201 points = 40,401 coords)
    print("Evaluating 100x100m detailed subset (0.5m spacing)...", flush=True)
    x_sub = np.linspace(200.0, 300.0, 201)
    y_sub = np.linspace(200.0, 300.0, 201)
    X_sub, Y_sub = np.meshgrid(x_sub, y_sub)
    
    Z_sub = fmu_query.query_profile_parallel(
        X_sub.ravel(), Y_sub.ravel(), num_threads=8, seed=seed, 
        Nf=Nf, Ntheta=Ntheta, f_min=f_min, f_max=f_max, road_class=road_class
    ).reshape(X_sub.shape)
    
    elapsed = time.time() - t0
    print(f"All grids evaluated in {elapsed:.2f} seconds.", flush=True)
    
    # Generate and save plots
    print("Generating terrain plots...", flush=True)
    fig = build_plots(X_macro, Y_macro, Z_macro, X_sub, Y_sub, Z_sub, seed, Nf, Ntheta, f_min, f_max)
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_path = os.path.join(script_dir, "road_terrain_500x500.png")
    plt.savefig(output_path, dpi=150)
    plt.close(fig)
    print(f"Terrain plot saved successfully to: {output_path}", flush=True)
    
    # Copy to artifact directory
    artifact_dir = os.environ.get("ANTIGRAVITY_ARTIFACT_DIR")
    if not artifact_dir:
        artifact_dir = r"C:\Users\novo\.gemini\antigravity\brain\eb2516c5-ab48-42c6-b6d8-0b90cc4ca6ca"
    try:
        os.makedirs(artifact_dir, exist_ok=True)
        fig_art = build_plots(X_macro, Y_macro, Z_macro, X_sub, Y_sub, Z_sub, seed, Nf, Ntheta, f_min, f_max)
        plt.savefig(os.path.join(artifact_dir, "road_terrain_500x500.png"), dpi=150)
        plt.close(fig_art)
        print("Copied terrain plot to artifact directory.", flush=True)
    except Exception as e:
        print(f"Could not copy terrain plot to artifact directory: {e}")

if __name__ == "__main__":
    main()
