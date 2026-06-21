# Advanced PSD and Cumulative PSD Analysis Report

This report validates the deterministic 2D isotropic road profile generator under extreme frequency scales, varying exponent parameters ($w \in [1.5, 4.5]$), roughness indices ($G \in [4, 16384] \mu\text{m}^3$), frequency bin densities ($N_f$), and angular directional bin densities ($N_\theta$).

---

## 1. Parameter Grid Validation ($w$ and $G$)

### Slicing Setup:
- **Slice Length:** 500m
- **Spatial Resolution ($dx$):** 0.002m (sampling rate $500\text{ Hz}$)
- **Bandwidth Limits:** $f_{\min} = 0.002\text{ cycles/m}$ to $f_{\max} = 2000.0\text{ cycles/m}$
- **Fit Window:** $[0.02, 200.0]\text{ cycles/m}$, corresponding to wavelengths from $50\text{m}$ down to $0.005\text{m}$.

### Parameter Plots (Infinite vs. Truncated Cumulative PSD)
The subplots below show how the cumulative PSD tracks the target curves for each parameter configuration:

![Parameter Comparison](C:/Users/novo/.gemini/antigravity/brain/bdb3005b-89b5-4dd1-a5e2-fedf6fb87855/psd_multi_params.png)

> [!NOTE]
> **Tail Drop-Off Effect:** Slicing through the 2D surface truncates the high-frequency integration tail. As a result, the physical cumulative PSD (blue) drops below the infinite theoretical model (red) near the high-frequency cutoff. The exact model (green) matches the physical curve perfectly because it integrates up to $2000.0\text{ cycles/m}$ only.

---

## 2. Sensitivity Analysis to Number of Freq Bins ($N_f$)

The number of radial frequency bins $N_f$ defines how densely the log-spaced wave rings are generated. We evaluated the sensitivity of the raw PSD and cumulative PSD for $N_f \in \{16, 64, 256, 1024\}$ (keeping $N_\theta = 16$ constant):

![Nf Sensitivity](C:/Users/novo/.gemini/antigravity/brain/bdb3005b-89b5-4dd1-a5e2-fedf6fb87855/psd_sensitivity_Nf_raw_cum.png)

### Key Observations:
1. **Raw PSD (Left):**
   - At $N_f = 16$, the raw PSD consists of a few isolated spikes with deep empty valleys. 
   - As $N_f$ increases to $256$ and $1024$, the frequency grid becomes denser, filling in the valleys and forming a much more continuous-looking spectrum that tracks the Target 1D PSD (red dashed line).
2. **Cumulative PSD (Right):**
   - At $N_f = 16$, the cumulative PSD displays huge, coarse "stairs" because the integration sums only a few large discrete steps.
   - At $N_f = 64$ (FMU default), the steps are smaller but still visible as minor ripples.
   - At $N_f = 256$, the stairs disappear, and the curve tracks the exact theoretical model with high precision.
   - At $N_f = 1024$, the cumulative PSD is a smooth, continuous line matching the target projection curve perfectly.

---

## 3. Sensitivity Analysis to Number of Directions ($N_\theta$)

The number of angular bins $N_\theta$ defines how many directions are used to distribute the wave components over the $2\pi$ circle. We evaluated the sensitivity of the raw PSD and cumulative PSD for $N_\theta \in \{4, 8, 16, 64\}$ (keeping $N_f = 64$ constant):

![Ntheta Sensitivity](C:/Users/novo/.gemini/antigravity/brain/bdb3005b-89b5-4dd1-a5e2-fedf6fb87855/psd_sensitivity_Ntheta_raw_cum.png)

### Key Observations:
1. **Raw PSD (Left):**
   - Since a 1D slice projects 2D waves onto a line ($f_p = f_c \cos(\theta_j - \theta_s)$), having a small number of angular bins $N_\theta$ means the projected frequencies are highly clumped. For $N_\theta = 4$, waves are spaced at $90^\circ$, causing some components to project near zero, creating massive gaps in the 1D spectrum.
   - Increasing $N_\theta \to 64$ provides a dense, continuous angular distribution, ensuring that waves project smoothly onto the 1D slice, tracking the target PSD profile.
2. **Cumulative PSD (Right):**
   - For $N_\theta = 4$ and $N_\theta = 8$, the cumulative PSD displays significant local deviations from the exact model because the directional sparsity causes power to be concentrated in only a few projected frequencies.
   - For $N_\theta = 16$ (FMU default) and $N_\theta = 64$, the cumulative PSD tracks the exact projection model smoothly, showing that $N_\theta = 16$ is sufficient to achieve isotropic power distribution along the 1D slice.

---

## Conclusion & Design Guidelines

> [!IMPORTANT]
> - **Simulating Isotropic Surfaces:** The number of frequency bins $N_f$ and angular bins $N_\theta$ represent the spatial grid discretization of the 2D surface.
> - **Discretization Ripples:** For high-fidelity vehicle simulations where the vehicle suspension is sensitive to minor road profile ripples, using **$N_f \ge 256$** and **$N_\theta \ge 16$** is recommended to completely eliminate discretization step-ripples in the heights and PSD.
