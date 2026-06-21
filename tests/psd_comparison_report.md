# Road Profile PSD Curves Comparison

This report displays the comparison of the target analytical road PSD, the Welch-estimated PSD from a single 1D slice, and the cumulative PSD curves used in our unbiased fitting methodology.

## PSD Visualizations

![PSD Comparison Curves](C:/Users/novo/.gemini/antigravity/brain/bdb3005b-89b5-4dd1-a5e2-fedf6fb87855/psd_comparison_curves.png)

### Key Observations from the Curves:

1. **Power Spectral Density Comparison (Left):**
   - **Target 1D PSD (Red dashed line):** Shows the continuous target scaling $S_{1D}(f) = C_1 f^{-2.0}$.
   - **Welch PSD (Gray line):** Slicing through the 2D surface (which is composed of $1024$ discrete wave sinusoids) yields a **discrete line spectrum**. The Welch estimate smears these discrete lines due to the Hanning window, but leaves many bins nearly empty (low values where the spectrum drops between the discrete components). A direct logarithmic fit to this curve is heavily biased by these empty valleys.

2. **Cumulative PSD Comparison (Right):**
   - **Cumulative Welch PSD (Blue solid line):** Integrating the Welch PSD from high to low frequencies calculates the residual variance $\Phi(f) = \sum_{f_k \ge f} S_{1D}(f_k) \cdot df$. This integration naturally smooths out the discrete line valleys, forming a robust and smooth monotonic curve.
   - **Exact Isotropic Cumulative Model (Red dashed line):** Represents the analytical isotropic projection integral. It perfectly matches the cumulative Welch PSD, showing how the cumulative integration successfully recovers the target isotropic power law parameters without any discretization bias.
   - **Fit Window (Orange shading):** The fit window $[0.1, 3.0]$ cycles/m completely avoids the low-frequency detrending attenuation (Welch segment length limit) and high-frequency sparsity (Nyquist limit), ensuring a clean, calibrated fit.
