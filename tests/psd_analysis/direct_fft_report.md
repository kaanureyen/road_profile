# Road Profile Direct FFT PSD Report

This report presents the discrete PSD calculated using a **direct FFT (Periodogram)** on the entire equally spaced 1D slice of 40,000 points ($dx=0.025$m, length $1000$m).

## Direct FFT PSD Visualizations

![Direct FFT PSD Comparison](psd_direct_fft_comparison.png)

### Key Insights:

1. **Direct FFT PSD (Left):**
   - **Direct FFT PSD (Gray line):** Calculated by taking the FFT of the entire 40,000-point signal (applying a single Hanning window across all 40,000 points).
   - Because the road profile is analytically a sum of discrete sinusoids, the true power spectrum consists of Dirac delta functions. 
   - At the extremely high resolution of $df = 0.001$ cycles/m, the direct FFT displays **extremely sharp spikes (the actual discrete sinusoids)** separated by deep valleys that drop down to the numerical noise floor ($<10^{-15}$ m³).
   - **Welch PSD (Blue line):** Welch's method averages over segments, which smooths out the variance of the spikes but smears their power over adjacent bins.

2. **Cumulative PSD (Right):**
   - **Cumulative Direct FFT PSD (Green line):** Obtained by integrating the direct FFT PSD from high to low frequencies.
   - **Cumulative Welch PSD (Blue dashed line):** The cumulative Welch PSD.
   - Despite the direct FFT PSD being extremely peaky and containing large spikes, its **cumulative PSD (green)** is incredibly smooth and tracks the cumulative Welch PSD (blue) almost perfectly. 
   - This validates that the total power in the discrete sinusoids is identical in both representations and aligns with our analytical model.
