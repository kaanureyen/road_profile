# Cumulative PSD Comparison: Infinite vs. Truncated Isotropic Model

This report compares the Welch cumulative PSD estimated from a **10,000m slice** down to $0.0001\text{ cycles/m}$ (1/20th of the lower cutoff frequency $f_{\min} = 0.002\text{ cycles/m}$) against:
1. **The Infinite Theoretical Model (No Cutoff):** $\Phi(f) = \frac{C_1}{w-1} f^{-(w-1)}$
2. **The Exact Model (Truncated at $20.0\text{ cycles/m}$):** Integrates the isotropic projected PSD up to the physical maximum generator limit.

## PSD Visualizations

![Infinite PSD Comparison](psd_infinite_comparison.png)

### Key Observations:

1. **Below the Lower Cutoff ($f < f_{\min} = 0.002$ cycles/m):**
   - **The physical lower limit** of the generator is $0.002\text{ cycles/m}$ (orange dashed line). No sinusoids are generated below this frequency.
   - Accordingly, the **Welch Cumulative PSD (Blue line)** is **completely flat** (horizontal) for $f \le 0.002$, indicating zero power is added as we integrate further to the left.
   - **The Infinite Theoretical Model (Red dashed line)** continues to grow towards infinity ($\Phi(f) \to \infty$ as $f \to 0$), severely deviating from the actual physical signal at low frequencies.
   - **The Exact Model (Green dotted line)** perfectly tracks this flat profile, showing that it correctly accounts for the band-limited nature of the physical system.

2. **At High Frequencies ($f > 3.0$ cycles/m):**
   - The **Welch Cumulative PSD (Blue line)** drops off steeply and approaches zero near the Nyquist frequency of $20.0\text{ cycles/m}$.
   - The **Infinite Theoretical Model (Red dashed line)** does *not* drop off because it assumes infinite tail energy ($f \to \infty$), causing it to lie above the true curve at high frequencies.
   - The **Exact Model (Green dotted line)** captures this truncation tail drop-off perfectly, tracking the blue Welch curve precisely to zero.
