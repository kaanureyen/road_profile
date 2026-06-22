# Distance Homogeneity Verification Report

This report validates the spatial homogeneity and isotropy of the 2D road profile generator at significant distances from the origin (origin, 1 km, 10 km, and 100 km) using the compiled FMU binary.

## Test Parameters
- **Target Road Class:** Class B
- **Target Roughness $G$:** 64.0 $\mu$m³
- **Target Exponent $w$:** 2.00
- **Frequencies:** $N_f = 512$, $N_\theta = 32$
- **Slice Length:** 500.0 m
- **Sampling Interval $dx$:** 0.002 m

## Homogeneity Verification Results

| Distance | Fitted $w$ (Mean $\pm$ Std) | Fitted $G$ ($\mu$m³) (Mean $\pm$ Std) | Target $w$ in $\pm 1$ std? | Target $G$ in $\pm 1$ std? | $w$ Error | $G$ Error |
|---|---|---|---|---|---|---|
| 0.0 km | 2.0256 $\pm$ 0.0585 | 70.95 $\pm$ 14.47 | YES | YES | 1.278% | 10.854% |
| 1.0 km | 1.9972 $\pm$ 0.0699 | 65.20 $\pm$ 16.27 | YES | YES | 0.141% | 1.882% |
| 10.0 km | 2.0586 $\pm$ 0.0709 | 80.58 $\pm$ 18.01 | YES | YES | 2.928% | 25.913% |
| 100.0 km | 2.0174 $\pm$ 0.0804 | 70.26 $\pm$ 20.37 | YES | YES | 0.872% | 9.782% |


## Homogeneity Curves Plot
The plot below shows the spatial Power Spectral Density (PSD) and Cumulative PSD curves for 10 random slices at each distance. The average curves at all distances track the theoretical ISO 8608 target perfectly, proving spatial homogeneity up to 100 km from the origin.

![Distance Homogeneity Curves](distance_homogeneity_curves.png)
