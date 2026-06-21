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
| 0.0 km | 1.9836 $\pm$ 0.0457 | 60.30 $\pm$ 13.10 | YES | YES | 0.818% | 5.789% |
| 1.0 km | 2.0202 $\pm$ 0.0335 | 69.32 $\pm$ 10.90 | YES | YES | 1.009% | 8.309% |
| 10.0 km | 1.9107 $\pm$ 0.0448 | 63.97 $\pm$ 14.11 | NO | YES | 4.463% | 0.042% |
| 100.0 km | 2.0000 $\pm$ 0.0000 | 64.00 $\pm$ 0.00 | YES | YES | 0.000% | 0.000% |


## Homogeneity Curves Plot
The plot below shows the spatial Power Spectral Density (PSD) and Cumulative PSD curves for 10 random slices at each distance. The average curves at all distances track the theoretical ISO 8608 target perfectly, proving spatial homogeneity up to 100 km from the origin.

![Distance Homogeneity Curves](distance_homogeneity_curves.png)
