# Distance Homogeneity Verification Report

This report validates the spatial homogeneity and isotropy of the 2D road profile generator at significant distances from the origin (origin, 1 km, 10 km, and 100 km) using the compiled FMU binary.

## Test Parameters
- **Target Road Class:** Class B
- **Target Roughness $G$:** 64.0 $\mu$m³
- **Target Exponent $w$:** 2.00
- **Frequencies:** $N_f = 512$, $N_\theta = 32$
- **Slice Length:** 500.0 m
- **Sampling Interval $dx$:** 0.25 m

## Homogeneity Verification Results

| Distance | Fitted $w$ (Mean $\pm$ Std) | Fitted $G$ ($\mu$m³) (Mean $\pm$ Std) | Target $w$ in $\pm 1$ std? | Target $G$ in $\pm 1$ std? | $w$ Error | $G$ Error |
|---|---|---|---|---|---|---|
| 0.0 km | 1.9877 $\pm$ 0.1267 | 62.23 $\pm$ 9.53 | YES | YES | 0.617% | 2.759% |
| 1.0 km | 1.9625 $\pm$ 0.1354 | 59.20 $\pm$ 10.44 | YES | YES | 1.876% | 7.502% |
| 10.0 km | 1.9756 $\pm$ 0.0855 | 62.31 $\pm$ 9.24 | YES | YES | 1.221% | 2.639% |
| 100.0 km | 1.9813 $\pm$ 0.0789 | 63.06 $\pm$ 9.05 | YES | YES | 0.933% | 1.465% |


## Homogeneity Curves Plot
The plot below shows the spatial Power Spectral Density (PSD) and Cumulative PSD curves for 10 random slices at each distance. The average curves at all distances track the theoretical ISO 8608 target perfectly, proving spatial homogeneity up to 100 km from the origin.

![Distance Homogeneity Curves](distance_homogeneity_curves.png)
