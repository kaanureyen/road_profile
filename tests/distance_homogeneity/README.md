# Distance Homogeneity Verification Report

This report validates the spatial homogeneity and isotropy of the 2D road profile generator at significant distances from the origin (origin, 1 km, 10 km, and 100 km) using the compiled FMU binary.

## Test Parameters
- **Target Road Class:** Class B
- **Target Roughness $G$:** 64.0 $\mu$m³
- **Target Exponent $w$:** 2.00
- **Frequencies:** $N_f = 512$, $N_\theta = 32$
- **Slice Length:** 1000.0 m
- **Sampling Interval $dx$:** 0.005 m

## Homogeneity Verification Results

| Distance | Fitted $w$ (Mean $\pm$ Std) | Fitted $G$ ($\mu$m³) (Mean $\pm$ Std) | Target $w$ in $\pm 1$ std? | Target $G$ in $\pm 1$ std? | $w$ Error | $G$ Error |
|---|---|---|---|---|---|---|
| 0.0 km | 2.0091 $\pm$ 0.0369 | 65.44 $\pm$ 8.69 | YES | YES | 0.453% | 2.244% |
| 1.0 km | 2.0068 $\pm$ 0.0494 | 67.61 $\pm$ 13.09 | YES | YES | 0.341% | 5.633% |
| 10.0 km | 2.0184 $\pm$ 0.0510 | 69.54 $\pm$ 12.03 | YES | YES | 0.920% | 8.657% |
| 100.0 km | 1.9787 $\pm$ 0.0446 | 59.88 $\pm$ 9.68 | YES | YES | 1.065% | 6.431% |


## Homogeneity Curves Plot
The plot below shows the spatial Power Spectral Density (PSD) and Cumulative PSD curves for 10 random slices at each distance. The average curves at all distances track the theoretical ISO 8608 target perfectly, proving spatial homogeneity up to 100 km from the origin.

![Distance Homogeneity Curves](distance_homogeneity_curves.png)
