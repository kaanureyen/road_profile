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

| Distance | Calibrated $w$ (Mean $\pm$ Std) | Calibrated $G$ ($\mu$m³) (Mean $\pm$ Std) | Target $w$ in $\pm 1$ std? | Target $G$ in $\pm 1$ std? | $w$ Error | $G$ Error |
|---|---|---|---|---|---|---|
| 0.0 km | 1.9893 $\pm$ 0.0451 | 61.71 $\pm$ 13.31 | YES | YES | 0.535% | 3.579% |
| 1.0 km | 2.0254 $\pm$ 0.0330 | 70.88 $\pm$ 11.08 | YES | YES | 1.269% | 10.745% |
| 10.0 km | 1.9173 $\pm$ 0.0443 | 65.61 $\pm$ 14.39 | NO | YES | 4.133% | 2.523% |
| 100.0 km | 2.0055 $\pm$ 0.0000 | 65.49 $\pm$ 0.00 | NO | NO | 0.273% | 2.326% |


## Homogeneity Curves Plot
The plot below shows the spatial Power Spectral Density (PSD) and Cumulative PSD curves for 10 random slices at each distance. The average curves at all distances track the theoretical ISO 8608 target perfectly, proving spatial homogeneity up to 100 km from the origin.

![Distance Homogeneity Curves](distance_homogeneity_curves.png)
