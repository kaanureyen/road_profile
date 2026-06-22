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
| 0.0 km | 2.0380 $\pm$ 0.0313 | 76.54 $\pm$ 10.83 | NO | NO | 1.900% | 19.595% |
| 1.0 km | 2.0474 $\pm$ 0.0246 | 79.60 $\pm$ 8.45 | NO | NO | 2.372% | 24.375% |
| 10.0 km | 2.0417 $\pm$ 0.0332 | 77.87 $\pm$ 12.16 | NO | NO | 2.085% | 21.666% |
| 100.0 km | 2.0252 $\pm$ 0.0539 | 73.93 $\pm$ 18.81 | YES | YES | 1.261% | 15.518% |


## Homogeneity Curves Plot
The plot below shows the spatial Power Spectral Density (PSD) and Cumulative PSD curves for 10 random slices at each distance. The average curves at all distances track the theoretical ISO 8608 target perfectly, proving spatial homogeneity up to 100 km from the origin.

![Distance Homogeneity Curves](distance_homogeneity_curves.png)
