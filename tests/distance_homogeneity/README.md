# Distance Homogeneity Verification Report

This report validates the spatial homogeneity and isotropy of the 2D road profile generator at significant distances from the origin (origin, 1 km, 10 km, and 100 km) using the compiled FMU binary.

## Test Configuration
- **Slice Length:** 500.0 m
- **Sampling Interval $dx$:** 0.002 m (250,000 points per slice)
- **Slices per Distance:** 10 random slices

## FMU Parameters Used
- **`seed`**: `int(dist) + slice_idx + 2026` (to keep random realizations deterministic and repeatable)
- **`road_class`**: `0` (custom $G_d(n_0)$ target)
- **`Gd_n0`**: `64e-6` ($64\ \mu\text{m}^3$, Class B target)
- **`w`**: `2.0` (target exponent)
- **`f_min`**: `0.002` (lower cutoff frequency)
- **`f_max`**: `2000.0` (upper cutoff frequency)
- **`Nf`**: `512` (radial frequency rings)
- **`Ntheta`**: `32` (angular divisions)

## Homogeneity Verification Results

| Distance | Calibrated $w$ (Mean $\pm$ Std) | Calibrated $G$ ($\\mu$m³) (Mean $\pm$ Std) | Target $w$ in $\pm 1$ std? | Target $G$ in $\pm 1$ std? | $w$ Error | $G$ Error |
|---|---|---|---|---|---|---|
| 0.0 km | 2.0169 $\pm$ 0.0465 | 64.09 $\pm$ 5.16 | YES | YES | 0.846% | 0.147% |
| 1.0 km | 2.0202 $\pm$ 0.0531 | 64.55 $\pm$ 6.30 | YES | YES | 1.009% | 0.865% |
| 10.0 km | 2.0175 $\pm$ 0.0494 | 64.12 $\pm$ 5.43 | YES | YES | 0.876% | 0.191% |
| 100.0 km | 2.0191 $\pm$ 0.0512 | 64.33 $\pm$ 5.92 | YES | YES | 0.957% | 0.518% |


## Homogeneity Curves Plot
The plot below shows the spatial Power Spectral Density (PSD) and Cumulative PSD curves for 10 random slices at each distance. The average curves at all distances track the theoretical ISO 8608 target perfectly, proving spatial homogeneity up to 100 km from the origin.

![Distance Homogeneity Curves](distance_homogeneity_curves.png)
