# corradetti-2026-roughness

MATLAB pipeline for PSD-based fault surface roughness characterization — Hurst exponent (H) and spectral pre-factor (alpha) — supporting:

> Corradetti, A., Romanelli, F., Cremonesi, A., Bonini, L., Agosta, F. "Beyond the Hurst Exponent: The Role of the Spectral Pre-factor in Fault Roughness Characterisation and Its Synthetic Validation" (under review).

## Overview

Two independent entry points depending on the input data:

- **Synthetic self-affine surfaces** (already on a regular grid, no pre-processing needed) → `PC_Master_synthetic.m`
- **Natural point clouds** (e.g. SfM-MVS models) requiring alignment and interpolation → `PC_Master_natural.m`

Both paths converge on the same PSD fitting routines (three variants: allmean, progmean, userbounded progmean), ensuring synthetic and natural surfaces are processed identically.

## Scripts

### Surface generation & noise
- `generateSelfAffineSurface.m` — generates anisotropic self-affine surfaces with target Hx, Hy and amplitude (Rq or P2V), via a separable filter in the Fourier domain.
- `noise_developer.m` — adds band-limited high-frequency noise (isotropic radial high-pass filter) to a generated surface, with a user-defined cutoff wavelength and target SNR.

### Alignment & interpolation (natural point clouds only)
- `import_PC.m` — generic point cloud importer (.csv, .xyz, .txt, .ply).
- `alignZ.m` — rotates the point cloud so Z is normal to the least-squares best-fit plane (SVD).
- `Figura_raw.m` — visualization utility for the raw, Z-aligned point cloud.
- `alignXY.m` — further rotation around Z to align X/Y with the slip direction.
- `Interpolate_PC.m` — interpolates the scattered, rotated point cloud onto a regular grid.
- `make_square_surface.m` — crops the interpolated grid to a square extent.

### FFT / PSD analysis
- `PrepareSurfaceForFFT.m` — direct FFT path for surfaces already on a regular grid (synthetic surfaces).
- `processFFT_XAxis.m` / `processFFT_YAxis.m` — computes the 1D power spectral density of every profile along X/Y (3% linear taper, PSD normalized by the frequency resolution).

### Spectral fitting variants
Applied independently along X (perpendicular to slip) and Y (parallel to slip):
- `XFit_all_data.m` / `YFit_all_data.m` — regression on all raw spectral estimates (**allmean**).
- `XFit_progre_mean.m` / `YFit_progre_mean.m` — geometric-progression binning before regression (**progmean**, following Candela et al. 2012).
- `XFit_userbounded_progmean.m` / `YFit_userbounded_progmean.m` — same binning restricted to a user-defined frequency window (**userbounded progmean**).
- `freqRangeGUI.m` — GUI helper for selecting the frequency window used by the userbounded variants.

Each fit routine also computes and exports per-frequency spectral dispersion (Std and IQR of log10 P(k)) and the theoretical RMS roughness (Rq_theoretical).

### Multi-surface comparison
- `compare_profiles_multizoom_stessi_assi.m` — compares central profiles from up to 4 surfaces at two zoom levels (100%, 10%) with shared Z axes, computing single-profile H, alpha and Rq_theoretical for each panel.

### Master scripts
- `PC_Master_synthetic.m` — end-to-end pipeline for synthetic surfaces.
- `PC_Master_natural.m` — end-to-end pipeline for natural point clouds.

## Requirements

MATLAB with:
- Curve Fitting Toolbox (`fit`)
- Computer Vision Toolbox (`pointCloud`, `pcdownsample`, `pcread`)
- Statistics and Machine Learning Toolbox (`ksdensity`, `iqr`)

## Citation

If you use this code, please cite the associated manuscript (citation and Zenodo DOI to be added upon publication).

## License

MIT License — see [LICENSE](LICENSE).
