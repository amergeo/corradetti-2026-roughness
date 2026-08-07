% PC_Master_natural.m
% Pipeline for natural point clouds (e.g. SfM-MVS models) requiring
% alignment, rotation and interpolation onto a regular grid before FFT.
clear all; close all; clc;

%% Step 1: Import point cloud data
[Raw_surface, filename] = import_PC();

%% Step 2: Align Z to the least-squares best-fit plane
[Rs, Rsnew, x0, a, d, normd, R, xyn, zn, r_deg] = alignZ(Raw_surface);
[gridSize] = Figura_raw(Rsnew);

%% Step 3: Rotate around Z (align X/Y to the slip direction)
[Rsnew, Rsnew2, Rz, xyn, zn, r_deg] = alignXY(Rsnew, xyn, zn, r_deg, gridSize);

%% Step 4: Interpolate onto a regular grid
borleft_m = 0.01; borright_m = 0.01; borup_m = 0.01; bordown_m = 0.01;
[Rsnew3, d2, d3, F, XX, YY, Xq, Yq, qz, samp, sampFFT] = ...
    Interpolate_PC(Rsnew2, xyn, zn, borleft_m, borright_m, ...
    bordown_m, borup_m, filename, x0, r_deg, gridSize);

%% Step 4b (optional): crop to a square extent before FFT
% make_square_surface

%% Step 5: FFT along Y (slip-parallel) and X (slip-perpendicular)
[Yx, freqy, absYsqall, freqy_cut] = processFFT_YAxis(qz, sampFFT);
[Xx, freqx, absXsqall, freqx_cut] = processFFT_XAxis(qz, sampFFT);

%% Step 6: Fit along Y
[absmeany, ymean]            = YFit_all_data(absYsqall, freqy_cut, filename);
[absmeany, ymean]            = YFit_progre_mean(absYsqall, freqy_cut, filename);
[absmeany, ymean, ymean_cut] = YFit_userbounded_progmean(absYsqall, freqy_cut, filename);

%% Step 7: Fit along X
[absmeanx, xmean]            = XFit_all_data(absXsqall, freqx_cut, filename);
[absmeanx, xmean]            = XFit_progre_mean(absXsqall, freqx_cut, filename);
[absmeanx, xmean, xmean_cut] = XFit_userbounded_progmean(absXsqall, freqx_cut, filename);
