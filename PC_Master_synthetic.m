% PC_Master_synthetic.m
% Pipeline for surfaces already on a regular grid with Y = slip direction
% (synthetic surfaces from generateSelfAffineSurface.m, optionally passed
% through noise_developer.m). Skips alignZ/alignXY/Interpolate_PC.
clear all; close all; clc;

%% Step 1: Import
[Raw_surface, filename] = import_PC();

%% Step 2: Direct FFT (no alignment/interpolation needed)
[qz, sampFFT, freqx, absXsqall, freqx_cut, freqy, absYsqall, freqy_cut] ...
    = PrepareSurfaceForFFT(Raw_surface);

%% Step 3: Fit along Y (slip-parallel)
[absmeany, ymean]                     = YFit_all_data(absYsqall, freqy_cut, filename);
[absmeany, ymean]                     = YFit_progre_mean(absYsqall, freqy_cut, filename);
[absmeany, ymean, ymean_cut]          = YFit_userbounded_progmean(absYsqall, freqy_cut, filename);

%% Step 4: Fit along X (slip-perpendicular)
[absmeanx, xmean]                     = XFit_all_data(absXsqall, freqx_cut, filename);
[absmeanx, xmean]                     = XFit_progre_mean(absXsqall, freqx_cut, filename);
[absmeanx, xmean, xmean_cut]          = XFit_userbounded_progmean(absXsqall, freqx_cut, filename);
