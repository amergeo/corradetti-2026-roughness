function noise_developer()
% Adds bandlimited high-frequency noise to a Run_XXXX.xyz surface file,
% simulating instrumental (photogrammetric/LiDAR) noise dominating at
% wavelengths shorter than a user-defined cutoff lambda_c.
%
% Noise is generated as white noise in the Fourier domain, then
% high-pass filtered (isotropic radial cutoff, k > k_c = 1/lambda_c)
% and scaled to a target RMS amplitude Rq_noise.
%
% Output: Run_noise_XXXX.xyz saved alongside the input file.

%% 1. Select input file
[fname, fpath] = uigetfile('*.xyz', 'Select a Run_XXXX.xyz surface file');
if isequal(fname, 0)
    disp('Cancelled by user.'); return;
end
fullpath = fullfile(fpath, fname);

%% 2. Parameter GUI
prompt = {'Noise cutoff wavelength lambda_c (m)  [noise added at lambda < lambda_c]', ...
    'Noise RMS amplitude  Rq_noise (m)', 'Random seed (integer)'};
def_ans = {'0.1', '1e-4', '42'};
answ = inputdlg(prompt, sprintf('Noise Developer - %s', fname), [1 68], def_ans);
if isempty(answ)
    disp('Cancelled by user.'); return;
end
lambda_c = str2double(answ{1});
Rq_noise = str2double(answ{2});
noise_seed = round(str2double(answ{3}));
if lambda_c <= 0; error('lambda_c must be > 0 m.'); end
if Rq_noise <= 0; error('Rq_noise must be > 0 m.'); end

%% 3. Load surface and reconstruct regular grid
fprintf('>>> Loading %s ...\n', fname);
raw = load(fullpath);
if size(raw, 2) < 3
    error('Expected 3 columns (x y z) in the XYZ file.');
end
x_col = raw(:, 1); y_col = raw(:, 2); z_col = raw(:, 3);
xu = unique(x_col, 'sorted'); yu = unique(y_col, 'sorted');
nx = numel(xu); ny = numel(yu);
if nx * ny ~= numel(z_col)
    error('Cannot reconstruct a regular grid: nx*ny (%d) ~= npoints (%d).', nx*ny, numel(z_col));
end
Z = zeros(ny, nx);
for k = 1:numel(z_col)
    ix = find(xu == x_col(k), 1);
    iy = find(yu == y_col(k), 1);
    Z(iy, ix) = z_col(k);
end
Lx = xu(end) - xu(1); Ly = yu(end) - yu(1);
dx = mean(diff(xu)); dy = mean(diff(yu));
fprintf('    Grid: %d x %d  |  Lx=%.4g m  Ly=%.4g m  dx=%.4g m  dy=%.4g m\n', nx, ny, Lx, Ly, dx, dy);
if lambda_c <= max(dx, dy)
    warning('lambda_c (%.4g m) <= grid spacing (%.4g m). No noise will be added.', lambda_c, max(dx, dy));
end

%% 4. Generate isotropic bandlimited noise (high-pass, k > k_c)
rng(noise_seed, 'twister');
fx = ifftshift((-floor(nx/2):ceil(nx/2)-1) / Lx);
fy = ifftshift((-floor(ny/2):ceil(ny/2)-1) / Ly);
[FX, FY] = meshgrid(fx, fy);
K = sqrt(FX.^2 + FY.^2);
k_c = 1 / lambda_c;
G = randn(ny, nx) + 1i * randn(ny, nx);
mask = double(K > k_c);
mask(K == 0) = 0;   % remove DC
Znoise_raw = ifft2(G .* mask, 'symmetric');
Znoise_raw = double(Znoise_raw - mean(Znoise_raw(:)));
current_Rq = sqrt(mean(Znoise_raw(:).^2));
if current_Rq == 0
    error('Noise RMS is zero after filtering. Try a different seed or lambda_c.');
end
Znoise = Znoise_raw * (Rq_noise / current_Rq);

%% 5. Add noise and compute statistics
Z_noisy = Z + Znoise;
Rq_orig = sqrt(mean(Z(:).^2));
Rq_final = sqrt(mean(Z_noisy(:).^2));
SNR_dB = 20 * log10(Rq_orig / Rq_noise);
fprintf('    Original  Rq = %.4g m\n', Rq_orig);
fprintf('    Noise     Rq = %.4g m  (target %.4g m)\n', sqrt(mean(Znoise(:).^2)), Rq_noise);
fprintf('    Noisy     Rq = %.4g m\n', Rq_final);
fprintf('    SNR          = %.1f dB\n', SNR_dB);
fprintf('    lambda_c     = %.4g m  |  k_c = %.4g m^-1\n', lambda_c, k_c);

%% 6. Save output file Run_noise_XXXX.xyz
tok = regexp(fname, 'Run_(\d+)\.xyz', 'tokens', 'ignorecase');
if isempty(tok)
    error('Input filename does not match expected pattern "Run_XXXX.xyz".');
end
runNum = tok{1}{1};
outName = sprintf('Run_noise_%s.xyz', runNum);
outPath = fullfile(fpath, outName);
[X, Y] = meshgrid(xu, yu);
fid = fopen(outPath, 'w');
if fid == -1; error('Cannot write output file "%s".', outPath); end
fprintf(fid, '%.8f\t%.8f\t%.8f\n', [X(:), Y(:), Z_noisy(:)]');
fclose(fid);
fprintf('>>> Saved: %s\n', outPath);

%% 7. Figures
clim_surf = [min(Z(:)), max(Z(:))];

figure(1); clf;
imagesc(xu, yu, Z); set(gca, 'YDir', 'normal'); axis image; clim(clim_surf);
xlabel('X (m) - perpendicular to slip'); ylabel('Y (m) - along slip');
title(sprintf('Original surface: %s', fname), 'Interpreter', 'none');
colorbar; colormap(parula);

figure(2); clf;
imagesc(xu, yu, Znoise); set(gca, 'YDir', 'normal'); axis image;
xlabel('X (m)'); ylabel('Y (m)');
title(sprintf('Added noise  |  \\lambda_c = %.4g m  |  Rq_{noise} = %.4g m', lambda_c, Rq_noise));
colorbar; colormap(parula);

figure(3); clf;
imagesc(xu, yu, Z_noisy); set(gca, 'YDir', 'normal'); axis image; clim(clim_surf);
xlabel('X (m) - perpendicular to slip'); ylabel('Y (m) - along slip');
title(sprintf('Noisy surface: %s  |  SNR = %.1f dB', outName, SNR_dB), 'Interpreter', 'none');
colorbar; colormap(parula);

end
