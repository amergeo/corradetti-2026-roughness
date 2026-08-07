function [qz, sampFFT, freqx, absXsqall, freqx_cut, freqy, absYsqall, freqy_cut] = PrepareSurfaceForFFT(Raw_surface)
% PrepareSurfaceForFFT - Direct path for surfaces already on a regular
% grid (e.g. synthetic surfaces), bypassing alignZ/alignXY/Interpolate_PC.
% Reshapes Raw_surface into a Z grid and runs the FFT analysis directly.

Rs = table2array(Raw_surface);
x_coords = double(Rs(:, 1));
y_coords = double(Rs(:, 2));
z_coords = double(Rs(:, 3));

if mod(length(unique(x_coords)), length(unique(y_coords))) ~= 0
    error('The point cloud is not regularly spaced. Check the input data.');
end

numX = length(unique(x_coords));
numY = length(unique(y_coords));
qz = reshape(z_coords, numY, numX);

x_unique = unique(x_coords, 'sorted');
y_unique = unique(y_coords, 'sorted');
sampX = mean(diff(x_unique));
sampY = mean(diff(y_unique));
if abs(sampX - sampY) > 1e-6
    error('Inconsistent grid spacing detected between X and Y.');
end
sampFFT = sampX;

[Xx, freqx, absXsqall, freqx_cut] = processFFT_XAxis(qz, sampFFT);
[Yx, freqy, absYsqall, freqy_cut] = processFFT_YAxis(qz, sampFFT);

save('FFT_Results_Direct.mat', 'qz', 'sampFFT', 'freqx', 'absXsqall', 'freqx_cut', ...
    'freqy', 'absYsqall', 'freqy_cut');

figure(31);
imagesc(x_unique, y_unique, qz); axis xy;
title('Surface', 'FontWeight', 'bold', 'FontSize', 16);
xlabel('X (m)', 'FontWeight', 'bold', 'FontSize', 12);
ylabel('Y (m)', 'FontWeight', 'bold', 'FontSize', 12);
colorbar; axis equal;
end
