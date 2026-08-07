function [Yx, freqy, absYsqall, freqy_cut] = processFFT_YAxis(qz, sampFFT)
% Computes the 1D power spectral density P(k) [m^3] of every profile
% along Y. Same taper, correction, and PSD normalization as
% processFFT_XAxis.m (see Methods for the density-normalization note).

Fs = 1 / sampFFT;
Ny = length(qz(:, 1));
ky = 0:Ny-1;
Ty = Ny / Fs;
freqy = ky / Ty;

%% 3% taper window
taper_length = round(0.03 * Ny);
taper_window = ones(Ny, 1);
taper_window(1:taper_length) = linspace(0, 1, taper_length);
taper_window(end-taper_length+1:end) = linspace(1, 0, taper_length);
correction_factor = 1 / sqrt(mean(taper_window.^2));   % RMS-based amplitude correction

%% FFT loop over all profiles
for indy = 1:length(qz(1, :))
    y = qz(:, indy);
    yw = y(:) .* taper_window;
    Yx = fft(yw) / Ny;

    if rem(Ny, 2) == 0
        cutOff = Ny / 2 + 1;
    else
        cutOff = (Ny + 1) / 2;
    end
    Yx_cut = Yx(1:cutOff);
    freqy_cut = freqy(1:cutOff);

    absY = abs(Yx_cut);
    if rem(Ny, 2) == 0
        absY(2:end-1) = absY(2:end-1) * 2;
    else
        absY(2:end) = absY(2:end) * 2;
    end
    absY = absY * correction_factor;
    absYsq = absY.^2;
    absYsq = absYsq * (Ny * sampFFT);   % power-per-bin -> PSD [m^3]

    if indy == 1
        absYsqall = zeros(length(absYsq), length(qz(1, :)));
        absYsqall(:, 1) = absYsq(:);
    else
        absYsqall(:, indy) = absYsq(:);
    end
end

%% Plot representative profiles
numProfiles = size(qz, 2);
profileIdx1 = round(numProfiles / 3);
profileIdx2 = round(2 * numProfiles / 3);
figure(52);
plot(1:Ny, qz(:, profileIdx1), 'b-', 'LineWidth', 1.5); hold on;
plot(1:Ny, qz(:, profileIdx2), 'g-', 'LineWidth', 1.5);
grid on;
title('Selected Profiles Along Y', 'FontWeight', 'bold', 'FontSize', 16);
xlabel('Index Along Y', 'FontWeight', 'bold', 'FontSize', 12);
ylabel('Amplitude', 'FontWeight', 'bold', 'FontSize', 12);
legend(['Profile at 1/3 (Index ', num2str(profileIdx1), ')'], ...
       ['Profile at 2/3 (Index ', num2str(profileIdx2), ')'], 'Location', 'best');
hold off;
end
