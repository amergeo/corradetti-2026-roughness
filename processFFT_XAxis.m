function [Xx, freqx, absXsqall, freqx_cut] = processFFT_XAxis(qz, sampFFT)
% Computes the 1D power spectral density P(k) [m^3] of every profile
% along X. A 3% linear taper is applied at both ends (reduces spectral
% leakage) and corrected by dividing by the RMS of the taper window.
% Power is normalized by the frequency resolution df = 1/(Nx*sampFFT) to
% give a proper one-sided PSD (not the Schuster/power-per-bin convention
% used in Corradetti et al. 2017); H is unaffected by this rescaling,
% only the intercept (prefactor a) changes. See Methods for details.

Fs = 1 / sampFFT;
Nx = length(qz(1, :));
kx = 0:Nx-1;
Tx = Nx / Fs;
freqx = kx / Tx;

%% 3% taper window
taper_length = round(0.03 * Nx);
taper_window = ones(Nx, 1);
taper_window(1:taper_length) = linspace(0, 1, taper_length);
taper_window(end-taper_length+1:end) = linspace(1, 0, taper_length);
correction_factor = 1 / sqrt(mean(taper_window.^2));   % RMS-based amplitude correction

%% FFT loop over all profiles
for indx = 1:length(qz(:, 1))
    x = qz(indx, :);
    xw = x(:) .* taper_window;
    Xx = fft(xw) / Nx;

    if rem(Nx, 2) == 0
        cutOff = Nx / 2 + 1;
    else
        cutOff = (Nx + 1) / 2;
    end
    Xx_cut = Xx(1:cutOff);
    freqx_cut = freqx(1:cutOff);

    absX = abs(Xx_cut);
    if rem(Nx, 2) == 0
        absX(2:end-1) = absX(2:end-1) * 2;
    else
        absX(2:end) = absX(2:end) * 2;
    end
    absX = absX * correction_factor;
    absXsq = absX.^2;
    absXsq = absXsq * (Nx * sampFFT);   % power-per-bin -> PSD [m^3]

    if indx == 1
        absXsqall = zeros(length(absXsq), length(qz(:, 1)));
        absXsqall(:, 1) = absXsq(:);
    else
        absXsqall(:, indx) = absXsq(:);
    end
end

%% Plot representative profiles
numProfiles = size(qz, 1);
profileIdx1 = round(numProfiles / 3);
profileIdx2 = round(2 * numProfiles / 3);
figure(51);
plot(1:Nx, qz(profileIdx1, :), 'b-', 'LineWidth', 2); hold on;
plot(1:Nx, qz(profileIdx2, :), 'g-', 'LineWidth', 1.5);
grid on;
title('Selected Profiles Along X', 'FontWeight', 'bold', 'FontSize', 16);
xlabel('Index Along X', 'FontWeight', 'bold', 'FontSize', 12);
ylabel('Amplitude', 'FontWeight', 'bold', 'FontSize', 12);
legend(['Profile at 1/3 (Index ', num2str(profileIdx1), ')'], ...
       ['Profile at 2/3 (Index ', num2str(profileIdx2), ')'], 'Location', 'best');
hold off;
end
