function [absmeany, ymean, ymean_cut] = YFit_userbounded_progmean(absYsqall, freqy_cut, filename)
% Same geometric-progression binning as YFit_progre_mean.m, restricted to
% a user-defined wavelength/frequency window (via freqRangeGUI). Mirrors
% XFit_userbounded_progmean.m.

absmeany = mean(absYsqall, 2);
if freqy_cut(1) == 0
    disp('Zero DC component detected and removed');
    ymean = [freqy_cut(2:end)', absmeany(2:end)];
else
    ymean = [freqy_cut(:), absmeany(:)];
end

%% User-selected frequency range
min_freq = min(ymean(:,1));
max_freq = max(ymean(:,1));
[min_freq_user, max_freq_user, ~] = freqRangeGUI(min_freq, max_freq, false);
disp(['User selected Min Frequency: ', num2str(min_freq_user)]);
disp(['User selected Max Frequency: ', num2str(max_freq_user)]);
ymean_cut = ymean(ymean(:, 1) >= min_freq_user & ymean(:, 1) <= max_freq_user, :);

%% Geometric-progression binning within the selected range
nPoints = size(ymean_cut, 1);
max_freq_cut = max(ymean_cut(:,1));
progmeany = zeros(3, 2);
progmeany(1:3, :) = ymean_cut(1:3, :);
strep = log10(ymean_cut(4, 1)) - log10(ymean_cut(3, 1));
eligible = false(nPoints, 1);
eligible(4:end) = true;
ymeanlog1 = log10(ymean_cut(:, 1));
binIndex = 3;
t = log10(ymean_cut(3, 1));
while t < log10(max_freq_cut)
    binIndex = binIndex + 1;
    inBin = eligible & (ymeanlog1 >= t) & (ymeanlog1 < t + strep);
    bin_rows = ymean_cut(inBin, :);
    if size(bin_rows, 1) == 1
        progmeany(binIndex, :) = bin_rows(1, :);
    elseif size(bin_rows, 1) > 1
        bin_min = 10^t; bin_max = 10^(t + strep);
        progmeany(binIndex, 1) = sqrt(bin_min * bin_max);
        progmeany(binIndex, 2) = mean(bin_rows(:, 2));
    else
        progmeany(binIndex, :) = [NaN, NaN];
    end
    t = t + strep;
end
yprogmean = progmeany(~isnan(progmeany(:, 2)), :);

%% Fit
[fitobject, gof] = fit(log10(yprogmean(:, 1)), log10(yprogmean(:, 2)), 'poly1');
coefy = coeffvalues(fitobject);
Hy = (coefy(1) + 1) / (-2);
ay = 10^(coefy(2));
y1 = (10^(coefy(2))) * ((min(yprogmean(:, 1)))^coefy(1));
y2 = (10^(coefy(2))) * ((max(yprogmean(:, 1)))^coefy(1));
maty = [min(yprogmean(:, 1)), y1; max(yprogmean(:, 1)), y2];

%% Per-frequency density and dispersion
xAll = repmat(freqy_cut(:), size(absYsqall,2), 1);
yAll = absYsqall(:);
validMask = xAll > 0 & yAll > 0 & isfinite(xAll) & isfinite(yAll);
nFreq = length(freqy_cut);
densityMat = nan(size(absYsqall));
freqDispersion = nan(nFreq, 4);
nome1 = filename(1:end-4);
for fi = 1:nFreq
    rowVals = absYsqall(fi, :);
    rowMask = rowVals > 0 & isfinite(rowVals);
    freqDispersion(fi,1) = freqy_cut(fi);
    if sum(rowMask) > 1
        logRow = log10(rowVals(rowMask));
        densityMat(fi, rowMask) = ksdensity(logRow, logRow);
        freqDispersion(fi,2) = std(logRow);
        freqDispersion(fi,3) = iqr(logRow);
        freqDispersion(fi,4) = sum(rowMask);
    elseif sum(rowMask) == 1
        densityMat(fi, rowMask) = 1;
        freqDispersion(fi,2:4) = [0, 0, 1];
    end
end
pointDensity = densityMat(:);
nome_disp = nome1 + "_y_dispersion_userbounded_progmean.txt";
fileID2 = fopen(nome_disp, 'w');
fprintf(fileID2, 'Frequency_m-1\tStd_log10P\tIQR_log10P\tnProfiles\r\n');
fprintf(fileID2, '%.6f\t%.6f\t%.6f\t%d\r\n', freqDispersion');
fclose(fileID2);

%% Plot
figure(12);
scatter(xAll(validMask), yAll(validMask), 6, pointDensity(validMask), 'filled'); hold on;
plot(yprogmean(:,1), yprogmean(:,2), 'ks-', 'MarkerFaceColor', 'k', 'LineWidth', 1.5);
plot(maty(:,1), maty(:,2), 'r-', 'LineWidth', 2);
colormap(gca, turbo);
clim(gca, [prctile(pointDensity(validMask), 2), prctile(pointDensity(validMask), 98)]);
cb = colorbar; cb.Label.String = 'Density along Y (per frequency)'; cb.Label.FontWeight = 'bold';
set(gca, 'XScale', 'log', 'YScale', 'log');
set(gca, 'Color', [1 1 1], 'XColor', [0 0 0], 'YColor', [0 0 0]);
set(gca, 'FontWeight', 'bold', 'FontSize', 12, 'FontName', 'Times New Roman');
set(gcf, 'Color', [1 1 1], 'Name', 'Power spectra along Y');
title('Power spectra along Y', 'FontWeight', 'bold', 'FontSize', 16, 'FontName', 'Times New Roman', 'Color', [0 0 0]);
xlabel('Spatial frequency, k [m^{-1}]', 'FontWeight', 'bold', 'FontSize', 12, 'FontName', 'Times New Roman', 'Color', [0 0 0]);
ylabel('Power amplitude density, P(k) [m^3]', 'FontWeight', 'bold', 'FontSize', 12, 'FontName', 'Times New Roman', 'Color', [0 0 0]);
hold off;

%% Theoretical Rq (eq. 5)
fmin_fit = min(yprogmean(:,1)); fmax_fit = max(yprogmean(:,1));
Rq_theoretical = sqrt( ay/(2*Hy) * (fmin_fit^(-2*Hy) - fmax_fit^(-2*Hy)) );

%% Write results
nome_data = nome1 + "_y_userbounded_progmean_.txt";
dlmwrite(nome_data, yprogmean, '\t');
restext = string({'ycoef1', 'ycoef2', 'Hy', 'ay', 'Rq_theoretical', ...
    'ygof.sse', 'ygof.rsquare', 'ygof.dfe', 'ygof.adjrsquare', 'ygof.rmse'});
resvars = [coefy(1), coefy(2), Hy, ay, Rq_theoretical, gof.sse, gof.rsquare, gof.dfe, gof.adjrsquare, gof.rmse]';
nome_results = nome1 + "_y_results_userbounded_progmean_.txt";
fileID = fopen(nome_results, 'w');
fprintf(fileID, 'name\t%s\r\n', nome1);
for ci = 1:numel(restext)
    fprintf(fileID, '%s\t', restext(ci));
    if strcmp(restext(ci), 'ay')
        fprintf(fileID, '%.7e\r\n', resvars(ci));
    else
        fprintf(fileID, '%.7f\r\n', resvars(ci));
    end
end
fclose(fileID);
end
