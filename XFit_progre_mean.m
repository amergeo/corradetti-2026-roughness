function [absmeanx, xmean] = XFit_progre_mean(absXsqall, freqx_cut, filename)
% Fits the mean PSD along X after geometric-progression binning
% (Corradetti et al. 2017): first 3 points kept raw, subsequent points
% averaged into log-spaced bins (step = log10(f4)-log10(f3); bin centre
% = geometric mean of bin edges, following Candela et al. 2012).

absmeanx = mean(absXsqall, 2);
if freqx_cut(1) == 0
    disp('Zero DC component detected and removed');
    xmean = [freqx_cut(2:end)', absmeanx(2:end)];
else
    xmean = [freqx_cut(:), absmeanx(:)];
end

%% Geometric-progression binning
nPoints = size(xmean, 1);
max_freq = max(xmean(:,1));
progmeanx = zeros(3, 2);
progmeanx(1:3, :) = xmean(1:3, :);
strep = log10(xmean(4, 1)) - log10(xmean(3, 1));
eligible = false(nPoints, 1);
eligible(4:end) = true;
xmeanlog1 = log10(xmean(:, 1));
binIndex = 3;
t = log10(xmean(3, 1));
while t < log10(max_freq)
    binIndex = binIndex + 1;
    inBin = eligible & (xmeanlog1 >= t) & (xmeanlog1 < t + strep);
    bin_rows = xmean(inBin, :);
    if size(bin_rows, 1) == 1
        progmeanx(binIndex, :) = bin_rows(1, :);
    elseif size(bin_rows, 1) > 1
        bin_min = 10^t; bin_max = 10^(t + strep);
        progmeanx(binIndex, 1) = sqrt(bin_min * bin_max);
        progmeanx(binIndex, 2) = mean(bin_rows(:, 2));
    else
        progmeanx(binIndex, :) = [NaN, NaN];
    end
    t = t + strep;
end
xprogmean = progmeanx(~isnan(progmeanx(:, 2)), :);

%% Fit
[fitobject, gof] = fit(log10(xprogmean(:, 1)), log10(xprogmean(:, 2)), 'poly1');
coefx = coeffvalues(fitobject);
Hx = (coefx(1) + 1) / (-2);
ax = 10^(coefx(2));
x1 = (10^(coefx(2))) * ((min(xprogmean(:, 1)))^coefx(1));
x2 = (10^(coefx(2))) * ((max(xprogmean(:, 1)))^coefx(1));
matx = [min(xprogmean(:, 1)), x1; max(xprogmean(:, 1)), x2];

%% Per-frequency density and dispersion
xAll = repmat(freqx_cut(:), size(absXsqall,2), 1);
yAll = absXsqall(:);
validMask = xAll > 0 & yAll > 0 & isfinite(xAll) & isfinite(yAll);
nFreq = length(freqx_cut);
densityMat = nan(size(absXsqall));
freqDispersion = nan(nFreq, 4);
nome1 = filename(1:end-4);
for fi = 1:nFreq
    rowVals = absXsqall(fi, :);
    rowMask = rowVals > 0 & isfinite(rowVals);
    freqDispersion(fi,1) = freqx_cut(fi);
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
nome_disp = nome1 + "_x_dispersion_progmean.txt";
fileID2 = fopen(nome_disp, 'w');
fprintf(fileID2, 'Frequency_m-1\tStd_log10P\tIQR_log10P\tnProfiles\r\n');
fprintf(fileID2, '%.6f\t%.6f\t%.6f\t%d\r\n', freqDispersion');
fclose(fileID2);

%% Plot
figure(10);
scatter(xAll(validMask), yAll(validMask), 6, pointDensity(validMask), 'filled'); hold on;
plot(xprogmean(:,1), xprogmean(:,2), 'ks-', 'MarkerFaceColor', 'k', 'LineWidth', 1.5);
plot(matx(:,1), matx(:,2), 'r-', 'LineWidth', 2);
colormap(gca, turbo);
clim(gca, [prctile(pointDensity(validMask), 2), prctile(pointDensity(validMask), 98)]);
cb = colorbar; cb.Label.String = 'Density along X (per frequency)'; cb.Label.FontWeight = 'bold';
set(gca, 'XScale', 'log', 'YScale', 'log');
set(gca, 'Color', [1 1 1], 'XColor', [0 0 0], 'YColor', [0 0 0]);
set(gca, 'FontWeight', 'bold', 'FontSize', 12, 'FontName', 'Times New Roman');
set(gcf, 'Color', [1 1 1], 'Name', 'Power spectra along X');
title('Power spectra along X', 'FontWeight', 'bold', 'FontSize', 16, 'FontName', 'Times New Roman', 'Color', [0 0 0]);
xlabel('Spatial frequency, k [m^{-1}]', 'FontWeight', 'bold', 'FontSize', 12, 'FontName', 'Times New Roman', 'Color', [0 0 0]);
ylabel('Power amplitude density, P(k) [m^3]', 'FontWeight', 'bold', 'FontSize', 12, 'FontName', 'Times New Roman', 'Color', [0 0 0]);
hold off;

%% Theoretical Rq (eq. 5)
fmin_fit = min(xprogmean(:,1)); fmax_fit = max(xprogmean(:,1));
Rq_theoretical = sqrt( ax/(2*Hx) * (fmin_fit^(-2*Hx) - fmax_fit^(-2*Hx)) );

%% Write results
nome_data = nome1 + "_x_progmean.txt";
dlmwrite(nome_data, xprogmean, '\t');
restext = string({'xcoef1', 'xcoef2', 'Hx', 'ax', 'Rq_theoretical', ...
    'xgof.sse', 'xgof.rsquare', 'xgof.dfe', 'xgof.adjrsquare', 'xgof.rmse'});
resvars = [coefx(1), coefx(2), Hx, ax, Rq_theoretical, gof.sse, gof.rsquare, gof.dfe, gof.adjrsquare, gof.rmse]';
nome_results = nome1 + "_x_results_progmean.txt";
fileID = fopen(nome_results, 'w');
fprintf(fileID, 'name\t%s\r\n', nome1);
for ci = 1:numel(restext)
    fprintf(fileID, '%s\t', restext(ci));
    if strcmp(restext(ci), 'ax')
        fprintf(fileID, '%.7e\r\n', resvars(ci));
    else
        fprintf(fileID, '%.7f\r\n', resvars(ci));
    end
end
fclose(fileID);
end
