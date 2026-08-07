function [absmeanx, xmean] = XFit_all_data(absXsqall, freqx_cut, filename)
% Fits the mean PSD along X using ALL raw spectral estimates (no binning).
% Derives H, alpha (prefactor), Rq_theoretical, and per-frequency
% dispersion (Std/IQR of log10 P(k)) used to assess data quality.

absmeanx = mean(absXsqall, 2);
if freqx_cut(1) == 0
    disp('Zero DC component detected and removed');
    xmean = [freqx_cut(2:end)', absmeanx(2:end)];
else
    xmean = [freqx_cut(1:end)', absmeanx(1:end)];
end
xmeanlog = log10(xmean);
if any(isinf(xmeanlog(:)))
    warning('Log of zero or near-zero values encountered.');
end

%% Fit
[fitobject,gof] = fit(xmeanlog(:,1), xmeanlog(:,2), 'poly1');
coefx = coeffvalues(fitobject);
Hx = (coefx(1) + 1) / (-2);
ax = 10^(coefx(2));
x1 = (10^(coefx(2))) * ((min(xmean(1:end,1)))^coefx(1));
x2 = (10^(coefx(2))) * ((max(xmean(1:end,1)))^coefx(1));
matx2 = [(min(xmean(1:end,1))) x1; (max(xmean(1:end,1))) x2];

%% Per-frequency density and dispersion (Std/IQR of log10 P(k))
xAll = repmat(freqx_cut(:), size(absXsqall,2), 1);
yAll = absXsqall(:);
validMask = xAll > 0 & yAll > 0 & isfinite(xAll) & isfinite(yAll);

nFreq = length(freqx_cut);
densityMat = nan(size(absXsqall));
freqDispersion = nan(nFreq, 4); % [frequency, std_log10P, IQR_log10P, nProfiles]
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

nome_disp = nome1 + "_x_dispersion_allmean.txt";
fileID2 = fopen(nome_disp, 'w');
fprintf(fileID2, 'Frequency_m-1\tStd_log10P\tIQR_log10P\tnProfiles\r\n');
fprintf(fileID2, '%.6f\t%.6f\t%.6f\t%d\r\n', freqDispersion');
fclose(fileID2);

%% Plot
figure(9);
scatter(xAll(validMask), yAll(validMask), 6, pointDensity(validMask), 'filled');
set(9, 'Color', [1 1 1], 'name', 'Power spectra along X');
colormap(gca, turbo);
clim(gca, [prctile(pointDensity(validMask), 2), prctile(pointDensity(validMask), 98)]);
cb = colorbar; cb.Label.String = 'Density along X (per frequency)'; cb.Label.FontWeight = 'bold';
set(gca, 'Color', [1 1 1], 'XColor', [0 0 0], 'YColor', [0 0 0]);
set(gca, 'fontweight', 'bold', 'fontsize', 12, 'fontname', 'Times New Roman');
set(gca, 'xscal', 'log'); set(gca, 'yscal', 'log');
hold on
plot(freqx_cut, absmeanx, 'ks-', 'MarkerFaceColor', 'k');
plot(matx2(:,1), matx2(:,2), 'r-', 'LineWidth', 2);
title('Power spectra along X', 'FontWeight', 'bold', 'FontSize', 16, 'FontName', 'Times New Roman', 'Color', [0 0 0]);
xlabel('Spatial frequency, k [m^-1]', 'FontWeight', 'bold', 'FontSize', 12, 'FontName', 'Times New Roman', 'Color', [0 0 0]);
ylabel('Power amplitude density, P(k) [m^3]', 'FontWeight', 'bold', 'FontSize', 12, 'FontName', 'Times New Roman', 'Color', [0 0 0]);
hold off;

%% Theoretical Rq (eq. 5)
fmin_fit = min(xmean(:,1)); fmax_fit = max(xmean(:,1));
Rq_theoretical = sqrt( ax/(2*Hx) * (fmin_fit^(-2*Hx) - fmax_fit^(-2*Hx)) );

%% Write results
nome = nome1 + "_x_allmean.txt";
dlmwrite(nome, xmean, '\t');
restext = string({'xcoef1', 'xcoef2', 'Hx', 'ax', 'Rq_theoretical', 'xgof.sse: Sum of Squares due to Error', ...
    'xgof.rsquare: R-squared', 'xgof.dfe: Degrees of Freedom for Error', ...
    'xgof.adjrsquare: Adjusted R-squared', 'xgof.rmse: Root Mean Squared Error'});
resvars = [coefx(1), coefx(2), Hx, ax, Rq_theoretical, gof.sse, gof.rsquare, gof.dfe, gof.adjrsquare, gof.rmse]';
nome = nome1 + "_x_results_allmean.txt";
fileID = fopen(nome,'w');
fprintf(fileID, '%s', 'name'); fprintf(fileID, ' \t'); fprintf(fileID, '%s \r\n', nome1);
for ci = 1:numel(restext)
    fprintf(fileID, '%s', restext(ci)); fprintf(fileID, '\t');
    if strcmp(restext(ci), 'ax')
        fprintf(fileID, '%.7e \r\n', resvars(ci));
    else
        fprintf(fileID, '%.7f \r\n', resvars(ci));
    end
end
fclose(fileID);
end
