function [absmeany, ymean] = YFit_all_data(absYsqall, freqy_cut, filename)
% Fits the mean PSD along Y using ALL raw spectral estimates (no binning).
% Mirrors XFit_all_data.m.

absmeany = mean(absYsqall,2);
if freqy_cut(1) == 0
    disp('Zero DC component detected and removed');
    ymean = [freqy_cut(2:end)', absmeany(2:end)];
else
    ymean = [freqy_cut(1:end)', absmeany(1:end)];
end
ymeanlog = log10(ymean);
if any(isinf(ymeanlog(:)))
    warning('Log of zero or near-zero values encountered.');
end

%% Fit
[fitobject,gof] = fit(ymeanlog(:,1), ymeanlog(:,2),'poly1');
coefy = coeffvalues(fitobject);
Hy = (coefy(1)+1)/(-2);
ay = 10^(coefy(2));
y1 = (10^(coefy(2)))*((min(ymean(1:end,1)))^coefy(1));
y2 = (10^(coefy(2)))*((max(ymean(1:end,1)))^coefy(1));
maty2 = [(min(ymean(1:end,1))) y1; (max(ymean(1:end,1))) y2];

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
nome_disp = nome1 + "_y_dispersion_allmean.txt";
fileID2 = fopen(nome_disp, 'w');
fprintf(fileID2, 'Frequency_m-1\tStd_log10P\tIQR_log10P\tnProfiles\r\n');
fprintf(fileID2, '%.6f\t%.6f\t%.6f\t%d\r\n', freqDispersion');
fclose(fileID2);

%% Plot
figure(7);
scatter(xAll(validMask), yAll(validMask), 6, pointDensity(validMask), 'filled');
set(7, 'Color',[1 1 1],'name','Power spectra along Y');
colormap(gca, turbo);
clim(gca, [prctile(pointDensity(validMask), 2), prctile(pointDensity(validMask), 98)]);
cb = colorbar; cb.Label.String = 'Density along Y (per frequency)'; cb.Label.FontWeight = 'bold';
set(gca,'Color',[1 1 1],'XColor',[0 0 0],'YColor',[0 0 0]);
set(gca,'fontweight','bold','fontsize',12,'fontname','Times New Roman');
set(gca,'xscal','log'); set(gca,'yscal','log');
hold on
plot(freqy_cut,absmeany,'ks-','MarkerFaceColor','k');
plot(maty2(:,1),maty2(:,2),'r-', 'LineWidth', 2);
title('Power spectra along Y', 'FontWeight','bold','FontSize',16,'FontName','Times New Roman','Color',[0 0 0])
xlabel('Spatial frequency, k [m^-1]', 'FontWeight','bold','FontSize',12,'FontName','Times New Roman','Color',[0 0 0])
ylabel('Power amplitude density, P(k) [m^3]', 'FontWeight','bold','FontSize',12,'FontName','Times New Roman','Color',[0 0 0])
hold off

%% Theoretical Rq (eq. 5)
fmin_fit = min(ymean(:,1)); fmax_fit = max(ymean(:,1));
Rq_theoretical = sqrt( ay/(2*Hy) * (fmin_fit^(-2*Hy) - fmax_fit^(-2*Hy)) );

%% Write results
nome1s = filename(1:end-4);
nome = nome1s+"_y_allmean_.txt";
dlmwrite(nome,ymean,'\t');
restext = string({'ycoef1' 'ycoef2' 'Hy' 'ay' 'Rq_theoretical' 'ygof.sse: Sum of Squares due to Error'...
    'ygof.rsquare: R-squared' 'ygof.dfe: Degrees of Freedom for Error' 'ygof.adjrsquare: Adjusted R-squared'...
    'ygof.rmse: Root Mean Squared Error'}');
resvars = [coefy(1) coefy(2) Hy ay Rq_theoretical gof.sse gof.rsquare gof.dfe gof.adjrsquare gof.rmse]';
nome = nome1s+"_y_results_allmean.txt";
fileID = fopen(nome,'w');
fprintf(fileID, '%s', 'name'); fprintf(fileID, ' \t'); fprintf(fileID, '%s \r\n', nome1s);
for ci = 1:numel(restext)
    fprintf(fileID, '%s', restext(ci)); fprintf(fileID, '\t');
    if strcmp(restext(ci), 'ay')
        fprintf(fileID, '%.7e \r\n', resvars(ci));
    else
        fprintf(fileID, '%.7f \r\n', resvars(ci));
    end
end
fclose(fileID);
end
