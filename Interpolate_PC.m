% Interpolate_PC.m - Creates a regularly spaced interpolated grid of Z
% values from a rotated, scattered point cloud. "samp" sets the FFT
% sampling rate; "bor*_m" (in meters) trims matrix edges affected by
% interpolation artifacts near the rotated boundary.
function [Rsnew3, d2, d3, F, XX, YY, Xq, Yq, qz, samp, sampFFT] = ...
    Interpolate_PC(Rsnew2, xyn, zn, borleft_m, borright_m, ...
    bordown_m, borup_m, filename, x0, r_deg, gridSize)

%% Estimate average spacing from bounding-box area / number of points
B = length(zn);
ca = 2*max(xyn(:,1));
fi = 2*max(xyn(:,2));
A = ca*fi;
samp = sqrt(A/B);

%% Regularly spaced query vectors
Xq = min(xyn(:,1)):samp:max(xyn(:,1));
Yq = min(xyn(:,2)):samp:max(xyn(:,2));

% Trim edges (artifacts occur where rotated axes do not meet unrotated
% ones; increase bor* or adjust r_deg if artifacts persist)
borleft = int32(borleft_m / samp);
borright = int32(borright_m / samp);
borup = int32(borup_m / samp);
bordown = int32(bordown_m / samp);
[Xqmin, Xqmax] = size(Xq);
Xqmin = Xqmin + borleft; Xqmax = Xqmax - borright;
[Yqmin, Yqmax] = size(Yq);
Yqmin = Yqmin + bordown; Yqmax = Yqmax - borup;
Xq = Xq(Xqmin:Xqmax);
Yq = Yq(Yqmin:Yqmax);
clear Xqmin Xqmax Yqmin Yqmax;

%% Interpolate onto regular grid
[XX, YY] = meshgrid(Xq,Yq);
F = scatteredInterpolant(xyn(:,1),xyn(:,2),zn,'linear','linear');
qz = F(XX,YY);
sampFFT = samp;

%% Figure: interpolated surface
figure(04), pcolor(XX,YY,qz), shading flat;
figure(4), imagesc(Xq, Yq, qz); axis xy;
set(04, 'Color',[1 1 1],'name','Interpolated point cloud');
set(gca,'fontweight','bold','fontsize',12,'fontname','Times New Roman');
title('Rotated data (interpolated)', 'FontWeight','bold','FontSize',16,'FontName','Times New Roman');
xlabel('x (m)', 'FontWeight','bold','FontSize',12,'FontName','Times New Roman');
ylabel('y (m)', 'FontWeight','bold','FontSize',12,'FontName','Times New Roman');
set(get(colorbar('EastOutside'),'ylabel'),'string','z (m)', ...
    'FontWeight','bold','FontSize',12,'FontName','Times New Roman');
axis equal;

%% Best-fit plane via SVD (post-cut, interpolated cloud)
Rsnew3 = [XX(:), YY(:), qz(:)];
x01 = mean(Rsnew3, 1);
Rsnew3 = Rsnew3 - x01;
[U3,S3,V3] = svd(Rsnew3, 0);
[s3, i3] = min(diag(S3));
a3 = V3(:, i3);
d3 = U3(:, i3)*s3;
if a3(3)<0
    d3 = -d3;
end
normd3 = norm(d3);
mean_3 = mean(d3);
median_3 = median(d3);
std_3 = std(d3);

figure(5),
Rsnew_pC3 = pointCloud(Rsnew3);
Rsnew_pCdown3 = pcdownsample(Rsnew_pC3, 'gridAverage', gridSize);
scatter(Rsnew_pCdown3.Location(:,1),Rsnew_pCdown3.Location(:,2), 2, Rsnew_pCdown3.Location(:,3))
set(5, 'name','RsnewpC3','Color',[1 1 1]);
set(gca,'fontweight','bold','fontsize',12,'fontname','Times New Roman');
title('RsnewpC3', 'FontWeight','bold','FontSize',16,'FontName','Times New Roman')
xlabel('x (m)', 'FontWeight','bold','FontSize',12,'FontName','Times New Roman')
ylabel('y (m)', 'FontWeight','bold','FontSize',12,'FontName','Times New Roman')
set(get(colorbar('EastOutside'),'ylabel'),'string','z (m)', ...
    'FontWeight','bold','FontSize',12,'FontName','Times New Roman');
grid off; ax.Layer = 'top';
xlim([min(Rsnew_pCdown3.Location(:,1)) max(Rsnew_pCdown3.Location(:,1))]);
ylim([min(Rsnew_pCdown3.Location(:,2)) max(Rsnew_pCdown3.Location(:,2))]);
axis equal;

%% Best-fit plane via SVD (pre-cut point cloud, for comparison)
[U2,S2,V2] = svd(Rsnew2,0);
[s2, i2] = min(diag(S2));
a2 = V2(:, i2);
d2 = U2(:, i2)*s2;
normd2 = norm(d2);
mean_2 = mean(d2);
median_2 = median(d2);
std_2 = std(d2);

%% Write interpolation results to file
nome1 = filename(1:end-4);
nome = nome1+"_resultsinterp.txt";
restext1 = string({'samp' 'sampFFT' 'r_deg' 'x0' 'y0' 'z0' ...
    'borleft_m' 'borright_m' 'bordown_m' 'borup_m'}');
restext2 = string({'normd2: Euclidean norm' 'median_2' 'mean_2' 'std_2'}');
restext3 = string({'normd3 : Euclidean norm' ...
    'median_3: Median value of residuals' ...
    'mean_3: Mean value of residuals' ...
    'std_3: Standard deviation of residuals'}');
resvars1 = [samp sampFFT r_deg x0(1) x0(2) x0(3) borleft_m borright_m bordown_m borup_m]';
resvars2 = [normd2 median_2 mean_2 std_2]';
resvars3 = [normd3 median_3 mean_3 std_3]';

fileID = fopen(nome,'w');
fprintf(fileID, '%s', 'name'); fprintf(fileID, ' \t'); fprintf(fileID, '%s \r\n', nome1);
for ci=1:numel(restext1)
    fprintf(fileID, '%s', restext1(ci)); fprintf(fileID, ' \t');
    fprintf(fileID, '%.7f \r\n', resvars1(ci));
end
fprintf(fileID, '\r\n');
fprintf(fileID, '%s ', 'Residual statistics for the POST-cut PC:'); fprintf(fileID, '\r\n');
for ci=1:numel(restext3)
    fprintf(fileID, '%s', restext3(ci)); fprintf(fileID, ' \t');
    fprintf(fileID, '%e \r\n', resvars3(ci));
end
fprintf(fileID, '\r\n');
fprintf(fileID, '%s ', 'Residual statistics for the PRE-cut PC:'); fprintf(fileID, '\r\n');
for ci=1:numel(restext2)
    fprintf(fileID, '%s', restext2(ci)); fprintf(fileID, ' \t');
    fprintf(fileID, '%e \r\n', resvars2(ci));
end
fclose(fileID);

%% Figure: residuals of post-cut point cloud
figure(6);
plot(d3, '.', 'MarkerSize', 1);
set(6, 'Color', [1 1 1], 'name', 'Residuals');
title('Residuals', 'FontWeight', 'bold', 'FontSize', 16, 'FontName', 'Times New Roman');
grid on;
yline(mean_3 + std_3, 'r--', 'LineWidth', 2);
yline(mean_3 - std_3, 'r--', 'LineWidth', 2);
xPos = length(d3) * 0.2;
yPos = max(d3) * 0.9;
text(xPos, yPos, sprintf('Standard Deviation: %.7f', std_3), ...
    'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Times New Roman', 'Color', 'black');

end
