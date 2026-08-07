% alignXY.m - Rotates a point cloud (already normal-aligned by alignZ)
% around Z by a user-specified angle, to match X/Y with the slip
% direction. A univocal x'/y' solution is required because alignZ only
% fixes z' (normal to the best-fit plane); this manual rotation avoids
% having to solve Wahba's problem for x'/y'.
function [Rsnew, Rsnew2, Rz, xyn, zn, r_deg] = alignXY(Rsnew, xyn, zn, r_deg, gridSize)

if ~exist('r_deg', 'var')
    r_deg = 0;
end
prompt = {'Enter clockwise rotation angle around z'};
r_deg2 = inputdlg(prompt, 'Rotation angle', [1 50], string(r_deg));
if isempty(r_deg2)
    uiwait(warndlg(sprintf('Empty input.\nAborted by user'))); return
else
    r_deg2 = str2double(r_deg2);
    if isnan(r_deg2)
        uiwait(warndlg('A value between -360 and +360 should be used.')); return
    else
        r_deg = r_deg2;
    end
end

%% Rotation
Rz = [cosd(r_deg) sind(r_deg) 0; -sind(r_deg) cosd(r_deg) 0; 0 0 1];
clear Rsnew2;
m = size(Rsnew, 1);
for i = 1:m
    Rsnew2(:,i) = Rz*(Rsnew(i,:))';
end
Rsnew2 = Rsnew2';
xyn = [Rsnew2(:,1), Rsnew2(:,2)];
zn = Rsnew2(:,3);

%% Visualize rotated and downsampled 2D grid
Rsnew_pC = pointCloud(Rsnew2);
Rsnew_pCdown = pcdownsample(Rsnew_pC, 'gridAverage', gridSize);
figure(03), scatter(Rsnew_pCdown.Location(:,1), Rsnew_pCdown.Location(:,2), ...
    2, Rsnew_pCdown.Location(:,3))
set(03, 'name','Rotated around z (downsampled)','Color',[1 1 1]);
set(gca,'fontweight','bold','fontsize',12,'fontname','Times New Roman');
title('Rotated around z', 'FontWeight','bold','FontSize',16,'FontName','Times New Roman')
xlabel('x (m)', 'FontWeight','bold','FontSize',12,'FontName','Times New Roman')
ylabel('y (m)', 'FontWeight','bold','FontSize',12,'FontName','Times New Roman')
set(get(colorbar('EastOutside'),'ylabel'),'string','z (m)', ...
    'FontWeight','bold','FontSize',12,'FontName','Times New Roman');
grid on; grid minor; ax.Layer = 'top';
xlim([min(Rsnew_pCdown.Location(:,1)) max(Rsnew_pCdown.Location(:,1))]);
ylim([min(Rsnew_pCdown.Location(:,2)) max(Rsnew_pCdown.Location(:,2))]);
axis equal;
end
