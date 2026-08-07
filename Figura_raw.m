function [gridSize] = Figura_raw(Rsnew)
% Asks the user for a downsampling grid size (visualization only, does
% not affect the point cloud file) and plots the raw, Z-aligned point
% cloud in 3D and as a 2D top-view scatter with orientation markers.
gridSize = NaN;

fig = uifigure('Name', 'Set Grid Size', 'Position', [400, 300, 300, 200]);
uitextarea(fig, 'Position', [20 120 260 80], 'Value', ...
    {'Enter a grid size for downsampling.', 'This will not affect the point cloud file.', ...
    'Value in meters (default: 0.001):'}, 'Editable', 'off');
gridSizeInput = uieditfield(fig, 'numeric', 'Position', [100 60 100 22], 'Value', 0.001);
uibutton(fig, 'Position', [100 10 100 30], 'Text', 'Proceed', ...
    'ButtonPushedFcn', @(btn, event) confirmGridSize());

    function confirmGridSize()
        gridSize = gridSizeInput.Value;
        close(fig);
    end

waitfor(fig);
fprintf('The selected grid size is: %.4f meters\n', gridSize);

Rsnew_pC = pointCloud(Rsnew);
Rsnew_pCdown = pcdownsample(Rsnew_pC, 'gridAverage', gridSize);

%% Downsampled 3D point cloud
figure(1);
pcshow(Rsnew_pCdown);
set(1, 'name', 'Raw data PC downsampled', 'Color', [1 1 1]);
set(gca, 'fontweight', 'bold', 'fontsize', 12, 'fontname', 'Times New Roman');
title('Raw Data Point Cloud', 'FontWeight', 'bold', 'FontSize', 16, 'FontName', 'Times New Roman');
xlabel('x (m)', 'FontWeight', 'bold', 'FontSize', 12, 'FontName', 'Times New Roman');
ylabel('y (m)', 'FontWeight', 'bold', 'FontSize', 12, 'FontName', 'Times New Roman');
zlabel('z (m)', 'FontWeight', 'bold', 'FontSize', 12, 'FontName', 'Times New Roman');
colorbar;

%% Downsampled 2D grid, with orientation markers
figure(2);
scatter(Rsnew_pCdown.Location(:,1), Rsnew_pCdown.Location(:,2), 2, Rsnew_pCdown.Location(:,3));
set(2, 'name', 'Raw data 2D grid downsampled', 'Color', [1 1 1]);
set(gca, 'fontweight', 'bold', 'fontsize', 12, 'fontname', 'Times New Roman');
title('Raw Data 2D Grid', 'FontWeight', 'bold', 'FontSize', 16, 'FontName', 'Times New Roman');
xlabel('x (m)', 'FontWeight', 'bold', 'FontSize', 12, 'FontName', 'Times New Roman');
ylabel('y (m)', 'FontWeight', 'bold', 'FontSize', 12, 'FontName', 'Times New Roman');
colorbar; grid on; axis equal;
xlim([min(Rsnew_pCdown.Location(:,1)) max(Rsnew_pCdown.Location(:,1))]);
ylim([min(Rsnew_pCdown.Location(:,2)) max(Rsnew_pCdown.Location(:,2))]);

hold on
r = max(Rsnew(:,1)) / 40;
fi = linspace(0, 2*pi, 360);
fx = r * cos(fi); fy = r * sin(fi);
plot(0, 0, '+');
plot(fx, fy);
fi1 = linspace(0, 0.5*pi, 10);
fx1 = r * cos(fi1); fy1 = r * sin(fi1);
quiver(fx1, fy1, fx1, fy1, 80, 'r');
fi2 = linspace(0, 2*pi, 5);
fx2 = r * cos(fi2); fy2 = r * sin(fi2);
quiver(fx2, fy2, fx2, fy2, 60, 'b');
hold off
end
