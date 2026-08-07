function [min_freq_user, max_freq_user, num_bins_user] = freqRangeGUI(min_freq, max_freq, askBins)
% Lets the user set a min/max spatial-frequency window (used by the
% userbounded_progmean fit variants). askBins (default true) toggles a
% "number of bins" field, not needed for the boundary-only use case;
% when omitted num_bins_user is returned as NaN.
if nargin < 3
    askBins = true;
end

fig = uifigure('Name', 'Frequency Range Selection', 'Position', [500, 300, 400, 350]);
fig.UserData.min_freq_user = min_freq;
fig.UserData.max_freq_user = max_freq;
fig.UserData.num_bins_user = 10;

uilabel(fig, 'Text', 'Frequency Range Selection', 'FontSize', 14, ...
    'FontWeight', 'bold', 'Position', [100, 310, 400, 30]);
uilabel(fig, 'Text', sprintf('Min Frequency: %.4f', min_freq), 'Position', [50, 260, 300, 30]);
uilabel(fig, 'Text', sprintf('Max Frequency: %.4f', max_freq), 'Position', [50, 220, 300, 30]);
uilabel(fig, 'Text', 'Enter new Min Frequency:', 'Position', [50, 180, 150, 30]);
min_freq_field = uieditfield(fig, 'numeric', 'Position', [200, 180, 150, 30], 'Value', min_freq);
uilabel(fig, 'Text', 'Enter new Max Frequency:', 'Position', [50, 140, 150, 30]);
max_freq_field = uieditfield(fig, 'numeric', 'Position', [200, 140, 150, 30], 'Value', max_freq);

if askBins
    uilabel(fig, 'Text', 'Enter number of bins:', 'Position', [50, 100, 150, 30]);
    num_bins_field = uieditfield(fig, 'numeric', 'Position', [200, 100, 150, 30], 'Value', 10);
else
    num_bins_field = [];
end

uibutton(fig, 'Text', 'OK', 'Position', [150, 50, 100, 30], ...
    'ButtonPushedFcn', @(btn, event) confirmSelection(min_freq_field, max_freq_field, num_bins_field, fig));

uiwait(fig);
min_freq_user = fig.UserData.min_freq_user;
max_freq_user = fig.UserData.max_freq_user;
num_bins_user = fig.UserData.num_bins_user;
close(fig);
end

function confirmSelection(min_freq_field, max_freq_field, num_bins_field, fig)
    fig.UserData.min_freq_user = min_freq_field.Value;
    fig.UserData.max_freq_user = max_freq_field.Value;
    if ~isempty(num_bins_field)
        fig.UserData.num_bins_user = num_bins_field.Value;
    else
        fig.UserData.num_bins_user = NaN;
    end
    uiresume(fig);
end
