function [Raw_surface, filename] = generateSelfAffineSurface()
% Generates an anisotropic self-affine surface with target Hx, Hy and
% amplitude (Rq or P2V). Saves Run_XXXX.xyz and appends a row to
% selfaffine_runs_log.csv (created if absent).
close all;

%% 1. Determine next run number
LOG_FILE = 'selfaffine_runs_log.csv';
if ~isfile(LOG_FILE)
    fid = fopen(LOG_FILE, 'w');
    if fid == -1
        error('Cannot create log file "%s" in:\n  %s', LOG_FILE, pwd);
    end
    fprintf(fid, 'RunID\tLx_m\tLy_m\tdx_m\tdy_m\tHx\tHy\tseed\tampMode\ttargetRq_m\ttargetP2V_m\tachievedRq_m\tachievedP2V_m\tscaleFactor\tnx\tny\ttimestamp\tdx_eff_m\tdy_eff_m\n');
    fclose(fid);
    nextRun = 1;
else
    try
        T = readtable(LOG_FILE, 'Delimiter', '\t', 'TextType', 'string');
        if isempty(T)
            nextRun = 1;
        else
            ids = T.RunID;
            nums = zeros(numel(ids), 1);
            for k = 1:numel(ids)
                tok = regexp(ids(k), 'Run_(\d+)', 'tokens');
                if ~isempty(tok); nums(k) = str2double(tok{1}{1}); end
            end
            nextRun = max(nums) + 1;
        end
    catch
        warning('Could not parse log file. Scanning folder for Run_XXXX.xyz files instead.');
        existing = dir('Run_*.xyz');
        if isempty(existing)
            nextRun = 1;
        else
            nums = zeros(numel(existing), 1);
            for k = 1:numel(existing)
                tok = regexp(existing(k).name, 'Run_(\d+)', 'tokens');
                if ~isempty(tok); nums(k) = str2double(tok{1}{1}); end
            end
            nextRun = max(nums) + 1;
        end
    end
end
runID = sprintf('Run_%04d', nextRun);
fprintf('>>> Starting %s\n', runID);

%% 2. GUI - parameter input
dflt.Lx = 0.5; dflt.Ly = 0.5;
dflt.dx = 5e-4; dflt.dy = 5e-4;
dflt.Hx = 0.80; dflt.Hy = 0.60;
dflt.seed = 1;
dflt.ampMode = 'Rq';
dflt.targetRq = 1e-3; dflt.targetP2V = 2e-3;

prompt = {'Lx (m) [perpendicular to slip]', 'Ly (m) [along slip]', ...
    'dx (m) [grid spacing along x]', 'dy (m) [grid spacing along y]', ...
    'Hx (0-1)', 'Hy (0-1)', 'Random seed (integer)', ...
    'Amplitude mode: Rq or P2V', 'Target Rq (m)', 'Target P2V (m)'};
def_ans = {num2str(dflt.Lx), num2str(dflt.Ly), num2str(dflt.dx), num2str(dflt.dy), ...
    num2str(dflt.Hx), num2str(dflt.Hy), num2str(dflt.seed), dflt.ampMode, ...
    num2str(dflt.targetRq), num2str(dflt.targetP2V)};
answ = inputdlg(prompt, sprintf('Self-affine surface generator - %s', runID), [1 75], def_ans);
if isempty(answ)
    disp('Cancelled by user. No run registered.');
    Raw_surface = table(); filename = ''; return;
end

%% 3. Parse and validate
params.Lx = str2double(answ{1}); params.Ly = str2double(answ{2});
params.dx = str2double(answ{3}); params.dy = str2double(answ{4});
params.Hx = str2double(answ{5}); params.Hy = str2double(answ{6});
params.seed = round(str2double(answ{7}));
params.ampMode = strtrim(answ{8});
params.targetRq = str2double(answ{9}); params.targetP2V = str2double(answ{10});

mustBeFiniteScalar(params.Lx,'Lx'); mustBePositive(params.Lx,'Lx');
mustBeFiniteScalar(params.Ly,'Ly'); mustBePositive(params.Ly,'Ly');
mustBeFiniteScalar(params.dx,'dx'); mustBePositive(params.dx,'dx');
mustBeFiniteScalar(params.dy,'dy'); mustBePositive(params.dy,'dy');
mustBeFiniteScalar(params.Hx,'Hx'); mustBeInRange01(params.Hx,'Hx');
mustBeFiniteScalar(params.Hy,'Hy'); mustBeInRange01(params.Hy,'Hy');
mustBeFiniteScalar(params.seed,'seed');
if ~any(strcmpi(params.ampMode, {'Rq','P2V'}))
    error('Amplitude mode must be "Rq" or "P2V".');
end
if strcmpi(params.ampMode,'Rq')
    mustBeFiniteScalar(params.targetRq,'Target Rq'); mustBePositive(params.targetRq,'Target Rq');
else
    mustBeFiniteScalar(params.targetP2V,'Target P2V'); mustBePositive(params.targetP2V,'Target P2V');
end

%% 4. Grid
nx = max(2, floor(params.Lx / params.dx) + 1);
ny = max(2, floor(params.Ly / params.dy) + 1);
x = linspace(0, params.Lx, nx);
y = linspace(0, params.Ly, ny);
[X, Y] = meshgrid(x, y);
% Actual spacing used (may differ slightly from requested dx/dy)
params.dx_eff = params.Lx / (nx - 1);
params.dy_eff = params.Ly / (ny - 1);

%% 5. Anisotropic self-affine surface - Fourier domain
rng(params.seed, 'twister');
fx = ifftshift((-floor(nx/2):ceil(nx/2)-1) / params.Lx);
fy = ifftshift((-floor(ny/2):ceil(ny/2)-1) / params.Ly);
[FX, FY] = meshgrid(fx, fy);
kx = max(abs(FX), 1/params.Lx);
ky = max(abs(FY), 1/params.Ly);
A = (kx).^(-(params.Hx + 0.5)) .* (ky).^(-(params.Hy + 0.5));  % separable anisotropic filter
A(FX==0 & FY==0) = 0;   % remove DC
G = randn(ny, nx) + 1i*randn(ny, nx);
Z0 = ifft2(G .* A, 'symmetric');
if any(~isfinite(Z0(:)))
    error('Generated surface contains NaN/Inf. Try different H values or grid.');
end
Z0 = double(Z0 - mean(Z0(:)));

%% 6. Impose physical amplitude
switch upper(params.ampMode)
    case 'RQ'
        currentRq = sqrt(mean(Z0(:).^2));
        if currentRq == 0; error('RMS is zero. Try a different seed.'); end
        s = params.targetRq / currentRq;
    case 'P2V'
        currentP2V = max(Z0(:)) - min(Z0(:));
        if currentP2V == 0; error('P2V is zero. Try a different seed.'); end
        s = params.targetP2V / currentP2V;
end
Z = Z0 * s;
params.achievedRq = sqrt(mean(Z(:).^2));
params.achievedP2V = max(Z(:)) - min(Z(:));
params.scaleFactorApplied = s;

%% 7. Build Raw_surface table (import_PC compatible)
Raw_surface = table(X(:), Y(:), Z(:), 'VariableNames', {'x_col','y_col','z_col'});
filename = runID;

%% 8. Save surface as XYZ
xyzName = [runID, '.xyz'];
try
    fid = fopen(xyzName, 'w');
    if fid == -1; error('Cannot open "%s" for writing in:\n  %s', xyzName, pwd); end
    fprintf(fid, '%.8f\t%.8f\t%.8f\n', [X(:), Y(:), Z(:)]');
    fclose(fid);
    fprintf('>>> Surface saved: %s\n', xyzName);
catch ME
    warning(ME.identifier, '%s', ME.message);
end

%% 9. Append row to log CSV
try
    fid = fopen(LOG_FILE, 'a');
    if fid == -1; error('Cannot open log file for appending.'); end
    timestamp = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
    fprintf(fid, '%s\t%.6g\t%.6g\t%.2e\t%.2e\t%.4f\t%.4f\t%d\t%s\t%.4e\t%.4e\t%.4e\t%.4e\t%.6g\t%d\t%d\t%s\t%.6e\t%.6e\n', ...
        runID, params.Lx, params.Ly, params.dx, params.dy, params.Hx, params.Hy, params.seed, ...
        upper(params.ampMode), params.targetRq, params.targetP2V, params.achievedRq, params.achievedP2V, ...
        params.scaleFactorApplied, nx, ny, timestamp, params.dx_eff, params.dy_eff);
    fclose(fid);
    fprintf('>>> Log updated: %s\n', LOG_FILE);
catch ME
    warning(ME.identifier, '%s', ME.message);
end

%% 10. Figures (requires makeTopoColormap on path; falls back to parula)
try
    topoCmap = makeTopoColormap();
catch
    topoCmap = parula(256);
end

figure(1); clf;
surf(X, Y, Z, 'EdgeColor','none');
axis tight; view(45,35); grid on; shading interp;
camlight headlight; lighting gouraud;
xlabel('X (m) - perpendicular to slip'); ylabel('Y (m) - along slip'); zlabel('Z (m) - elevation');
title(sprintf('%s | Hx=%.2f, Hy=%.2f | %s target', runID, params.Hx, params.Hy, upper(params.ampMode)), 'Interpreter','none');
colormap(topoCmap); colorbar;

figure(2); clf;
imagesc(x, y, Z); set(gca,'YDir','normal'); axis image;
xlabel('X (m) - perpendicular to slip'); ylabel('Y (m) - along slip');
title(sprintf('%s | achieved Rq=%.3g m, P2V=%.3g m', runID, params.achievedRq, params.achievedP2V), 'Interpreter','none');
colormap(topoCmap); colorbar;

%% 11. Push to base workspace
assignin('base','Raw_surface', Raw_surface);
assignin('base','filename', filename);
assignin('base','X', X); assignin('base','Y', Y); assignin('base','Z', Z);
assignin('base','params', params); assignin('base','runID', runID);
fprintf('>>> Done. runID=%s | Hx=%.3f Hy=%.3f Rq=%.4g m P2V=%.4g m\n', ...
    runID, params.Hx, params.Hy, params.achievedRq, params.achievedP2V);
end

%% Local validation helpers
function mustBeFiniteScalar(v, name)
    if ~(isscalar(v) && isfinite(v)); error('%s must be a finite scalar.', name); end
end
function mustBePositive(v, name)
    if v <= 0; error('%s must be > 0.', name); end
end
function mustBeInRange01(v, name)
    if v <= 0 || v >= 1; error('%s must be strictly in (0, 1).', name); end
end
