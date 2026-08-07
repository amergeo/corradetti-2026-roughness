%% compare_profiles_multizoom_stessi_assi.m
% Compares central profiles from up to 4 surfaces at two zoom levels
% (100% and 10% of profile length), with a shared Z scale within each
% zoom row so amplitude is directly comparable across surfaces. For each
% profile, computes H, alpha and Rq_theoretical via single-profile PSD
% with Corradetti (2017) geometric-progression binning (Fig. 5).
% If both directions are requested: Figure 1 = Y profiles, Figure 2 = X profiles.

%% 0. User settings
nSurfaces = input('How many surfaces to compare? (1-4): ');
if nSurfaces < 1 || nSurfaces > 4
    error('Please enter a number between 1 and 4.');
end
dirChoice = input('Direction(s) to analyse - type X, Y, or B (both): ', 's');
dirChoice = upper(strtrim(dirChoice));
if ~ismember(dirChoice, {'X','Y','B'})
    error('Please enter X, Y, or B.');
end
if strcmp(dirChoice, 'B')
    directions = {'Y','X'};
elseif strcmp(dirChoice, 'Y')
    directions = {'Y'};
else
    directions = {'X'};
end

zoomWindows = [0.00, 1.00;
               0.00, 0.10];
nZooms = size(zoomWindows, 1);
taperFrac = 0.03;

%% 1. Load surfaces (once), store full grid for both directions
grids  = cell(nSurfaces, 1);
xus    = cell(nSurfaces, 1);
yus    = cell(nSurfaces, 1);
labels = cell(nSurfaces, 1);
for s = 1:nSurfaces
    fprintf('\n--- Surface %d ---\n', s);
    [fname, fpath] = uigetfile({'*.xyz;*.txt', 'Surface files (*.xyz, *.txt)'}, ...
        sprintf('Select surface %d of %d', s, nSurfaces));
    if isequal(fname, 0); error('User cancelled at surface %d.', s); end
    raw = load(fullfile(fpath, fname));
    x_col = raw(:,1); y_col = raw(:,2); z_col = raw(:,3);
    xu = unique(x_col, 'sorted'); yu = unique(y_col, 'sorted');
    nx = numel(xu); ny = numel(yu);
    if nx * ny ~= numel(z_col)
        error('Cannot reconstruct regular grid for %s.', fname);
    end
    grids{s} = reshape(z_col, ny, nx);
    xus{s} = xu; yus{s} = yu;
    labels{s} = fname(1:end-4);
    fprintf('    Loaded: %d x %d  |  Rq_true = %.4e m\n', ny, nx, sqrt(mean(z_col.^2)));
end

%% 2. Global Z limits per zoom row (shared across all surfaces and directions)
zlims = zeros(nZooms, 2);
for z = 1:nZooms
    zmin_all = Inf; zmax_all = -Inf;
    for s = 1:nSurfaces
        Z = grids{s}; xu = xus{s}; yu = yus{s};
        [~, iy] = min(abs(yu - (yu(end)-yu(1))/2));
        profX = Z(iy,:) - mean(Z(iy,:));
        Lx = xu(end)-xu(1);
        x0 = xu(1) + zoomWindows(z,1)*Lx; x1 = xu(1) + zoomWindows(z,2)*Lx;
        maskX = xu >= x0 & xu <= x1;
        zmin_all = min(zmin_all, min(profX(maskX))); zmax_all = max(zmax_all, max(profX(maskX)));

        [~, ix] = min(abs(xu - (xu(end)-xu(1))/2));
        profY = Z(:,ix)' - mean(Z(:,ix));
        Ly = yu(end)-yu(1);
        y0 = yu(1) + zoomWindows(z,1)*Ly; y1 = yu(1) + zoomWindows(z,2)*Ly;
        maskY = yu >= y0 & yu <= y1;
        zmin_all = min(zmin_all, min(profY(maskY))); zmax_all = max(zmax_all, max(profY(maskY)));
    end
    zlims(z,:) = [zmin_all, zmax_all];
end

%% 3. Loop over requested directions
for dIdx = 1:numel(directions)
    thisDir = directions{dIdx};
    figNum = dIdx;
    profiles = cell(nSurfaces, 1);
    coordVecs = cell(nSurfaces, 1);
    Hfit = zeros(nSurfaces, 1);
    afit = zeros(nSurfaces, 1);
    Rq_theo = zeros(nSurfaces, 1);
    Rq_true = zeros(nSurfaces, 1);

    for s = 1:nSurfaces
        Z = grids{s}; xu = xus{s}; yu = yus{s};
        [ny, nx] = size(Z);
        Rq_true(s) = sqrt(mean(Z(:).^2));

        if strcmp(thisDir, 'X')
            [~, idx] = min(abs(yu - (yu(end)-yu(1))/2));
            prof = Z(idx, :); coord = xu; dx = mean(diff(xu));
        else
            [~, idx] = min(abs(xu - (xu(end)-xu(1))/2));
            prof = Z(:, idx)'; coord = yu; dx = mean(diff(yu));
        end
        prof = prof - mean(prof);
        profiles{s} = prof; coordVecs{s} = coord;

        %% Single-profile PSD with Corradetti (2017) binning
        N = numel(prof);
        nTap = round(taperFrac * N);
        win = ones(1, N);
        win(1:nTap) = linspace(0, 1, nTap);
        win(end-nTap+1:end) = linspace(1, 0, nTap);
        cf = 1 / mean(win);
        prof_w = prof .* win * cf;
        F = fft(prof_w) / N;
        if rem(N,2)==0; cut=N/2+1; else; cut=(N+1)/2; end
        Fc = F(1:cut);
        A = abs(Fc);
        if rem(N,2)==0; A(2:end-1)=A(2:end-1)*2; else; A(2:end)=A(2:end)*2; end
        A = A * cf;
        freq = (0:cut-1) / (N*dx);
        psd = (A.^2) * (N*dx);
        freq = freq(2:end)'; psd = psd(2:end)';

        nPoints = numel(freq);
        maxFreq = max(freq);
        progmean = zeros(3, 2);
        progmean(1:3, 1) = freq(1:3);
        progmean(1:3, 2) = psd(1:3);
        strep = log10(freq(4)) - log10(freq(3));
        eligible = false(nPoints, 1);
        eligible(4:end) = true;
        freqLog = log10(freq);
        binIdx = 3;
        t = log10(freq(3));
        while t < log10(maxFreq)
            binIdx = binIdx + 1;
            inBin = eligible & (freqLog >= t) & (freqLog < t + strep);
            rows = [freq(inBin), psd(inBin)];
            if size(rows,1) == 1
                progmean(binIdx,:) = rows(1,:);
            elseif size(rows,1) > 1
                progmean(binIdx,1) = sqrt(10^t * 10^(t+strep));
                progmean(binIdx,2) = mean(rows(:,2));
            else
                progmean(binIdx,:) = [NaN, NaN];
            end
            t = t + strep;
        end
        pmean = progmean(~isnan(progmean(:,2)), :);
        fitobj = fit(log10(pmean(:,1)), log10(pmean(:,2)), 'poly1');
        coef = coeffvalues(fitobj);
        Hfit(s) = (coef(1)+1)/(-2);
        afit(s) = 10^coef(2);
        fmin = min(pmean(:,1)); fmax = max(pmean(:,1));
        Rq_theo(s) = sqrt( afit(s)/(2*Hfit(s)) * (fmin^(-2*Hfit(s)) - fmax^(-2*Hfit(s))) );
        fprintf('    [%s] H=%.4f  a=%.4e  Rq_th=%.4e  Rq_true=%.4e\n', ...
            thisDir, Hfit(s), afit(s), Rq_theo(s), Rq_true(s));
    end

    %% Figure
    fig = figure(figNum); clf(fig);
    set(fig, 'Color', [1 1 1], 'Name', sprintf('Profile comparison - along %s - shared Z per zoom', thisDir));
    fig.Position = [50 + (figNum-1)*30, 50, 300*nSurfaces, 250*nZooms];
    axLabel = thisDir;
    for z = 1:nZooms
        for s = 1:nSurfaces
            ax = subplot(nZooms, nSurfaces, (z-1)*nSurfaces + s);
            cv = coordVecs{s}; pro = profiles{s};
            Lp = cv(end) - cv(1);
            x0 = cv(1) + zoomWindows(z,1) * Lp;
            x1 = cv(1) + zoomWindows(z,2) * Lp;
            mask = cv >= x0 & cv <= x1;
            plot(cv(mask), pro(mask), 'k-', 'LineWidth', 0.8);
            ylim(zlims(z,:));
            set(ax, 'Color', [1 1 1], 'XColor', [0 0 0], 'YColor', [0 0 0]);
            set(ax, 'FontName', 'Times New Roman', 'FontSize', 10, 'FontWeight', 'bold');
            xlabel(sprintf('%s (m)', axLabel), 'FontName', 'Times New Roman', 'FontSize', 10, 'FontWeight', 'bold');
            if z == 1
                title(sprintf('%s\nH=%.2f | \\alpha=%.2e m^3 | Rq_{th}=%.2e m', ...
                    labels{s}, Hfit(s), afit(s), Rq_theo(s)), ...
                    'FontName', 'Times New Roman', 'FontSize', 9, 'FontWeight', 'bold', ...
                    'Color', [0 0 0], 'Interpreter', 'tex');
            end
            if s == 1
                zoomPct = (zoomWindows(z,2) - zoomWindows(z,1)) * 100;
                ylabel(sprintf('Z (m)  [%.0f%% of L_%s]', zoomPct, axLabel), ...
                    'FontName', 'Times New Roman', 'FontSize', 9, 'FontWeight', 'bold');
            else
                ylabel('Z (m)', 'FontName', 'Times New Roman', 'FontSize', 9, 'FontWeight', 'bold');
            end
        end
    end
    sgtitle(sprintf('Profile comparison along %s - shared Z scale within each zoom level', thisDir), ...
        'FontName', 'Times New Roman', 'FontSize', 13, 'FontWeight', 'bold');
end
