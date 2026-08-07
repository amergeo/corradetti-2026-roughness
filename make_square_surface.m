%% make_square_surface.m
% Crops the interpolated surface grid (qz, sampFFT, Xq, Yq, filename -
% already in workspace after Interpolate_PC) to a square by trimming the
% longer dimension symmetrically about its centre, and saves the result
% as a tab-separated xyz .txt file compatible with import_PC.m.
if ~exist('qz','var') || ~exist('sampFFT','var') || ~exist('Xq','var') || ~exist('Yq','var')
    error(['Required variables not found in workspace. ' ...
           'Run the alignment/interpolation steps first.']);
end
[ny, nx] = size(qz);
fprintf('Current grid: %d (Y) x %d (X) points  |  dx = %.4g m\n', ny, nx, sampFFT);

%% Centre-symmetric crop
nSq = min(ny, nx);
y_start = floor((ny - nSq) / 2) + 1;
y_end   = y_start + nSq - 1;
x_start = floor((nx - nSq) / 2) + 1;
x_end   = x_start + nSq - 1;
qz_sq = qz(y_start:y_end, x_start:x_end);

x_vec = (0:nSq-1) * sampFFT;
y_vec = (0:nSq-1) * sampFFT;
[Xq_sq, Yq_sq] = meshgrid(x_vec, y_vec);
Lx_sq = (nSq - 1) * sampFFT;
fprintf('Square grid:  %d x %d points  |  side = %.4g m\n', nSq, nSq, Lx_sq);

%% Save as tab-separated xyz txt
outName = sprintf('%s_square.txt', filename);
[saveName, savePath] = uiputfile('*.txt', 'Save square surface as', outName);
if isequal(saveName, 0)
    disp('Save cancelled by user.'); return
end
outPath = fullfile(savePath, saveName);
fid = fopen(outPath, 'w');
if fid == -1
    error('Cannot open output file for writing: %s', outPath);
end
fprintf(fid, '%.8f\t%.8f\t%.8f\n', [Xq_sq(:), Yq_sq(:), qz_sq(:)]');
fclose(fid);
fprintf('>>> Saved: %s\n', outPath);
fprintf('    Grid: %d x %d  |  Lx = Ly = %.4g m  |  dx = %.4g m\n', nSq, nSq, Lx_sq, sampFFT);

%% Push to workspace
assignin('base', 'qz_square',  qz_sq);
assignin('base', 'Xq_square',  Xq_sq);
assignin('base', 'Yq_square',  Yq_sq);
assignin('base', 'filename_square', saveName(1:end-4));
