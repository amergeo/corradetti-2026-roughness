% alignZ.m - Rotates a point cloud so that Z is normal to its
% least-squares best-fit plane (parallel to the smallest eigenvector).
function [Rs, Rsnew, x0, a, d, normd, R, xyn, zn, r_deg] = alignZ(Raw_surface)

Rs = table2array(Raw_surface);
[x0, a, d, normd] = lsplane(Rs);

Vz = [0 0 1];
if a(3) < 0
    a = -a;   % ensure normal points upward
end
Vnz = [a(1), a(2), a(3)];
Vnx = cross(Vz, Vnz) / norm(cross(Vz, Vnz));
Vny = cross(Vnz, Vnx) / norm(cross(Vnz, Vnx));
R = [Vnx; Vny; Vnz];

m = size(Rs, 1);
for i = 1:m
    Rsnew(:, i) = R * (Rs(i, :) - x0')';
end
Rsnew = Rsnew';

xyn = [Rsnew(:, 1), Rsnew(:, 2)];
zn = Rsnew(:, 3);
r_deg = 0;  % initial value; updated later by alignXY

end

% Least-squares plane via SVD (orthogonal distance regression).
% Based on LSPLANE.M by I M Smith, 2002.
function [x0, a, d, normd] = lsplane(Rs)
    m = size(Rs, 1);
    if m < 3
        error('At least 3 data points are required.');
    end
    x0 = mean(Rs)';
    A = [(Rs(:, 1) - x0(1)), (Rs(:, 2) - x0(2)), (Rs(:, 3) - x0(3))];
    [~, S, V] = svd(A, 0);
    [~, i] = min(diag(S));
    a = V(:, i);
    if nargout > 2
        d = A * a;
        normd = norm(d);
    end
end
