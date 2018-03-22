function dnaTop = dnaInfo2dnaTop(dnaInfo)

% Parameters
param.diameterHX = 2.0;                 % [nm] diameter of a helix
param.angMinor = 120;                   % [degree] angle of the minor groove

dnaTop = dnaInfo.dnaTop;
dNode = dnaInfo.dnaGeom.dNode / 10;     % unit: nm
triad = dnaInfo.dnaGeom.triad;
id_nt = dnaInfo.dnaGeom.id_nt;

% Sizes
assert(size(dnaTop, 1) == 1);
n_nt = size(dnaTop, 2);
n_bp = size(dNode, 1);
assert(size(dNode, 2) == 3);
assert(size(triad, 1) == 3 && size(triad, 2) == 3 && size(triad, 3) == n_bp);
assert(size(id_nt, 1) == n_bp && size(id_nt, 2) == 2);

% Initialize the field dnaTop(i).xyz
for i = 1 : n_nt
    assert(~isempty(dnaTop(i).id) && ~isempty(dnaTop(i).up) && ...
           ~isempty(dnaTop(i).down) && ~isempty(dnaTop(i).across) && ...
           ~isempty(dnaTop(i).seq));
%     dnaTop(i).xyz = NaN(3, 1);
%     dnaTop(i).xyz = zeros(3, 1);
    dnaTop(i).xyz = randn(3, 1);
end

for i = 1 : n_bp
    i_prefer = id_nt(i, 1);
    i_other = id_nt(i, 2);
    
    [xyz_prefer, xyz_other] = nt_pos(dNode(i, :), triad(:, :, i), param);
%     assert(sum(isnan([dnaTop(i_prefer).xyz; dnaTop(i_other).xyz])) == 6);
    dnaTop(i_prefer).xyz = xyz_prefer;
    dnaTop(i_other).xyz = xyz_other;
end

end


function [xyz_prefer, xyz_other] = nt_pos(dNode, triad, param)

assert(size(dNode, 1) == 1 && size(dNode, 2) == 3);
assert(size(triad, 1) == 3 && size(triad, 2) == 3);
r = param.diameterHX / 2;
a = param.angMinor;
a1 = deg2rad(180 - a / 2);
a2 = deg2rad(180 + a / 2);

% In the local frame
xyz_1 = [r * cos(a1), r * sin(a1), 0]';
xyz_2 = [r * cos(a2), r * sin(a2), 0]';

% In the global frame
xyz_prefer = triad * xyz_1 + dNode';
xyz_other = triad * xyz_2 + dNode';

end