function dnaInfo = cndo2dnaInfo(cndo_FN)

tol = 1e-10;

%% Initialization
conn = zeros(0,5);      % id, up, down, across
seq = cell(0,1);        % seq
dNode = zeros(0,4);
triad = zeros(0,10);
id_nt = zeros(0,3);


%% Read the file
fid = fopen(cndo_FN);

% Ignore 2 lines
for i = 1 : 2
    ss = strtrim(fgetl(fid));
    assert(ischar(ss));
end

% Read the field 'dnaTop'
ss = strtrim(fgetl(fid));
assert(strcmp(ss(1:6), 'dnaTop'));
ss = strtrim(fgetl(fid));
while(~isempty(ss));
    a = strsplit(ss, ',');
    assert(numel(a) == 6);
    conn = cat(1, conn, [str2double(a{1}), ...
                         str2double(a{2}), ...
                         str2double(a{3}), ...
                         str2double(a{4}), ...
                         str2double(a{5})]);
    seq = cat(1, seq, a(6));
    ss = strtrim(fgetl(fid));
end

% Read the field 'dNode'
ss = strtrim(fgetl(fid));
assert(strcmp(ss(1:5), 'dNode'));
ss = strtrim(fgetl(fid));
while(~isempty(ss))
    a = strsplit(ss, ',');
    assert(numel(a) == 4);
    dNode = cat(1, dNode, [str2double(a{1}), ...
                           str2double(a{2}), ...
                           str2double(a{3}), ...
                           str2double(a{4})]);
    ss = strtrim(fgetl(fid));
end

% Read the field 'triad'
ss = strtrim(fgetl(fid));
assert(strcmp(ss(1:5), 'triad'));
ss = strtrim(fgetl(fid));
while(~isempty(ss))
    a = strsplit(ss, ',');
    assert(numel(a) == 10);
    triad = cat(1, triad, [str2double(a{1}), ...
                           str2double(a{2}), ...
                           str2double(a{3}), ...
                           str2double(a{4}), ...
                           str2double(a{5}), ...
                           str2double(a{6}), ...
                           str2double(a{7}), ...
                           str2double(a{8}), ...
                           str2double(a{9}), ...
                           str2double(a{10})]);
    ss = strtrim(fgetl(fid));
end

% Read the field 'id_nt'
ss = strtrim(fgetl(fid));
assert(strcmp(ss(1:5), 'id_nt'));
ss = strtrim(fgetl(fid));
while(~isempty(ss))
    a = strsplit(ss, ',');
    assert(numel(a) == 3);
    id_nt = cat(1, id_nt, [str2double(a{1}), ...
                           str2double(a{2}), ...
                           str2double(a{3})]);
    ss = fgetl(fid);
    if(~ischar(ss))
        break;
    end
    ss = strtrim(ss);
end

fclose(fid);


%% Generate the MATLAB script 'dnaInfo'
n_nt = size(conn,1);
assert(n_nt == numel(seq));
assert(norm(conn(:,1) - (1:n_nt)') < tol);
conn(:,1) = [];

n_bp = size(dNode,1);
assert(n_bp == size(triad,1) && n_bp == size(id_nt,1));
dNode(:,1) = [];
triad(:,1) = [];
id_nt(:,1) = [];

triad2 = zeros(3,3,n_bp);
for i = 1 : n_bp
    triad2(:,:,i) = reshape(triad(i,:), 3, 3);
end

for i = 1 : n_nt
    dnaInfo.dnaTop(i) = struct('id', conn(i,1), ...
                               'up', conn(i,2), ...
                               'down', conn(i,3), ...
                               'across', conn(i,4), ...
                               'seq', seq{i});
    dnaInfo.dnaGeom.dNode = dNode;
    dnaInfo.dnaGeom.triad = triad2;
    dnaInfo.dnaGeom.id_nt = id_nt;
end

end

