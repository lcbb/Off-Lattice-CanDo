function dnaTop = Tiamat2topology(filename)

% Open the file
fid_top = fopen(filename);

% Read the number of bases
tline = fgetl(fid_top);
nBase = sscanf(tline, '%d');

% Raw data from Tiamat
raw.id = -ones(nBase,1);
raw.up = -ones(nBase,1);         % 5' neighbor
raw.down = -ones(nBase,1);       % 3' neighbor
raw.across = -ones(nBase,1);     % Watson-Crick neighbor

% DNA topology
dnaTop(nBase).id = -1;
dnaTop(nBase).up = -1;           % 5' neighbor
dnaTop(nBase).down = -1;         % 3' neighbor
dnaTop(nBase).across = -1;       % Watson-Crick neighbor
dnaTop(nBase).xyz = zeros(3,1);  % coordinate in a Euclidean space

% Read the file
% See http://stackoverflow.com/questions/9440592/fastest-matlab-file-reading
for i = 1:nBase
    tline = fgetl(fid_top);
    buffer = sscanf(tline, '%d%x%x%x%x%f%f%f');
    %fprintf('%s\n', tline);
    
    raw.id(i) = buffer(2);
    raw.up(i) = buffer(3);
    raw.down(i) = buffer(4);
    raw.across(i) = buffer(5);
    
    dnaTop(i).id = buffer(1);
    dnaTop(i).xyz = buffer(6:8);
end

% Fix the ID for each base
for i = 1:nBase
    tmp = find(raw.up == raw.id(i));
    if(numel(tmp)==1)
        dnaTop(i).up = tmp;
    elseif(numel(tmp)==0)
        dnaTop(i).up = -1;
    else
        error('Duplication.');
    end
    
    tmp = find(raw.down == raw.id(i));
    if(numel(tmp)==1)
        dnaTop(i).down = tmp;
    elseif(numel(tmp)==0)
        dnaTop(i).down = -1;
    else
        error('Duplication.');
    end

    tmp = find(raw.across == raw.id(i));
    if(numel(tmp)==1)
        dnaTop(i).across = tmp;
    elseif(numel(tmp)==0)
        dnaTop(i).across = -1;
    else
        error('Duplication.');
    end
end

% Clean up
fclose(fid_top);

end