function [R, d, isMain] = generateBulgeDOF(R, d, isMain, isCircular)

n_nt = numel(R);
assert(n_nt == numel(d));
% assert(n_nt == numel(isMain));

% Find all single-stranded region
ss = find_ss(R, isCircular);
% Remove single-stranded regions with free ends
if(~isCircular)
    for i = numel(ss) : -1 : 1
        if(ss{i}(1) == 1 || ss{i}(end) == n_nt)
            ss(i) = [];
        end
    end
end

% Find positions and orientations of the single-stranded regions
for i = 1 : numel(ss)
    i_1 = ss{i}(1) - 1;
    if(i_1 < 1)
        assert(isCircular);
        i_1 = n_nt;
    end
    R_1 = R{i_1};
    if(~isMain(i_1))
        R_1 = R_1 * vrrotvec2mat([0 0 1 pi]);
    end
    d_1 = d{i_1};
    
    i_2 = ss{i}(end) + 1;
    if(i_2 > n_nt)
        assert(isCircular);
        i_2 = 1;
    end
    R_2 = R{i_2};
    if(~isMain(i_2))
        R_2 = R_2 * vrrotvec2mat([0 0 1 pi]);
    end
    d_2 = d{i_2};
    
    % Find the positions and orientations
    [R_fit, d_fit] = fit_R_d(R_1, d_1, R_2, d_2, numel(ss{i}));
    
    % Update the data
    for j = 1 : numel(ss{i})
        isMain(ss{i}(j)) = true;
        R{ss{i}(j)} = R_fit(:,:,j);
        d{ss{i}(j)} = d_fit(:,j);
    end
end

end


function ss = find_ss(R, isCircular)

loop = 0;
ss = [];

n_nt = numel(R);
isVisited = false(size(R));
for i = 1 : n_nt
    if(~isempty(R{i}))
        isVisited(i) = true;
    end
end
assert(~isempty(find(isVisited,1)));   % Exist at least one dsDNA region

% Scan for single-stranded regions
while(find(~isVisited,1))
    loop = loop + 1;
    ss{loop} = [];
    i_0 = find(~isVisited,1);
    
    % Scen in the 5'-direction
    i = i_0;
    while(~isVisited(i))
        isVisited(i) = true;
        ss{loop} = cat(2, i, ss{loop});
        if(i > 1)
            i = i - 1;
        elseif(isCircular)
            i = n_nt;
        else
            break;
        end 
    end
    ss{loop}(end) = [];
    isVisited(i_0) = false;

    % Scan in the 3'-direction
    i = i_0;
    while(~isVisited(i))
        isVisited(i) = true;
        ss{loop} = cat(2, ss{loop}, i);
        if(i < n_nt)
            i = i + 1;
        elseif(isCircular)
            i = 1;
        else
            break;
        end
    end
end 

end


function [R_fit, d_fit] = fit_R_d(R_1, d_1, R_2, d_2, n_ss)

R_fit = zeros(3,3,n_ss);
d_fit = zeros(3,n_ss);

% Rotation matrix R
% R * R_1 = R_2
R = R_2 / R_1;
v = vrrotmat2vec(R);
a = v(1:3);
theta = v(4);

% Calculate for R_fit and d_fit
for i = 1 : n_ss
    R_fit(:,:,i) = vrrotvec2mat([a, theta*i/(n_ss+1)]) * R_1;
    d_fit(:,i) = (d_1*(n_ss+1-i) + d_2*i) / (n_ss+1);
end

end