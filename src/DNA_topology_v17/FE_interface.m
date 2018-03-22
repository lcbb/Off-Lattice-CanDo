function [HJE, connectivity, strand4FE, nick4FE, bulge4FE, connSS4FE, connBB4FE] ...
         = FE_interface(dnaTop, strand, dsDNA, conn_es, conn_gap, ...
                        junction, junction_non4way, checktable)


%% Calculate the type, position, and orientation of each HJE
n_HJE = numel(junction);
HJE = struct('nWay', {}, 'type', {}, ...
             'c', {}, 'e1', {}, 'e2', {}, 'e3', {}, ...
             'armLength', {}, 'idSeq_main', {}, 'idSeq_comp', {});
i_dsxover = [];
i_ssxover = [];
for i = 1 : n_HJE
    if(junction(i).isCircular && junction(i).nWay == 4)
        i_dsxover = cat(2, i_dsxover, i);
        HJE(i).type = 'ds';
    elseif(~junction(i).isCircular && junction(i).nWay == 4)
        i_ssxover = cat(2, i_ssxover, i);
        HJE(i).type = 'ss';
    else
        error('Exception.');
    end
        
    HJE(i).nWay = 4;
    [c, e1, e2, e3] = junctionTriad(dnaTop, dsDNA, junction(i));
    HJE(i).c = c;
    HJE(i).e1 = e1;
    HJE(i).e2 = e2;
    HJE(i).e3 = e3;
end
assert(n_HJE == numel([i_dsxover, i_ssxover]));

arm_list = cell(n_HJE, 4);
duplex_list = cell(0, 1);
nick_list = cell(n_HJE, 4);
bulge_list = cell(n_HJE, 4);

nick4FE = zeros(0, 3);
bulge4FE = zeros(0, 3);
connSS4FE = zeros(0, 5);
connBB4FE = zeros(0, 4);


%% Build the connectivity between HJEs
% Use (1) nick_real, (2) bulge to "glue" dsDNA into ds_conn
connectivity = zeros(0, 7);
isVisited = false(size(checktable, 1), 1);

% For each junction arm, find its neighbor as another junction arm.
n_junction = numel(junction);
for i = 1 : n_junction
    assert(junction(i).nWay == 4);
    for j = 1 : junction(i).nWay
        % Start from junction i, arm j
        i_left = find(cellfun(@(x) find_cell(x, i), checktable(:, 1)) & ...
                      cellfun(@(y) find_cell(y, j), checktable(:, 2)));
        i_right = find(cellfun(@(x) find_cell(x, i), checktable(:, 4)) & ...
                       cellfun(@(y) find_cell(y, j), checktable(:, 5)));
        assert(numel(i_left) + numel(i_right) == 1);
        if(isVisited([i_left; i_right]))
            % This junction arm has already been visited from its neighbor.
            continue;
        end
        
        if(~isempty(i_left))
            % Start from the left end of dsDNA i_left
            [tour, isVisited, nick_pos, bulge_pos, i_junction, i_arm] = find_dsDNA_tour(dsDNA, checktable, isVisited, i_left, 'left');
        elseif(~isempty(i_right))
            % Start from the right end of dsDNA i_right
            [tour, isVisited, nick_pos, bulge_pos, i_junction, i_arm] = find_dsDNA_tour(dsDNA, checktable, isVisited, i_right, 'right');
        else
            error('Exception.');
        end
        
        % Record the connectivity, incorporate the arm length information into HJE
        assert(isempty(arm_list{i, j}) && isempty(nick_list{i, j}) && isempty(bulge_list{i, j}));
        if(isempty(i_junction))
            arm_list{i, j} = tour;
            nick_list{i, j} = nick_pos;
            bulge_list{i, j} = bulge_pos;
        else
            connectivity = cat(1, connectivity, [i, j, i_junction, i_arm, 0, 0, 0]);
            assert(isempty(arm_list{i_junction, i_arm}) && isempty(nick_list{i_junction, i_arm}) && isempty(bulge_list{i_junction, i_arm}));
            mid_pos = floor(size(tour, 2) / 2 + 1);
            
            arm_list{i, j} = tour(:, 1 : mid_pos);
            nick_list{i, j} = nick_pos(nick_pos < mid_pos);
            bulge_list{i, j} = bulge_pos(bulge_pos < mid_pos);
            
            arm_list{i_junction, i_arm} = flipud(tour(:, end : -1 : mid_pos));
            nick_list{i_junction, i_arm} = size(tour, 2) - fliplr(nick_pos(nick_pos >= mid_pos));
            bulge_list{i_junction, i_arm} = size(tour, 2) - fliplr(bulge_pos(bulge_pos >= mid_pos));
        end
    end
end


%% Build the arms for each HJE
for i = 1 : n_HJE
    % Find the number of base-pairs in each arm
    HJE(i).armLength = zeros(4, 1);
    for j = 1 : 4
        HJE(i).armLength(j) = size(arm_list{i, j}, 2);
    end
    
    % Find the sequence of bases for the current Holliday junction
    HJE(i).idSeq_main = [fliplr(arm_list{i, 1}(2, :)), ...
                                arm_list{i, 2}(1, :), ...
                         fliplr(arm_list{i, 3}(2, :)), ...
                                arm_list{i, 4}(1, :)]';
    HJE(i).idSeq_comp = [fliplr(arm_list{i, 1}(1, :)), ...
                                arm_list{i, 2}(2, :), ...
                         fliplr(arm_list{i, 3}(1, :)), ...
                                arm_list{i, 4}(2, :)]';
end


%% Build the nicks and bulges
for i = 1 : n_HJE
    curr_nick = [HJE(i).armLength(1) - nick_list{i, 1}'; ...
                                       nick_list{i, 2}' + HJE(i).armLength(1); ...
                 HJE(i).armLength(3) - nick_list{i, 3}' + sum(HJE(i).armLength(1:2)); ...
                                       nick_list{i, 4}' + sum(HJE(i).armLength(1:3))];
    curr_bulge = [HJE(i).armLength(1) - bulge_list{i, 1}'; ...
                                        bulge_list{i, 2}' + HJE(i).armLength(1); ...
                  HJE(i).armLength(3) - bulge_list{i, 3}' + sum(HJE(i).armLength(1:2)); ...
                                        bulge_list{i, 4}' + sum(HJE(i).armLength(1:3))];

    nick4FE = cat(1, nick4FE, [i * ones(size(curr_nick)), curr_nick, zeros(size(curr_nick))]);
    bulge4FE = cat(1, bulge4FE, [i * ones(size(curr_bulge)), curr_bulge, zeros(size(curr_bulge))]);
end


%% Find isolated duplexes
i_duplex = 0;
while(~isempty(find(~isVisited, 1)))
    i_duplex = i_duplex + 1;
    i_curr = find(~isVisited, 1);
    tour_l = [];
    nick_pos_l = [];
    bulge_pos_l = [];
    tour_r = [];
    nick_pos_r = [];
    bulge_pos_r = [];
    
    % Start from the right end of the current dsDNA and move to the left
    [tour_l, isVisited, nick_pos_l, bulge_pos_l, i_junction, i_arm] = find_dsDNA_tour(dsDNA, checktable, isVisited, i_curr, 'right');
    assert(isempty(i_junction) && isempty(i_arm));
    
    tour_l = rot90(tour_l, 2);
    nick_pos_l = size(tour_l, 2) - fliplr(nick_pos_l);
    bulge_pos_l = size(tour_l, 2) - fliplr(bulge_pos_l);
    % In case of multiple duplexes forming a ring via nicks
    nick_pos_l(nick_pos_l == 0) = [];
    bulge_pos_l(bulge_pos_l == 0) = [];
    isVisited(i_curr) = false;
    
    % Start from the left end of the current dsDNA and move to the right
    [tour_r, isVisited, nick_pos_r, bulge_pos_r, i_junction, i_arm] = find_dsDNA_tour(dsDNA, checktable, isVisited, i_curr, 'left');
    assert(isempty(i_junction) && isempty(i_arm));
    
    len_curr = numel(dsDNA(i_curr).tour_1);
    tour_r(:, 1 : len_curr) = [];
    nick_pos_r = nick_pos_r - len_curr + size(tour_l, 2);
    bulge_pos_r = bulge_pos_r - len_curr + size(tour_l, 2);
    % In case of multiple duplexes forming a ring via nicks
    nick_pos_r(nick_pos_r == 0) = [];
    bulge_pos_r(bulge_pos_r == 0) = [];
    
    tour = [tour_l, tour_r];
    
    nick_pos = [nick_pos_l, nick_pos_r];
    bulge_pos = [bulge_pos_l, bulge_pos_r];
    nick4FE = cat(1, nick4FE, [(n_HJE + i_duplex) * ones(size(nick_pos')), nick_pos', zeros(size(nick_pos'))]);
    bulge4FE = cat(1, bulge4FE, [(n_HJE + i_duplex) * ones(size(bulge_pos')), bulge_pos', zeros(size(bulge_pos'))]);
    
    HJE(n_HJE + i_duplex).nWay = 2;
    HJE(n_HJE + i_duplex).type = 'duplex';
    
    HJE(n_HJE + i_duplex).armLength = size(tour, 2);
    HJE(n_HJE + i_duplex).idSeq_main = tour(1, :)';
    HJE(n_HJE + i_duplex).idSeq_comp = tour(2, :)';
    
    [c, e1, e2, e3] = duplexTriad(dnaTop, HJE(n_HJE + i_duplex));
    HJE(n_HJE + i_duplex).c = c;
    HJE(n_HJE + i_duplex).e1 = e1;
    HJE(n_HJE + i_duplex).e2 = e2;
    HJE(n_HJE + i_duplex).e3 = e3;
end


%% Build ssDNA connections
% arm_end = cell(numel(HJE), 4);
% for i = 1 : numel(HJE)
%     if(HJE(i).nWay == 4)
%         len_BH = sum(HJE(i).armLength(1:2));
%         arm_end{i, 1} = [HJE(i).idSeq_main(1), HJE(i).idSeq_comp(1)];
%         arm_end{i, 2} = [HJE(i).idSeq_comp(len_BH), HJE(i).idSeq_main(len_BH)];
%         arm_end{i, 3} = [HJE(i).idSeq_main(len_BH + 1), HJE(i).idSeq_comp(len_BH + 1)];
%         arm_end{i, 4} = [HJE(i).idSeq_comp(end), HJE(i).idSeq_main(end)];
%     elseif(HJE(i).nWay == 2)
%         arm_end{i, 1} = [HJE(i).idSeq_main(1), HJE(i).idSeq_comp(1)];
%         arm_end{i, 2} = [HJE(i).idSeq_comp(end), HJE(i).idSeq_main(end)];
%     else
%         error('Exception.');
%     end
% end

conn_ssDNA = [conn_es, conn_gap];
for i = 1 : numel(conn_ssDNA)
    s1 = conn_ssDNA(i).end5(1);
    s2 = conn_ssDNA(i).end5(2);
    t1 = conn_ssDNA(i).end3(1);
    t2 = conn_ssDNA(i).end3(2);
    
    if(s2 == 2)
        e1 = dsDNA(s1).end3_1;
    elseif(s2 == 4)
        e1 = dsDNA(s1).end3_2;
    else
        error('Exception.');
    end
    
    if(t2 == 1)
        e2 = dsDNA(t1).end5_1;
    elseif(t2 == 3)
        e2 = dsDNA(t1).end5_2;
    else
        error('Exception.');
    end
    
%     i_HJE_1 = [];
%     i_arm_1 = [];
%     i_HJE_2 = [];
%     i_arm_2 = [];
%     
%     for i1 = 1 : size(arm_end, 1)
%         for i2 = 1 : size(arm_end, 2)
%             if(~isempty(arm_end{i1, i2}))
%                 if(e1 == arm_end{i1, i2}(2))
%                     i_HJE_1 = cat(2, i_HJE_1, i1);
%                     i_arm_1 = cat(2, i_arm_1, i2);
%                 end
%                 if(e2 == arm_end{i1, i2}(1))
%                     i_HJE_2 = cat(2, i_HJE_2, i1);
%                     i_arm_2 = cat(2, i_arm_2, i2);
%                 end
%             end
%         end
%     end
% 
%     assert(numel(i_HJE_1) == 1 && numel(i_arm_1) == 1 && ...
%            numel(i_HJE_2) == 1 && numel(i_arm_2) == 1);
%     connectivity = cat(1, connectivity, [i_HJE_1, i_arm_1, i_HJE_2, i_arm_2, conn_ssDNA(i).len, 0, 0]);
    
    for j = 1 : numel(HJE)
        j1 = find(HJE(j).idSeq_main == e1 | HJE(j).idSeq_comp == e1);
        if(~isempty(j1))
            assert(numel(j1) <= 2);
            break;
        end
    end
    connSS4FE(i, 1 : 2) = [j, j1(1)];
    
    for j = 1 : numel(HJE)
        j1 = find(HJE(j).idSeq_main == e2 | HJE(j).idSeq_comp == e2);
        if(~isempty(j1))
            assert(numel(j1) <= 2);
            break;
        end
    end
    connSS4FE(i, 3 : 4) = [j, j1(1)];
    
    connSS4FE(i, 5) = conn_ssDNA(i).len;
end


%% Build backbone connections
for i = 1 : numel(junction_non4way)
    nWay = junction_non4way(i).nWay;
    assert(nWay >= 2 && nWay ~= 4);
    tour = junction_non4way(i).tour;
    if(junction_non4way(i).isCircular)
        tour = [tour; tour(1:2, :)];
    end
    assert(mod(size(tour, 1), 2) == 0);
    
    for j = 1 : (size(tour, 1) - 2) / 2
        connBB4FE = cat(1, connBB4FE, zeros(1, 4));
        
        % Find the nucleotide IDs of the nucleotides at both ends of the
        % current backbone connection.
        s1 = tour(j*2, 1);
        s2 = tour(j*2, 2);
        t1 = tour(j*2+1, 1);
        t2 = tour(j*2+1, 2);
        
        if(s2 == 2)
            e1 = dsDNA(s1).end3_1;
        elseif(s2 == 4)
            e1 = dsDNA(s1).end3_2;
        else
            error('Exception.');
        end
        
        if(t2 == 1)
            e2 = dsDNA(t1).end5_1;
        elseif(t2 == 3)
            e2 = dsDNA(t1).end5_2;
        else
            error('Exception.');
        end
        
        % Find the HJE IDs and node IDs of the nucleotides at both ends of the
        % current backbone connection.
        for k = 1 : numel(HJE)
            k1 = find(HJE(k).idSeq_main == e1 | HJE(k).idSeq_comp == e1);
            if(~isempty(k1))
                assert(numel(k1) <= 2);
                break;
            end
        end
        
        if(HJE(k).nWay == 2)
            assert(k1(1) == 1 || k1(1) == sum(HJE(k).armLength));
        elseif(HJE(k).nWay == 4)
            assert(k1(1) == 1 || k1(1) == sum(HJE(k).armLength(1:2)) || ...
                   k1(1) == sum(HJE(k).armLength(1:2)) + 1 || k1(1) == sum(HJE(k).armLength));
        else
            error('Exception.');
        end
        
        connBB4FE(end, 1 : 2) = [k, k1(1)];
        
        for k = 1 : numel(HJE)
            k1 = find(HJE(k).idSeq_main == e2 | HJE(k).idSeq_comp == e2);
            if(~isempty(k1))
                assert(numel(k1) <= 2);
                break;
            end
        end
        
        if(HJE(k).nWay == 2)
            assert(k1(1) == 1 || k1(1) == sum(HJE(k).armLength));
        elseif(HJE(k).nWay == 4)
            assert(k1(1) == 1 || k1(1) == sum(HJE(k).armLength(1:2)) || ...
                   k1(1) == sum(HJE(k).armLength(1:2)) + 1 || k1(1) == sum(HJE(k).armLength));
        else
            error('Exception.');
        end
        
        connBB4FE(end, 3 : 4) = [k, k1(1)];
    end
end


%% Final adjustment
% For each connection, identify its polarity
for i = 1 : size(connectivity, 1)
    if(connectivity(i, 5) > 0 || connectivity(i, 6) > 0)  % connected by ssDNA
        connectivity(i, 7) = Inf;
    else
        connectivity(i, 7) = 1 - mod(connectivity(i, 2) + connectivity(i, 4), 2);
    end
end

% Assign each nt as main/comp
strand4FE = strand;
% Initialize the field isMain
for i = 1 : numel(strand4FE)
    strand4FE(i).isMain = NaN(1, numel(strand4FE(i).tour));
end
% Calculate for the field isMain
for i = 1 : numel(HJE)
    for j = 1 : numel(HJE(i).idSeq_main)
        iStrand = dnaTop(HJE(i).idSeq_main(j)).strand;
        iResidue = dnaTop(HJE(i).idSeq_main(j)).residue;
        assert(~isinf(strand4FE(iStrand).isMain(iResidue)));
        if(isnan(strand4FE(iStrand).isMain(iResidue)))
            strand4FE(iStrand).isMain(iResidue) = true;
        else
            strand4FE(iStrand).isMain(iResidue) = Inf;  % some bases could be in two HJE's
        end
    end
    
    for j = 1 : numel(HJE(i).idSeq_comp)
        iStrand = dnaTop(HJE(i).idSeq_comp(j)).strand;
        iResidue = dnaTop(HJE(i).idSeq_comp(j)).residue;
        assert(~isinf(strand4FE(iStrand).isMain(iResidue)));
        if(isnan(strand4FE(iStrand).isMain(iResidue)))
            strand4FE(iStrand).isMain(iResidue) = false;
        else
            strand4FE(iStrand).isMain(iResidue) = Inf;  % some bases could be in two HJE's
        end
    end
end

% Check if the nicks are stacked.
% nick4FE(i,3)==0 ---> stacked
% nick4FE(i,3)==1 ---> non-stacked, connected at idSeq_main
% nick4FE(i,3)==2 ---> non-stacked, connected at idSeq_comp
thres_nickStack = 3;    % Unit: nm
for i = 1 : size(nick4FE, 1)
    i_HJE = nick4FE(i, 1);
    i_pos = nick4FE(i, 2);
    i_pos_next = mod(i_pos, size(HJE(i_HJE).idSeq_main, 1)) + 1;
    v1 = dnaTop(HJE(i_HJE).idSeq_main(i_pos_next)).xyz - dnaTop(HJE(i_HJE).idSeq_main(i_pos)).xyz;
    v2 = dnaTop(HJE(i_HJE).idSeq_comp(i_pos_next)).xyz - dnaTop(HJE(i_HJE).idSeq_comp(i_pos)).xyz;
    if(norm(v1)>thres_nickStack || norm(v2)>thres_nickStack)
        % non-stacked nick
        
        if(dnaTop(HJE(i_HJE).idSeq_comp(i_pos)).up ~= HJE(i_HJE).idSeq_comp(i_pos_next))
            % connected as the idSeq_main
            nick4FE(i, 3) = 1;
        elseif(dnaTop(HJE(i_HJE).idSeq_main(i_pos)).down ~= HJE(i_HJE).idSeq_main(i_pos_next))
            % connected as the idSeq_comp
            nick4FE(i, 3) = 2;
        else
            error('Exception.');
        end
    end
end

end


% This function finds the tour from one junction arm to the neighboring
% junction arm, or a dsDNA end.
function [tour, isVisited, nick_pos, bulge_pos, i_junction, i_arm] = find_dsDNA_tour(dsDNA, checktable, isVisited, i_dsDNA, i_polarity)

%     tour
%     Row 1: tour start -------------------> tour end
%     Row 2: tour start <------------------- tour end
tour = zeros(2, 0);
nick_pos = [];
bulge_pos = [];
isEnd = false;
prev_type = [];

while(~isEnd)
    curr_ds = dsDNA(i_dsDNA);
    assert(~isVisited(i_dsDNA));
    isVisited(i_dsDNA) = true;
    
    if(strcmp(i_polarity, 'left'))
        % Enter from the left side of the dsDNA, exit on the right side
        i_start = 1;
        i_end = 4;
        tour = cat(2, tour, [curr_ds.tour_1, flipud(curr_ds.tour_2)]');
    elseif(strcmp(i_polarity, 'right'))
        % Enter from the right side of the dsDNA, exit on the left side
        i_start = 4;
        i_end = 1;
        tour = cat(2, tour, [curr_ds.tour_2, flipud(curr_ds.tour_1)]');
    else
        error('Exception.');
    end
    
    if(~isempty(prev_type))
        assert(strcmp(checktable{i_dsDNA, i_start + 2}, prev_type));
    end

    if(isempty(checktable{i_dsDNA, i_end + 2}))
        % End of a dsDNA without connecting to another dsDNA
        i_junction = [];
        i_arm = [];
        isEnd = true;
    elseif(strcmp(checktable{i_dsDNA, i_end + 2}, 'junction'))
        % Find another junction
        i_junction = checktable{i_dsDNA, i_end};
        i_arm = checktable{i_dsDNA, i_end + 1};
        assert(isnumeric(i_junction) && isnumeric(i_arm));
        isEnd = true;
    elseif(strcmp(checktable{i_dsDNA, i_end + 2}, 'nick'))
        % Find a nick
        i_polarity = checktable{i_dsDNA, i_end + 1};
        assert(ischar(i_polarity));
        i_dsDNA = checktable{i_dsDNA, i_end};
        prev_type = 'nick';
        nick_pos = cat(2, nick_pos, size(tour, 2));
        i_junction = [];
        i_arm = [];
        isEnd = isVisited(i_dsDNA);
    elseif(strcmp(checktable{i_dsDNA, i_end + 2}, 'bulge'))
        % Find a bulge
        i_polarity = checktable{i_dsDNA, i_end + 1};
        assert(ischar(i_polarity));
        i_dsDNA = checktable{i_dsDNA, i_end};
        prev_type = 'bulge';
        bulge_pos = cat(2, bulge_pos, size(tour, 2));
        i_junction = [];
        i_arm = [];
        isEnd = isVisited(i_dsDNA);
    else
        error('Exception.');
    end
end

end


% This function helps to find a number of string in a cell array
function x = find_cell(val, key)

if(isempty(val) || isempty(key))
    x = false;
elseif(isnumeric(val) && isnumeric(key))
    x = (val == key);
elseif(ischar(val) && ischar(key))
    x = strcmp(val, key);
elseif((isnumeric(val) && ischar(key)) || (ischar(val) && isnumeric(key)))
    x = false;
else
    error('Exception.');
end

end


% This function build three axes e1, e2, and e3 for a double-stranded or
% single-stranded crossover
function [c, e1, e2, e3] = junctionTriad(dnaTop, dsDNA, junction)

if(numel(junction) ~= 1)
    error('Multiple junctions');
end
if(junction.nWay~=4)
    error('Not a 4-way junction');
end

t = zeros(4, 1);
for i = 1 : numel(t)
    if(junction.tour(i * 2, 2) == 2)
        t(i) = dsDNA(junction.tour(i * 2, 1)).end3_1;
    elseif(junction.tour(i * 2, 2) == 4)
        t(i) = dsDNA(junction.tour(i * 2, 1)).end3_2;
    else
        error('Exception.');
    end
end

% Position of 8 bases at the branch point
b1 = dnaTop(t(1)).xyz;
b2 = dnaTop(t(2)).xyz;
b3 = dnaTop(t(3)).xyz;
b4 = dnaTop(t(4)).xyz;
b1across = dnaTop(dnaTop(t(1)).across).xyz;
b2across = dnaTop(dnaTop(t(2)).across).xyz;
b3across = dnaTop(dnaTop(t(3)).across).xyz;
b4across = dnaTop(dnaTop(t(4)).across).xyz;

% Junction center
c = mean([b1, b2, b3, b4, b1across, b2across, b3across, b4across], 2);

% Vectors from bases to their complementary bases
bpv1 = b1across - b1;
bpv2 = b2across - b2;
bpv3 = b3across - b3;
bpv4 = b4across - b4;
assert(norm(bpv1)>0 && norm(bpv2)>0 && norm(bpv3)>0 && norm(bpv4)>0);

% Axis e1 (axis of arm 1-2, pointing to arm 2)
e1 = cross(bpv1, bpv2);
% if(norm(e1) <= 0)
if(norm(e1) <= 0.05)
    e1 = (b2 + b2across) / 2 - (b1 + b1across) / 2;
end
assert(norm(e1) > 0 && norm(b2-b1) > 0);
if(e1'*(b2-b1) < 0)
    e1 = -e1;
end
e1 = e1/norm(e1);

% Axis e2
a12 = mean([b1, b2, b1across, b2across], 2);
a34 = mean([b3, b4, b3across, b4across], 2);
e3b = a34 - a12;
e2 = cross(e3b, e1);
assert(norm(e2) > 0);
e2 = e2/norm(e2);

% Axis e3
e3 = cross(e1, e2);
assert(norm(e3) > 0);
e3 = e3/norm(e3);

end


% This function build three axes e1, e2, and e3 for a DNA duplex
function [c, e1, e2, e3] = duplexTriad(dnaTop, duplex)

if(numel(duplex) ~= 1)
    error('Multiple duplexes');
end

% Find the triad center
b1 = dnaTop(duplex.idSeq_main(1)).xyz;
b1across = dnaTop(duplex.idSeq_comp(1)).xyz;
c = mean([b1, b1across], 2);

if(numel(duplex.idSeq_main) > 1)
    b2 = dnaTop(duplex.idSeq_main(2)).xyz;
    b2across = dnaTop(duplex.idSeq_comp(2)).xyz;
    
    % Vectors from bases to their complementary bases
    bpv1 = b1across - b1;
    bpv2 = b2across - b2;
    assert(norm(bpv1)>0 && norm(bpv2)>0);
    
    % Axis e1 (axis of the duplex, pointing to the 3'-end of the main strand)
    e1 = cross(bpv1, bpv2);
    assert(norm(e1)>0 && norm(b2-b1)>0);
    if(e1'*(b2-b1) < 0)
        e1 = -e1;
    end
    e1 = e1/norm(e1);
    
    % Axis e3
    e3 = bpv1/norm(bpv1);
    
    % Axis e2
    e2 = cross(e3, e1);
    assert(norm(e2) > 0);
    e2 = e2/norm(e2);
else
    e1 = [1,0,0]';
    e2 = [0,1,0]';
    e3 = [0,0,1]';
end

end