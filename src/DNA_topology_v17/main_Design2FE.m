function [] = main_Design2FE(DesignPATH, DesignType, latticeType, outputDIR, seqPATH)

% If a ssDNA contains >= N_ss nts, then model it as a truss
N_ss = 5;   % Used in the on-lattice CanDo


[~,outputNAME,~] = fileparts(DesignPATH);
assert(~isempty(outputNAME));
outputPATH = fullfile(outputDIR,outputNAME);

% Create a log file
fid_log = fopen(strcat(outputPATH,'.log'), 'w');


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1. Create the directed graph model for the DNA assembly
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write the log file
fprintf(fid_log, '%s\n', repmat('#',1,80));
fprintf(fid_log, '# Step 1. Create the directed graph model for the DNA assembly\n');
fprintf(fid_log, '%s\n', repmat('#',1,80));

%% Step 1.1 Read the directed graph model from a design file in Tiamat or
% caDNAno format
if strcmp(DesignType,'tiamat')
    % Tiamat design file
    dnaTop = Tiamat2topology(DesignPATH);      % get the topology information
    dnaTop = modify_dnaTop(dnaTop);                 % swap 'up' and 'down' in Tiamat
elseif strcmp(DesignType,'cndo')
	% CNDO design file
	dnaInfo = cndo2dnaInfo(DesignPATH);
	dnaTop = modify_dnaTop(dnaInfo.dnaTop);
else
    error('Unrecognizable design type.');
end

% Write the log file
fprintf(fid_log, '\n[STEP 1.1] Read the directed graph model from a design file in %s format...\n', DesignType);
fprintf(fid_log, '\tNumber of nucleotides = %d\n', numel(dnaTop));


%% Step 1.2 Organize graph nodes (nucleotides) into DNA strands
[dnaTop, strand] = buildStrand(dnaTop);         % get strand information, routing the whole design

% Write the log file
fprintf(fid_log, '\n[STEP 1.2] Organize graph nodes (nucleotides) into DNA strands...\n');
fprintf(fid_log, '\tNumber of strands = %d\n', numel(strand));

% For CNDO design only
if(strcmp(DesignType,'cndo'))
	% Add a field to dnaTop as the Cartesian coordinate of each nucleotide
	for i = 1 : numel(dnaTop)
		dnaTop(i).xyz = zeros(3,1);
	end
	% Compute the Cartesian coordinate of each nucleotide
	[xyz_prefer, xyz_alt] = nt_coord(dnaInfo.dnaGeom.dNode, dnaInfo.dnaGeom.triad);
	for i = 1 : size(dnaInfo.dnaGeom.id_nt,1)
		i_prefer = dnaInfo.dnaGeom.id_nt(i,1);
		i_alt = dnaInfo.dnaGeom.id_nt(i,2);
		% Convert the unit of Cartesian coordinates from Angstrom (in dnaInfo) to nm (in dnaTop)
		dnaTop(i_prefer).xyz = xyz_prefer(i,:)' / 10;
		dnaTop(i_alt).xyz = xyz_alt(i,:)' / 10;
	end
end

%% Step 1.3 Assign each graph node (nucleotide) as A, T, G, or C.
%% Tiamat design only, CNDO design already includes sequence 
if(strcmp(DesignType,'tiamat'))
	if(~isempty(seqPATH))
		%dnaTop = assignSeq(dnaTop, strand, seqPATH);    % assign the pre-defined sequences
		dnaTop = assignSeqFromTiamat(dnaTop, strand, seqPATH);    % assign the pre-defined sequences
	else
		dnaTop = assignSeqRand(dnaTop);                 % randomly assign the sequences
	end
end

% Write the log file
fprintf(fid_log, '\n[STEP 1.3] Assign each graph node (nucleotide) as A, T, G, or C...\n');
if(~isempty(seqPATH))
    fprintf(fid_log, '\tSequences of strands:\t%s\n', seqPATH);
else
    fprintf(fid_log, '\tNo sequence information is available.\n\tRandomly generate nucleotide sequences for all strands.\n');
end
    

%% Step 1.4 Visualization
generateBILD(dnaTop, strand, [], [], sprintf('%s_strand.bild',outputPATH), 'strand');   % colored according to strand IDs
generateBILD(dnaTop, strand, [], [], sprintf('%s_sequence.bild',outputPATH), 'sequence');   % colored according to ATGC

% Write the log file
fprintf(fid_log, '\n[STEP 1.4] Visualize the topology and sequence information...\n');
fprintf(fid_log, '\tGenerate a .bild file: %s_strand.bild\n', outputPATH);
fprintf(fid_log, '\tGenerate a .bild file: %s_sequence.bild\n', outputPATH);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 2. Divide the directed graph model into:
%         (1) dsDNA, (2) ssDNA,
%         and find the connectivity between them
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write the log file
fprintf(fid_log, '\n\n\n%s\n', repmat('#',1,80));
fprintf(fid_log, '# Step 2. Divide the directed graph model into:\n');
fprintf(fid_log, '#         (1) dsDNA segments, (2) ssDNA segments,\n');
fprintf(fid_log, '#         and find the connectivity between them\n');
fprintf(fid_log, '%s\n', repmat('#',1,80));

%% Step 2.1 Divide the DNA assembly into dsDNA and ssDNA
% Definition of ends
%
%                    dsDNA
%       (4)                        (3)
%     end3_2 <------------------- end5_2
%     end5_1 -------------------> end3_1
%       (1)                        (2)
%                    ssDNA
%       end5 -------------------> end3
%        (1)                      (2)
%
[dsDNA, ssDNA] = divideTopology(dnaTop);    % divide the directed graph into dsDNA and ssDNA

% Write the log file
fprintf(fid_log, '\n[STEP 2.1] Find all dsDNA segment(s) and ssDNA segment(s) in the directed graph model\n');
fprintf(fid_log, '\tNumber of dsDNA segment(s) = %d\n', numel(dsDNA));
fprintf(fid_log, '\tNumber of ssDNA segment(s) = %d\n', numel(ssDNA));
fprintf(fid_log, '\n\t%s\n', repmat('-',1,40));
for i = 1:numel(dsDNA);
    assert(numel(dsDNA(i).tour_1) == numel(dsDNA(i).tour_2));
    fprintf(fid_log, '\n\tdsDNA %d (%d basepairs)\n', i, numel(dsDNA(i).tour_1));
    fprintf(fid_log, '\tStrand 1%6d -------------------> %-6d\n', dsDNA(i).end5_1, dsDNA(i).end3_1);
    fprintf(fid_log, '\tStrand 2%6d <------------------- %-6d\n', dsDNA(i).end3_2, dsDNA(i).end5_2);
end
for i = 1:numel(ssDNA);
    fprintf(fid_log, '\n\tssDNA %d (%d basepairs)\n', i, numel(ssDNA(i).tour));
    fprintf(fid_log, '\tStrand 1%6d -------------------> %-6d\n', ssDNA(i).end5, ssDNA(i).end3);
end
fprintf(fid_log, '\n\t%s\n', repmat('-',1,40));


%% Step 2.2 Find the connectivity between the dsDNA and ssDNA
conn_dsDNA_ssDNA = findConn(dnaTop, dsDNA, ssDNA);  % find connections between these dsDNA and ssDNA

% Write the log file
fprintf(fid_log, '\n[STEP 2.2] Find the phosphodiester bonds that connect the dsDNA/ssDNA segments\n');
fprintf(fid_log, '\tNumber of such phosphodiester bonds = %d\n', size(conn_dsDNA_ssDNA,1));
log_print_pd_bond(fid_log, conn_dsDNA_ssDNA);


%% Step 2.3 Visualization
generateBILD(dnaTop, strand, conn_dsDNA_ssDNA, [], sprintf('%s_conn_dsDNA_ssDNA.bild',outputPATH), 'conn_dsDNA_ssDNA');

% Write the log file
fprintf(fid_log, '\n[STEP 2.3] Visualize the phosphodiester bonds that connect the dsDNA/ssDNA segments...\n');
fprintf(fid_log, '\tGenerate a .bild file: %s_conn_dsDNA_ssDNA.bild\n', outputPATH);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 3. Cluster the connections between dsDNA/ssDNA segments
%         into (1) junction, (2) nick, and (3) unpaired regions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write the log file
fprintf(fid_log, '\n\n\n%s\n', repmat('#',1,80));
fprintf(fid_log, '# Step 3. Cluster the connections between dsDNA/ssDNA segments\n');
fprintf(fid_log, '#         into junctions and nicks.\n');
fprintf(fid_log, '%s\n', repmat('#',1,80));

%% Step 3.1 Find topological motifs consisting of dsDNAs by clustering the 
% phosphodiester bond connections between dsDNA segments into
% - Nicks (or equivalently, 2-way junctions)
% - N-way junctions (N = 3, 4, 5, ...)
% junction.tour: (nWay*2)-by-2 matrix, each row is (dsDNA ID, end ID)
% nick.tour: 4-by-2 matrix, each row is (dsDNA ID, end ID)
[junction, nick, conn_cluster] = find_junction_nick(conn_dsDNA_ssDNA);

% Write the log file
fprintf(fid_log, '\n[STEP 3.1] Find topological motifs consisting of dsDNAs by clustering the\n');
fprintf(fid_log, '\tphosphodiester bond connections between dsDNA segments into\n');
fprintf(fid_log, '\t(1) nicks (or equivalently, 2-way junctions) and (2) n-way junctions (n = 3, 4, 5, ...)\n');
fprintf(fid_log, '\n\t%s\n', repmat('-',1,40));
fprintf(fid_log, '\tNumber of n-way junctions = %d\n', numel(junction));
fprintf(fid_log, '\tNumber of nicks = %d\n', numel(nick));
fprintf(fid_log, '\t%s\n', repmat('-',1,40));


%% Step 3.2 Find connections by ssDNA, which could be
% - ssDNA connecting two dsDNA (unpaired_conn)
% - ssDNA with a free end (unpaired_free)
[unpaired_conn, unpaired_free] = find_unpaired(conn_dsDNA_ssDNA);

% Write the log file
fprintf(fid_log, '\n[STEP 3.2] Find connections by ssDNA, which could be\n');
fprintf(fid_log, '\t(1) ssDNA connecting two dsDNA, or (2) ssDNA with a free end\n');
fprintf(fid_log, '\n\tTotal number of ssDNA segments = %d\n', numel(ssDNA));
fprintf(fid_log, '\n\t%s\n', repmat('-',1,40));
fprintf(fid_log, '\tNumber of ssDNA connecting two dsDNA = %d\n', numel(unpaired_conn));
fprintf(fid_log, '\tNumber of ssDNA with a free end = %d\n', numel(unpaired_free));
fprintf(fid_log, '\t%s\n', repmat('-',1,40));


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 4. Create topology motifs
%         - Double-stranded crossovers (4-way junctions)
%         - Single-stranded crossovers (4-way junctions)
%         - Nicks
%         - Single-stranded regions (ssDNA with length >= N_ss)
%         - Gaps (ssDNA with length < N_ss)
%         - Bulges (nick + gap)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write the log file
fprintf(fid_log, '\n\n\n%s\n', repmat('#',1,80));
fprintf(fid_log, '# Step 4. Create topology motifs\n');
fprintf(fid_log, '#         - Double-stranded crossovers (4-way junctions)\n');
fprintf(fid_log, '#         - Single-stranded crossovers (4-way junctions)\n');
fprintf(fid_log, '#         - Nicks\n');
fprintf(fid_log, '#         - Single-stranded region (ssDNA with length >= N_ss)\n');
fprintf(fid_log, '#         - Gaps (ssDNA with length < N_ss)\n');
fprintf(fid_log, '#         - Bulges (nick + ssDNA)\n');
fprintf(fid_log, '%s\n', repmat('#',1,80));


%% Step 4.1 Categorize all junctions into:
%      (1) 4-way junctions as double-stranded crossovers,
%      (2) 4-way junctions as single-stranded crossovers,
%      (3) non-4-way junctions
n_junction = numel(junction);
i_dsxover = [];
i_ssxover = [];
i_non4way = [];
for i = 1:n_junction
    if(junction(i).nWay == 4)
        if(junction(i).isCircular == 1)
            i_dsxover = cat(2, i_dsxover, i);
        elseif(junction(i).isCircular == 0)
            i_ssxover = cat(2, i_ssxover, i);
        else
            error('Exception.');
        end
    else
        i_non4way = cat(2, i_non4way, i);
    end
end

% Write the log file
fprintf(fid_log, '\n[STEP 4.1] Categorize all junctions into:\n');
fprintf(fid_log, '\t(1) 4-way junctions as double-stranded crossovers\n');
fprintf(fid_log, '\t(2) 4-way junctions as single-stranded crossovers\n');
fprintf(fid_log, '\t(3) non-4-way junctions\n');

fprintf(fid_log, '\n\t(1) Number of 4-way junctions as double-stranded crossover: %d\n', numel(i_dsxover));
log_print_junction_nick(fid_log, junction, i_dsxover);
fprintf(fid_log, '\n\t(2) Number of 4-way junctions as single-stranded crossover: %d\n', numel(i_ssxover));
log_print_junction_nick(fid_log, junction, i_ssxover);
fprintf(fid_log, '\n\t(3) Number of non-4-way junctions: %d\n', numel(i_non4way));
log_print_junction_nick(fid_log, junction, i_non4way);


%% Step 4.2 Re-categorize (1) all previously identified nicks and (2) ssDNA connecting two dsDNA into:
% (1) real nicks (nick_real)
% (2) bulges
% (3) Single-stranded regions (ssDNA with length >= N_ss)
% (4) Gaps (ssDNA with length < N_ss)
[nick_real, bulge, conn_es, conn_gap] = find_bulge_es_gap(nick, unpaired_conn, ssDNA, N_ss);

% Write the log file
fprintf(fid_log, '\n[STEP 4.2] Categorize all previously identified nicks and ssDNAs into:\n');
fprintf(fid_log, '\t(1) real nicks\n');
fprintf(fid_log, '\t(2) bulges\n');
fprintf(fid_log, '\t(3) single-stranded regions (ssDNA with length >= N_ss)\n');
fprintf(fid_log, '\t(4) gaps (ssDNA with length < N_ss)\n');
fprintf(fid_log, '\t(Note: N_ss = %d)\n', N_ss);

fprintf(fid_log, '\n\t(1) Number of real nicks: %d\n', numel(nick_real));
log_print_junction_nick(fid_log, nick_real, 1 : numel(nick_real));
fprintf(fid_log, '\n\t(2) Number of bulges: %d\n', numel(bulge));
log_print_junction_nick(fid_log, bulge, 1 : numel(bulge));
fprintf(fid_log, '\n\t(3) Number of single-stranded regions: %d\n', numel(conn_es));
log_print_conn_es_gap(fid_log, conn_es);
fprintf(fid_log, '\n\t(4) Number of gaps: %d\n', numel(conn_gap));
log_print_conn_es_gap(fid_log, conn_gap);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 5. Create structural motifs for the finite element calculation
%         - Double-stranded crossovers (4-way junctions)
%         - Single-stranded crossovers (4-way junctions)
%         - Nicks
%         - Single-stranded regions (ssDNA with length >= N_ss)
%         - Gaps (ssDNA with length < N_ss, or bulges)
%         - Free duplexes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write the log file
fprintf(fid_log, '\n\n\n%s\n', repmat('#',1,80));
fprintf(fid_log, '# Step 5. Create structural motifs for the finite element calculation\n');
fprintf(fid_log, '#         - Double-stranded crossovers (4-way junctions)\n');
fprintf(fid_log, '#         - Single-stranded crossovers (4-way junctions)\n');
fprintf(fid_log, '#         - Nicks\n');
fprintf(fid_log, '#         - Single-stranded region (ssDNA with length >= N_ss)\n');
fprintf(fid_log, '#         - Gaps (ssDNA with length < N_ss, or bulges)\n');
fprintf(fid_log, '#         - Free non-circular duplexes\n');
fprintf(fid_log, '%s\n', repmat('#',1,80));


%% Step 5.1 Search for topological motifs that will not be modeled in the finite element calculation
% Write the log file
fprintf(fid_log, '\n[STEP 5.1] Search for topological motifs that will not be modeled in the finite element calculation:\n');
fprintf(fid_log, '\t(1) Number of non-4-way junctions = %d\n', numel(i_non4way));
% Fine the number of circular DNA duplexes
n_circular_dsDNA = 0;
for i = 1 : numel(dsDNA)
    if(dsDNA(i).isCircular)
        n_circular_dsDNA = n_circular_dsDNA + 1;
    end
end
fprintf(fid_log, '\t(2) Number of circular DNA duplexes = %d\n', n_circular_dsDNA);
fprintf(fid_log, '\t[WARNING] Circular DNA duplexes are ignored in the downstream steps.\n');

% if(numel(i_non4way) > 0 || n_circular_dsDNA > 0)
%     fprintf(fid_log, '\tNon-4-way junction(s) or circular DNA duplex(es) are found. The program is going to terminate.\n');
%     error('Non-4-way junction(s) or circular DNA duplex(es) are found.');
% end


%% Step 5.2 Determine the isomeric states of all double-stranded crossovers
% Write the log file
fprintf(fid_log, '\n[STEP 5.2] Determine the isomeric states of all double-stranded crossovers:\n');
fprintf(fid_log, '\n\t%s\n', repmat('-',1,40));
for i = i_dsxover
    score_isoI = valStack(dnaTop, dsDNA, junction(i).tour(1,:), junction(i).tour(3,:)) + ...
                 valStack(dnaTop, dsDNA, junction(i).tour(5,:), junction(i).tour(7,:));
    score_isoII = valStack(dnaTop, dsDNA, junction(i).tour(3,:), junction(i).tour(5,:)) + ...
                  valStack(dnaTop, dsDNA, junction(i).tour(7,:), junction(i).tour(1,:));
    fprintf(fid_log, '\n\tJunction %d as a double-stranded crossover: iso-I score = %f, iso-II score = %f\n', i, score_isoI, score_isoII);
    if(score_isoI < score_isoII)
        fprintf(fid_log, '\tAdjust the isomeric form.\n');
        junction(i).tour = junction(i).tour([3:8 1:2], :);
        score_isoI = valStack(dnaTop, dsDNA, junction(i).tour(1,:), junction(i).tour(3,:)) + ...
                     valStack(dnaTop, dsDNA, junction(i).tour(5,:), junction(i).tour(7,:));
        score_isoII = valStack(dnaTop, dsDNA, junction(i).tour(3,:), junction(i).tour(5,:)) + ...
                      valStack(dnaTop, dsDNA, junction(i).tour(7,:), junction(i).tour(1,:));
        fprintf(fid_log, '\tJunction %d as a double-stranded crossover: iso-I score = %f, iso-II score = %f\n', i, score_isoI, score_isoII);
    end
end
fprintf(fid_log, '\n\t%s\n', repmat('-',1,40));


%% Step 5.3 Create a checktable showing the connectivity of the dsDNAs via 
%           - nicks
%           - bulges
%           - single-stranded regions
%           - single-strained and double-strained crossovers

% Write the log file
fprintf(fid_log, '\n[STEP 5.3] Create a checktable showing the connectivity of the dsDNAs via:\n');
fprintf(fid_log, '\t- nicks\n');
fprintf(fid_log, '\t- bulges\n');
fprintf(fid_log, '\t- single-stranded regions\n');
fprintf(fid_log, '\t- single-strained and double-strained crossovers\n');

% Convert bulges to 2-way junctions
n_junction = numel(junction);
n_bulge = numel(bulge);
for i = 1 : n_bulge
    junction(n_junction + i).isCircular = false;
    junction(n_junction + i).nWay = 2;
    junction(n_junction + i).tour = bulge(i).tour;
end
i_non4way = cat(2, i_non4way, n_junction + (1:n_bulge));
bulge = [];

% Create the checktable
checktable = find_checktable(nick_real, bulge, junction([i_dsxover, i_ssxover]), numel(dsDNA));

% Write the log file
fprintf(fid_log, '\n\tChecktable\n');
fprintf(fid_log, '\t%s\n', repmat('-',1,40));
for i = 1 : size(checktable, 1)
    fprintf(fid_log, '\tdsDNA %d |', i);
    if(strcmp(checktable{i,3}, 'junction'))
        fprintf(fid_log, '\t%d\t%d\t%s\t |', checktable{i,1}, checktable{i,2}, checktable{i,3});
    else
        fprintf(fid_log, '\t%d\t%s\t%s\t\t |', checktable{i,1}, checktable{i,2}, checktable{i,3});
    end
    if(strcmp(checktable{i,6}, 'junction'))
        fprintf(fid_log, '\t%d\t%d\t%s\n', checktable{i,4}, checktable{i,5}, checktable{i,6});
    else
        fprintf(fid_log, '\t%d\t%s\t%s\n', checktable{i,4}, checktable{i,5}, checktable{i,6});
    end
end
fprintf(fid_log, '\t%s\n', repmat('-',1,40));


%% Step 5.4 Create the data structures for the finite element solver
% Write the log file
fprintf(fid_log, '\n[STEP 5.4] Create the data structures for the finite element solver:\n');

[HJE, connectivity, strand4FE, nick4FE, bulge4FE, connSS4FE, connBB4FE] ...
	= FE_interface(dnaTop, strand, dsDNA, conn_es, conn_gap, ...
                   junction([i_dsxover, i_ssxover]), junction(i_non4way), checktable);

strand = strand4FE;
nick = nick4FE;
bulge = bulge4FE;
connSSE = connSS4FE;
connBBE = connBB4FE;
% fprintf(fid_log, '\n\tNumber of finite element motifs for double-stranded crossovers = %d\n', n_dsHJE);
% fprintf(fid_log, '\tNumber of finite element motifs for double-stranded crossovers = %d\n', n_ssHJE);


%% Step 5.5 Visualization
generateBILD(dnaTop, strand, [], HJE, sprintf('%s_HJE.bild',outputPATH), 'HJE');   % colored according to ATGC
fprintf(fid_log, '\n[STEP 5.5] Visualize the topology and sequence information...\n');
fprintf(fid_log, '\tGenerate a .bild file: %s_HJE.bild\n', outputPATH);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Final step. Save the results and clean up
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write the log file
fprintf(fid_log, '\n\n\n%s\n', repmat('#',1,80));
fprintf(fid_log, '# Final step. Save the results and clean up\n');
fprintf(fid_log, '%s\n', repmat('#',1,80));

% Save the results
save(strcat(outputPATH, '.mat'), 'dnaTop', 'strand', 'HJE', 'connectivity', 'nick', 'bulge', 'connSSE', 'connBBE', 'n_bulge');
% save(strcat(outputPATH, '.mat'));
fprintf(fid_log, '\n\tSave the results to %s.mat\n', outputPATH);

% Close the log file and terminate the parser
fprintf(fid_log, '\n[The parser has terminated successfully.]\n');
fclose(fid_log);

end



function [] = log_print_pd_bond(fid_log, conn_dsDNA_ssDNA)

fprintf(fid_log, '\n\t%s\n', repmat('-',1,40));
for i = 1:size(conn_dsDNA_ssDNA,1)
    if(conn_dsDNA_ssDNA(i,1) == 1)
        s11 = 'dsDNA';
        if(conn_dsDNA_ssDNA(i,3) == 2)
            s12 = 'strand 1';
        elseif(conn_dsDNA_ssDNA(i,3) == 4)
            s12 = 'strand 2';
        else
            error('Exception.');
        end
    elseif(conn_dsDNA_ssDNA(i,1) == 0)
        s11 = 'ssDNA';
        if(conn_dsDNA_ssDNA(i,3) == 2)
            s12 = 'strand 1';
        else
            error('Exception.');
        end
    else
        error('Exception.');
    end
    if(conn_dsDNA_ssDNA(i,5) == 1)
        s21 = 'dsDNA';
        if(conn_dsDNA_ssDNA(i,7) == 1)
            s22 = 'strand 1';
        elseif(conn_dsDNA_ssDNA(i,7) == 3)
            s22 = 'strand 2';
        else
            error('Exception.');
        end
    elseif(conn_dsDNA_ssDNA(i,5) == 0)
        s21 = 'ssDNA';
        if(conn_dsDNA_ssDNA(i,7) == 1)
            s22 = 'strand 1';
        else
            error('Exception.');
        end
    else
        error('Exception.');
    end
    fprintf(fid_log, '\n\t%s %-6d | %s | 3''-end ---------> %s %-6d | %s | 5''-end\n', ...
        s11, conn_dsDNA_ssDNA(i,2), s12, s21, conn_dsDNA_ssDNA(i,6), s22);
end
fprintf(fid_log, '\n\t%s\n', repmat('-',1,40));

end



function [] = log_print_junction_nick(fid_log, junction, ind)

fprintf(fid_log, '\t%s\n', repmat('-',1,40));
for i = ind
    nWay = size(junction(i).tour, 1) / 2;
    assert(mod(nWay, 1) == 0);
    for j = 1 : nWay;
        j1 = j * 2 - 1;
        j2 = j * 2;
        assert(junction(i).tour(j1,1) == junction(i).tour(j2,1));
        fprintf(fid_log, '\tdsDNA %d ', junction(i).tour(j1,1));
        if(junction(i).tour(j1,2) == 1 && junction(i).tour(j2,2) == 4)
            fprintf(fid_log, 'left end ');
        elseif(junction(i).tour(j1,2) == 3 && junction(i).tour(j2,2) == 2)
            fprintf(fid_log, 'right end ');
        else
            error('Exception.');
        end
        if(j < nWay)
            fprintf(fid_log, '---> ');
        else
            fprintf(fid_log, '\n');
        end
    end
end
fprintf(fid_log, '\t%s\n', repmat('-',1,40));

end



function [] = log_print_conn_es_gap(fid_log, conn_es_gap)

s11 = 'dsDNA';
s21 = 'dsDNA';

fprintf(fid_log, '\t%s\n', repmat('-',1,40));
for i = 1 : numel(conn_es_gap)
    if(conn_es_gap(i).end5(2) == 2)
        s12 = 'strand 1';
    elseif(conn_es_gap(i).end5(2) == 4)
        s12 = 'strand 2';
    else
        error('Exception.');
    end
    if(conn_es_gap(i).end3(2) == 1)
        s22 = 'strand 1';
    elseif(conn_es_gap(i).end3(2) == 3)
        s22 = 'strand 2';
    else
        error('Exception.');
    end
    fprintf(fid_log, '\t%s %-6d | %s | 3''-end ---------> %s %-6d | %s | 5''-end\n', ...
        s11, conn_es_gap(i).end5(1), s12, s21, conn_es_gap(i).end3(1), s22);
end
fprintf(fid_log, '\t%s\n', repmat('-',1,40));

end



function vStack = valStack(dnaTop, dsDNA, tour_1, tour_2)

i1 = tour_1(1);
i2 = tour_1(2);
if(i2 == 1)
    b1 = dsDNA(i1).end5_1;
elseif(i2 == 2)
    b1 = dsDNA(i1).end3_1;
elseif(i2 == 3)
    b1 = dsDNA(i1).end5_2;
elseif(i2 == 4)
    b1 = dsDNA(i1).end3_2;
else
    error('Exception.');
end

i1 = tour_2(1);
i2 = tour_2(2);
if(i2 == 1)
    b2 = dsDNA(i1).end5_1;
elseif(i2 == 2)
    b2 = dsDNA(i1).end3_1;
elseif(i2 == 3)
    b2 = dsDNA(i1).end5_2;
elseif(i2 == 4)
    b2 = dsDNA(i1).end3_2;
else
    error('Exception.');
end

bAcross1 = dnaTop(b1).across;
bAcross2 = dnaTop(b2).across;

if(~isempty(find([b1,b2,bAcross1,bAcross2] < 0, 1)))
    error('Single-stranded DNA');
end

vecBp1 = dnaTop(bAcross1).xyz - dnaTop(b1).xyz;
vecBp2 = dnaTop(bAcross2).xyz - dnaTop(b2).xyz;
vecB = cross(vecBp1,vecBp2);
if(norm(vecB) < 0.05)
    vecB = zeros(size(vecB));
end

coorC1 = (dnaTop(b1).xyz + dnaTop(bAcross1).xyz)/2;
coorC2 = (dnaTop(b2).xyz + dnaTop(bAcross2).xyz)/2;
vecC = coorC2 - coorC1;

vStack = vecB'*vecC / norm(vecB) / norm(vecC);
vStack = abs(vStack);
if(isnan(vStack))
    vStack = vecBp1'*vecBp2 / norm(vecBp1) / norm(vecBp2);
    assert(abs(vStack) > 0.5);
    vStack = 0.5 - vStack/2;     % convert the range from [-1, 1] to [0, 1]
end

end

function [xyz_prefer, xyz_alt] = nt_coord(dNode, triad)
% Constant parameters
diameterHX = 20;                  % [Angstrom] diameter of a helix
angMinor = 120;                   % [degree] angle of the minor groove
% Initialization
xyz_prefer = zeros(size(dNode));
xyz_alt = zeros(size(dNode));
a_scaf = deg2rad(180 - angMinor / 2);
a_stap = deg2rad(180 + angMinor / 2);
R_HX = diameterHX / 2;           % [Angstrom] radius of a helix
for i = 1 : size(dNode,1)
    v_prefer = triad(:,:,i) * [cos(a_scaf), sin(a_scaf), 0]';
    v_alt = triad(:,:,i) * [cos(a_stap), sin(a_stap), 0]';
    xyz_prefer(i,:) = dNode(i,:) + v_prefer' * R_HX;
    xyz_alt(i,:) = dNode(i,:) + v_alt' * R_HX;
end
end
