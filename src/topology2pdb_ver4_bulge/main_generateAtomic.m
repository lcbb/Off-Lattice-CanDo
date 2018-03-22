function main_generateAtomic(designPATH,postprocDIR,prefixFN,sysParam)

% [~,jobName] = fileparts(designPATH);
% jobName = jobName(1:end-4);     % Remove the suffix '_HJE'
jobName = prefixFN;
assert(~isempty(jobName));
pdbName = fullfile(postprocDIR, sprintf('%s.pdb',jobName));
% CanDoPath = fullfile(postprocDIR, prefixFN);
CanDoPath = postprocDIR;
NMAPath = fullfile(postprocDIR, sprintf('%s_NMA',jobName));


%% Create the topology of the DNA origami
load(designPATH,'dnaTop','strand','HJE');


%% Assign the sequences for the scaffold and staples
fprintf('Assign the sequences for the scaffold and staples ...');
strand = assignSeq(dnaTop, strand);
fprintf('Done.\n');
toc;


%% Translation and rotation of the reference base pair structure
fprintf('Translation and rotation of the reference base pair structure ...');
[strand,base2node] = transformMat(dnaTop,strand,HJE,CanDoPath,jobName);
for i = 1:numel(strand)
    for j = 1:numel(strand(i).tour)
        assert((isempty(strand(i).R{j})&&isempty(strand(i).d{j})) || (~isempty(strand(i).R{j})&&~isempty(strand(i).d{j})));
        if(~isempty(strand(i).R{j}))
            % This matrix convert from 3DNA convention to CanDo convention, see getTriad.m
            strand(i).R{j} = strand(i).R{j} * [0,0,1; 1,0,0; 0,1,0];
        end
    end
    [strand(i).R, strand(i).d, strand(i).isMain] = generateBulgeDOF(strand(i).R, strand(i).d, strand(i).isMain, strand(i).isCircular);
end
fprintf('Done.\n');
toc;


%% Write PDB file
fprintf('Generate PDB Data ...');
pdbFinal = pdbGenerate(strand);
% pdbFinal = pdbModify(pdbFinal);
pdbFinal_singleModel = pdbModify(pdbFinal);
% pdbFinal = pdbModify_chainA(pdbFinal);
pdbFinal = pdbModify_multimodel(pdbFinal);
fprintf('Done.\n');
toc;

fprintf('Write PDB file ...');
mypdbwrite_v2(fullfile(postprocDIR, sprintf('%s_singlemodel.pdb',jobName)),pdbFinal_singleModel);
mypdbwrite_v2(pdbName,pdbFinal);  % For MMB, the same function in version 2 is for visualization
fprintf('Done.\n');
toc;


%% NMA trajectory
if(sysParam.atomicMovie)
    if(exist(NMAPath,'dir') ~= 7)
        mkdir(NMAPath);
    end
    [R_NMA, D_NMA] = NMA_deformation(CanDoPath, jobName);
%     [R_NMA, D_NMA] = NMA_deformation(CanDoPath, jobName, sysParam.T_NMA, sysParam.nFrame, sysParam.nPeriod);
    
    for loop = 1:sysParam.nFrame
        fprintf('---------- Frame %d of %d ----------\n', loop, sysParam.nFrame);
        
        %% Translation and rotation of the reference base pair structure
        fprintf('Translation and rotation of the reference base pair structure ...');
        strand_NMA = transformMat_modify(strand,base2node,R_NMA(:,:,:,loop), D_NMA(:,:,loop));
        fprintf('Done.\n');
        toc;
        
        %% Write PDB file
        fprintf('Generate PDB Data ...');
        pdbFinal = pdbGenerate(strand_NMA);
        pdbFinal = pdbModify_multimodel(pdbFinal);
        fprintf('Done.\n');
        toc;
        
        fprintf('Write PDB file ...');
        pdbName = fullfile(NMAPath, sprintf('%s_%d.pdb',jobName,loop));
        mypdbwrite_v2(pdbName,pdbFinal);  % For MMB, the same function in version 2 is for visualization
        fprintf('Done.\n');
        toc;
    end
end


%% Visualization
fprintf('Visualizing the generated atomic structures ...');
visualDIR = fullfile(postprocDIR, sprintf('%s_visual',jobName));
atomicVisualization(postprocDIR, jobName, visualDIR, strand, sysParam);
fprintf('Done.\n');
toc;

end