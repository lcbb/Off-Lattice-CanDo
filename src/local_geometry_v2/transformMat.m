function [strand,base2node] = transformMat(dnaTop,strand,HJE,FEDIR,prefixFN)

%% Load data
load(fullfile(FEDIR,strcat(prefixFN,'_DOF.mat')), 'dNode','triad');

%% Map table from node IDs to base IDs
node2base = zeros(0,2);
for i = 1:numel(HJE)
    node2base = cat(1, node2base, [HJE(i).idSeq_main, HJE(i).idSeq_comp]);
end

%% Change the sequence of the node from ADINA convention to Tiamat topology base ID convention
nStrand = numel(strand);
nHJE = numel(HJE);
nBase = numel(dnaTop);
base2node = zeros(nBase,1);

iNode = 0;
for i = 1:nHJE
    for j = 1:numel(HJE(i).idSeq_main)
        iNode = iNode + 1;
        base2node(HJE(i).idSeq_main(j)) = iNode;
    end
end
iNode = 0;
for i = 1:nHJE
    for j = 1:numel(HJE(i).idSeq_comp)
        iNode = iNode + 1;
        base2node(HJE(i).idSeq_comp(j)) = iNode;
    end
end
for i = 1:numel(base2node)
    if(base2node(i)==0 && dnaTop(i).across>=0)
        error('Undefined base-pair position and orientation.');
    end
end

%% Get 6 DOFs for each base
for i = 1:nStrand
    for j = 1:numel(strand(i).tour)
        strand(i).R{j} = [];
        strand(i).d{j} = [];
        iNode = base2node(strand(i).tour(j));
        if(iNode > 0)
            strand(i).R{j} = triad{iNode};
            strand(i).d{j} = dNode(iNode,:)';
        end
        % Figure out if some bases are in 'main' or in 'comp'
        if(isinf(strand(i).isMain(j)))
            if(node2base(iNode,1) == strand(i).tour(j))
                strand(i).isMain(j) = true;
            elseif(node2base(iNode,2) == strand(i).tour(j))
                strand(i).isMain(j) = false;
            else
                error('Exception.');
            end
        end
    end
end

end