function pdb_multimodel = pdbModify_multimodel(pdbOld)

pdbNew = pdbOld;
nChain = length(pdbNew.Model.Terminal);
chainList = cell(nChain,3);

% Start atomSerNo and end atomSerNo
atomSerNo_start = 1;
for i = 1:nChain
    chainList{i,1} = atomSerNo_start;
    chainList{i,2} = pdbNew.Model.Terminal(i).SerialNo - i;
    chainList{i,3} = assignChainID(i);
    atomSerNo_start = chainList{i,2} + 1;
end

% Assign the new chain IDs
for i = 1:nChain
    AtomSerNo_start = pdbNew.Model.Atom(chainList{i,1}).AtomSerNo;
    for j = chainList{i,1}:chainList{i,2}
        pdbNew.Model.Atom(j).AtomSerNo = pdbNew.Model.Atom(j).AtomSerNo - AtomSerNo_start + 1;
        pdbNew.Model.Atom(j).chainID = ' A';
    end
    pdbNew.Model.Terminal(i).SerialNo = pdbNew.Model.Atom(chainList{i,2}).AtomSerNo + 1;
    pdbNew.Model.Terminal(i).chainID = ' A';
end

% Create multiple models
pdb_multimodel.Model(nChain) = struct('Atom',[], 'Terminal',[]);
for i = 1:nChain
    pdb_multimodel.Model(i).Atom = pdbNew.Model.Atom(chainList{i,1}:chainList{i,2});
    pdb_multimodel.Model(i).Terminal = pdbNew.Model.Terminal(i);
end

end


% Assign chain ID
function ID = assignChainID(index)

% if(index>26*27)
%     error('The number of scaffold/staple strands exceed the limit');
% end

start = 'A'-1;
if(index<=26)
    ID = sprintf(' %s', start+index);
else
    index1 = mod(index-1,26)+1;
    index2 = floor((index-1)/26);
    ID = sprintf('%c%c', start+index2, start+index1);
end

end