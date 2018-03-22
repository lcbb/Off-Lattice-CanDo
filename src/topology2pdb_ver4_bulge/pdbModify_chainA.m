function pdbNew = pdbModify_chainA(pdbOld)

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
    for j = chainList{i,1}:chainList{i,2}
        pdbNew.Model.Atom(j).chainID = chainList{i,3};
    end
    pdbNew.Model.Terminal(i).chainID = chainList{i,3};
end

end


% Assign chain ID
function ID = assignChainID(index)

if(index>26*27)
    error('The number of scaffold/staple strands exceed the limit');
end

start = 'A'-1;
% if(index<=26)
%     ID = sprintf(' %s', start+index);
% else
%     index1 = mod(index-1,26)+1;
%     index2 = floor((index-1)/26);
%     ID = sprintf('%c%c', start+index2, start+index1);
% end
ID = sprintf(' %s', start+1);

end