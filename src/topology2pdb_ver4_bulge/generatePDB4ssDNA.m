function colorList = generatePDB4ssDNA(refPDB, ssDNA_PDB, strand, strandColor)

nStrand = numel(strand);
pdbStruct = pdb2struct(refPDB);
bPts = [];
colorList = [];

% es is a 1xn cell array. es{i} is a 2xm matrix.
% es{i}(1,:) are the residue IDs of the starting bases.
% es{i}(2,:) are the residue IDs of the ending bases.
es = findEntSpring(strand);
    
for i = 1:nStrand
    currPts = findBezierPoints(es{i},pdbStruct,i);
    bPts = cat(1,bPts,currPts);
    colorList = cat(1,colorList,repmat(strandColor(i,:),numel(currPts),1));
end

% Generate the Bezier curves
fid = fopen(ssDNA_PDB,'w');
for i = 1:numel(bPts)
    p = cubicBezier(linspace(0,1,64)',bPts{i});
    fprintf(fid, 'MODEL     %4d\n',i);
    for j = 1:size(p,1)
        fprintf(fid, 'ATOM  %5d  CA  ALA A%4d    %8.3f%8.3f%8.3f%6.2f%6.2f           C  \n', ...
                j, j, p(j,1), p(j,2), p(j,3), 1, 0);
    end
    fprintf(fid, 'ENDMDL\n');
end
fclose(fid);

end


% Find residue IDs of start/end of ssDNAs
function es = findEntSpring(strand)

es = cell(size(strand));
nStrand = numel(strand);

% Find the Bezier points and colors
for i = 1:nStrand
    es{i} = zeros(2,0);
    nBase = numel(strand(i).tour);
    
    % Whether a base is single-straned
    isSS = false(nBase,1);
    for j = 1:nBase
        isSS(j) = isempty(strand(i).R{j});
    end
    
    % Find the starts and ends of ssDNAs
    iCurrES = 1;
    for j = 1:(nBase-1)
        if(~isSS(j) && isSS(j+1))
            es{i}(1,iCurrES) = j;
        elseif(isSS(j) && ~isSS(j+1))
            es{i}(2,iCurrES) = j+1;
            iCurrES = iCurrES+1;
        end 
    end
    
    % Is there any ssDNA from the end to the beginning of the current
    % strand?
    if(strand(i).isCircular)
        if(~isSS(nBase) && isSS(1))
            es{i}(1,iCurrES) = nBase;
        elseif(isSS(nBase) && ~isSS(1))
            es{i}(2,iCurrES) = 1;
        end
        
        if(~isempty(es{i}) && 0==es{i}(2,end))
            assert(es{i}(1,end)>0 && es{i}(1,1)==0 && size(es{i},2)>1);
            es{i}(1,1) = es{i}(1,end);
            es{i}(:,end) = [];
        end
    end
end

end


% Create control points for Bezier curves
function bPts = findBezierPoints(es,pdbStruct,chainSerNo)

assert(size(es,1)==2);
if(~isempty(es) && es(1,1)==0)
    es(:,1) = [];
end

startList = [];
endList = [];
for i = 1:size(es,2)
    if(es(1,i)>0)
        startList = cat(2,startList,es(1,i));
    end
    if(es(2,i)>0)
        endList = cat(2,endList,es(2,i));
    end
end
assert(numel(startList)==numel(endList) || numel(startList)==numel(endList)+1);

bPts = cell(numel(endList),1);
for i = 1:numel(endList)
    iP0 = find(pdbStruct.chainSerNo==chainSerNo & ...
               pdbStruct.resSeq==startList(i) & ...
               strcmp(pdbStruct.AtomName,'P'));
%     assert(numel(iP0)==1);
    if(numel(iP0)~=1)
        error('Cannot find P0 for ssDNA. chainSerNo=%d, resSeq=%d', chainSerNo, startList(i));
    end
    xyzP0 = pdbStruct.XYZ(iP0,:);
    
    iP1 = find(pdbStruct.chainSerNo==chainSerNo & ...
               pdbStruct.resSeq==startList(i) & ...
               strcmp(pdbStruct.AtomName,'C3'''));
    assert(numel(iP1)==1);
    xyzP1 = pdbStruct.XYZ(iP1,:);

    iP2 = find(pdbStruct.chainSerNo==chainSerNo & ...
               pdbStruct.resSeq==endList(i) & ...
               strcmp(pdbStruct.AtomName,'P'));
    assert(numel(iP2)==1);
    xyzP2 = pdbStruct.XYZ(iP2,:);
    
    iP3 = find(pdbStruct.chainSerNo==chainSerNo & ...
               pdbStruct.resSeq==endList(i) & ...
               strcmp(pdbStruct.AtomName,'C3'''));
    assert(numel(iP3)==1);
    xyzP3 = pdbStruct.XYZ(iP3,:);
    
    assert(norm(xyzP0-xyzP1)>0 && norm(xyzP2-xyzP3)>0);
    
    bPts{i} = zeros(4,3);
    bPts{i}(1,:) = xyzP1;
    bPts{i}(2,:) = xyzP1 + (xyzP1-xyzP0)/norm(xyzP1-xyzP0)*20;
    bPts{i}(3,:) = xyzP2 + (xyzP2-xyzP3)/norm(xyzP2-xyzP3)*20;
    bPts{i}(4,:) = xyzP2;
end

end


% 4-point (cubic) Bezier curve
function b = cubicBezier(t,P)

assert(size(t,2)==1 && size(P,1)==4 && size(P,2)==3);

a0 = (1-t).^3;
a1 = 3 * (1-t).^2 .* t;
a2 = 3 * (1-t) .* t.^2;
a3 = t.^3;

b0 = a0*P(1,:);
b1 = a1*P(2,:);
b2 = a2*P(3,:);
b3 = a3*P(4,:);

b = b0+b1+b2+b3;

end