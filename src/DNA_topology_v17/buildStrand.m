function [dnaTop,strand] = buildStrand(dnaTop)

nBase = numel(dnaTop);      % number of bases in the DNA nanostructure
iStrand = 0;
isVisited = false(nBase,1);

% Each loop creates a new strand
while(~isempty(find(~isVisited,1)))
    iStrand = iStrand + 1;
    currBase = dnaTop(find(~isVisited,1));
    initBaseID = currBase.id;
    
    % Find the first base in the current strand
    while(currBase.up>=0 && currBase.up~=initBaseID)
        currBase = dnaTop(currBase.up);
        if(isVisited(currBase.id))
            error('Reached a visited base.');
        end
    end
    if(currBase.up<0)   % currBase is at the 5'-end of the strand
        strand(iStrand).isCircular = false;
    elseif(currBase.up==initBaseID) % currBase goes back to the starting point
        strand(iStrand).isCircular = true;
        currBase = dnaTop(initBaseID);
    else
        error('Exception.');
    end
        
    % Walk through the current strand
    iResidue = 1;
    strand(iStrand).tour = currBase.id;
    dnaTop(currBase.id).strand = iStrand;
    dnaTop(currBase.id).residue = iResidue;
    isVisited(currBase.id) = true;
    % Each loop adds a new base
    while((~strand(iStrand).isCircular && currBase.down>=0) || ...
          (strand(iStrand).isCircular && currBase.down~=initBaseID))
        currBase = dnaTop(currBase.down);
        if(isVisited(currBase.id))
            error('Reached a visited base.');
        end
        iResidue = iResidue + 1;
        strand(iStrand).tour = cat(1, strand(iStrand).tour, currBase.id);
        dnaTop(currBase.id).strand = iStrand;
        dnaTop(currBase.id).residue = iResidue;
        isVisited(currBase.id) = true;
    end
end

end