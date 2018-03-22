function [dsDNA, ssDNA] = divideTopology(dnaTop)

nBase = numel(dnaTop);      % number of bases in the DNA nanostructure
i_dsDNA = 0;
i_ssDNA = 0;
isVisited = false(nBase,1);

dsDNA = struct('isCircular', {}, 'tour_1', {}, 'tour_2', {}, ...
               'end5_1', {}, 'end3_1', {}, 'end5_2', {}, 'end3_2', {}, ...
               'isConn_end5_1_end3_2', {}, 'isConn_end5_2_end3_1', {});
ssDNA = struct('isCircular', {}, 'tour', {}, 'end5', {}, 'end3', {});


% Each loop creates a new (1) dsDNA *OR* (2) ssDNA
while(~isempty(find(~isVisited,1)))
    currBase = dnaTop(find(~isVisited,1));
    
    if(currBase.across > 0)     % in dsDNA
        i_dsDNA = i_dsDNA + 1;
        [dsDNA(i_dsDNA), isVisited] = build_dsDNA(dnaTop, currBase, isVisited);
    elseif(currBase.across < 0) % in ssDNA
        i_ssDNA = i_ssDNA + 1;
        [ssDNA(i_ssDNA), isVisited] = build_ssDNA(dnaTop, currBase, isVisited);
    else
        error('Exception.');
    end
end

end


function [dsDNA, isVisited] = build_dsDNA(dnaTop, currBase, isVisited)

initBaseID = currBase.id;
initAcrossBaseID = currBase.across;
dsDNA.isCircular = false;

% Find the 5'-end of tour_1
i_1 = initBaseID;
i_2 = initAcrossBaseID;
assert(dnaTop(i_1).across == i_2 && dnaTop(i_2).across == i_1);
while(dnaTop(i_1).up>0 && dnaTop(i_2).down>0 && ...
      dnaTop(i_1).up~=i_2 && ...
      dnaTop(dnaTop(i_1).up).across == dnaTop(i_2).down && ...
      dnaTop(dnaTop(i_2).down).across == dnaTop(i_1).up)
    i_1 = dnaTop(i_1).up;
    i_2 = dnaTop(i_2).down;
    if(i_1 == initBaseID)
        assert(i_2 == initAcrossBaseID);
        dsDNA.isCircular = true;
        break;
    end
end
initBaseID = i_1;
initAcrossBaseID = i_2;

% Create the current dsDNA
dsDNA.tour_1 = i_1;
dsDNA.tour_2 = i_2;
assert(dnaTop(i_1).across == i_2 && dnaTop(i_2).across == i_1);
assert(~isVisited(i_1) && ~isVisited(i_2));
isVisited(i_1) = true;
isVisited(i_2) = true;
while(dnaTop(i_1).down>0 && dnaTop(i_2).up>0 && ...
      dnaTop(i_1).down~=i_2 && ...
      dnaTop(dnaTop(i_1).down).across == dnaTop(i_2).up && ...
      dnaTop(dnaTop(i_2).up).across == dnaTop(i_1).down)
    i_1 = dnaTop(i_1).down;
    i_2 = dnaTop(i_2).up;
    if(i_1 == initBaseID)
        assert(i_2 == initAcrossBaseID && dsDNA.isCircular);
        break;
    end
    
    assert(~isVisited(i_1) && ~isVisited(i_2));
    isVisited(i_1) = true;
    isVisited(i_2) = true;
    
    dsDNA.tour_1 = cat(1,dsDNA.tour_1,i_1);
    dsDNA.tour_2 = cat(1,i_2,dsDNA.tour_2);
end

% Further annotation
dsDNA.end5_1 = dsDNA.tour_1(1);
dsDNA.end3_1 = dsDNA.tour_1(end);
dsDNA.end5_2 = dsDNA.tour_2(1);
dsDNA.end3_2 = dsDNA.tour_2(end);

dsDNA.isConn_end5_1_end3_2 = false;
if(dnaTop(dsDNA.end5_1).up == dsDNA.end3_2)
    assert(dsDNA.end5_1 == dnaTop(dsDNA.end3_2).down);
    dsDNA.isConn_end5_1_end3_2 = true;
    dsDNA.end5_1 = -1;
    dsDNA.end3_2 = -1;
end

dsDNA.isConn_end5_2_end3_1 = false;
if(dnaTop(dsDNA.end5_2).up == dsDNA.end3_1)
    assert(dsDNA.end5_2 == dnaTop(dsDNA.end3_1).down);
    dsDNA.isConn_end5_2_end3_1 = true;
    dsDNA.end5_2 = -1;
    dsDNA.end3_1 = -1;
end

if(dsDNA.isCircular)
    assert(~dsDNA.isConn_end5_1_end3_2 && ~dsDNA.isConn_end5_2_end3_1);
    dsDNA.end5_1 = -1;
    dsDNA.end3_2 = -1;
    dsDNA.end5_2 = -1;
    dsDNA.end3_1 = -1;
end

end


function [ssDNA, isVisited] = build_ssDNA(dnaTop, currBase, isVisited)

initBaseID = currBase.id;
assert(dnaTop(initBaseID).across < 0);
ssDNA.isCircular = false;

% Find the 5'-end of tour
i_1 = initBaseID;
while(dnaTop(i_1).up>0 && dnaTop(dnaTop(i_1).up).across < 0)
    i_1 = dnaTop(i_1).up;
    if(i_1 == initBaseID)
        ssDNA.isCircular = true;
        break;
    end
end
initBaseID = i_1;

% Create the current ssDNA
ssDNA.tour = i_1;
assert(~isVisited(i_1));
isVisited(i_1) = true;
while(dnaTop(i_1).down>0 && dnaTop(dnaTop(i_1).down).across<0)
    i_1 = dnaTop(i_1).down;
    if(i_1 == initBaseID)
        assert(ssDNA.isCircular);
        break;
    end
    
    assert(~isVisited(i_1));
    isVisited(i_1) = true;
    
    ssDNA.tour = cat(1,ssDNA.tour,i_1);
end

% Further annotation
ssDNA.end5 = ssDNA.tour(1);
ssDNA.end3 = ssDNA.tour(end);

if(ssDNA.isCircular)
    ssDNA.end5 = -1;
    ssDNA.end3 = -1;
end

end