function conn_dsDNA_ssDNA = findConn(dnaTop, dsDNA, ssDNA)

% Format of conn_dsDNA_ssDNA:
% Column 1: dsDNA (1) or ssDNA (0) in the 5'-direction
% Column 2: ID of dsDNA or ssDNA in the 5'-direction
% Column 3: ID of ends (1-4 for dsDNA, 1-2 for ssDNA) in the 5'-direction
% Column 4: base ID
% Column 5: dsDNA (1) or ssDNA (0) in the 3'-direction
% Column 6: ID of dsDNA or ssDNA in the 3'-direction
% Column 7: ID of ends (1-4 for dsDNA, 1-2 for ssDNA) in the 3'-direction
% Column 8: base ID
conn_dsDNA_ssDNA = zeros(0,8);

n_dsDNA = numel(dsDNA);
n_ssDNA = numel(ssDNA);

% Format of end5_list & end3_list:
% Column 1: dsDNA (1) or ssDNA (0)
% Column 2: ID of dsDNA or ssDNA
% Column 3: ID of ends (1-4 for dsDNA, 1-2 for ssDNA)
% Column 4: base ID
end5_list = zeros(0,4);
end3_list = zeros(0,4);

for i = 1:n_dsDNA
    if(dsDNA(i).end5_1 > 0)
        assert(dsDNA(i).end3_2 > 0);
        end5_list = cat(1, end5_list, [1,i,1,dsDNA(i).end5_1]);
    end
    if(dsDNA(i).end3_1 > 0)
        assert(dsDNA(i).end5_2 > 0);
        end3_list = cat(1, end3_list, [1,i,2,dsDNA(i).end3_1]);
    end
    if(dsDNA(i).end5_2 > 0)
        assert(dsDNA(i).end3_1 > 0);
        end5_list = cat(1, end5_list, [1,i,3,dsDNA(i).end5_2]);
    end
    if(dsDNA(i).end3_2 > 0)
        assert(dsDNA(i).end5_1 > 0);
        end3_list = cat(1, end3_list, [1,i,4,dsDNA(i).end3_2]);
    end
end

for i = 1:n_ssDNA
    if(ssDNA(i).end5 > 0)
        assert(ssDNA(i).end3 > 0);
        end5_list = cat(1, end5_list, [0,i,1,ssDNA(i).end5]);
    end
    if(ssDNA(i).end3 > 0)
        assert(ssDNA(i).end5 > 0);
        end3_list = cat(1, end3_list, [0,i,2,ssDNA(i).end3]);
    end
end

assert(size(end5_list,1) == size(end3_list,1));
while(~isempty(end3_list))
    conn_curr = zeros(1,8);
    conn_curr(1:4) = end3_list(1,:);
    currBase = end3_list(1,4);    % the base ID in the 5'-direction
    end3_list(1,:) = [];

    nextBase = dnaTop(currBase).down;  % the base ID in the 3'-direction, if existing
    if(nextBase > 0)
        tmp = find(end5_list(:,4) == nextBase);
        assert(numel(tmp)==1);
        conn_curr(5:8) = end5_list(tmp,:);
        end5_list(tmp,:) = [];
        conn_dsDNA_ssDNA = cat(1,conn_dsDNA_ssDNA,conn_curr);   % add the connection to the results
    end
end

end