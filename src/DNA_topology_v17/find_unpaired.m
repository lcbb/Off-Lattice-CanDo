function [unpaired_conn, unpaired_free] = find_unpaired(conn_dsDNA_ssDNA)

unpaired_conn = [];
unpaired_free = [];

% Extract connections involving ssDNA
assert(isempty(conn_dsDNA_ssDNA(conn_dsDNA_ssDNA(:,1)==0 & conn_dsDNA_ssDNA(:,5)==0, :)));
conn = conn_dsDNA_ssDNA(conn_dsDNA_ssDNA(:,1)==0 | conn_dsDNA_ssDNA(:,5)==0, :);
isVisited = false(size(conn,1),1);

i_conn = 0;
i_free = 0;
list_ssDNA_ID = zeros(size(isVisited));
for i = 1:numel(isVisited)
    if(conn(i,1) == 0)
        list_ssDNA_ID(i) = conn(i,2);
    elseif(conn(i,1) == 1)
        assert(conn(i,5) == 0);
        list_ssDNA_ID(i) = conn(i,6);
    else
        error('Exception.');
    end
end

while(~isempty(find(~isVisited,1)))
    i1_curr = find(~isVisited,1);
    isVisited(i1_curr) = true;
    i2_curr = find(list_ssDNA_ID == list_ssDNA_ID(i1_curr) & ~isVisited);
    
    if(~isempty(i2_curr))    % the current ssDNA has no free end
        assert(numel(i2_curr) == 1);
        isVisited(i2_curr) = true;
        i_conn = i_conn + 1;
        if(conn(i1_curr,5) == 0)
            assert(conn(i2_curr,1) == 0);
            unpaired_conn(i_conn).end5 = conn(i1_curr, 2:3);    % [dsDNA ID, end ID]
            unpaired_conn(i_conn).end3 = conn(i2_curr, 6:7);    % [dsDNA ID, end ID]
        elseif(conn(i1_curr,1) == 0)
            assert(conn(i2_curr,5) == 0)
            unpaired_conn(i_conn).end5 = conn(i2_curr, 2:3);    % [dsDNA ID, end ID]
            unpaired_conn(i_conn).end3 = conn(i1_curr, 6:7);    % [dsDNA ID, end ID]
        else
            error('Exception.');
        end
        unpaired_conn(i_conn).ssDNA_ID = list_ssDNA_ID(i1_curr);
        
    else                     % the current ssDNA has a free end
        i_free = i_free + 1;
        if(conn(i1_curr,5) == 0)
            unpaired_free(i_free).end5 = conn(i1_curr, 2:3);    % [dsDNA ID, end ID]
            unpaired_free(i_free).end3 = [-1, -1];              % [dsDNA ID, end ID]
        elseif(conn(i1_curr,1) == 0)
            unpaired_free(i_free).end5 = [-1, -1];              % [dsDNA ID, end ID]
            unpaired_free(i_free).end3 = conn(i1_curr, 6:7);    % [dsDNA ID, end ID]
        else
            error('Exception.');
        end
        unpaired_free(i_free).ssDNA_ID = list_ssDNA_ID(i1_curr);
    end
end

end