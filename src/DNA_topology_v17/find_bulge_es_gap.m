function [nick_real, bulge, conn_es, conn_gap] = find_bulge_es_gap(nick, unpaired_conn, ssDNA, N_ss)

nick_real = [];
bulge = [];
conn_es = [];
conn_gap = [];

isBulge = false(size(unpaired_conn));

% Categorize nicks into (1) real nicks and (2) bulges
for i = 1 : numel(nick)
    tag = false;
    for j = 1 : numel(unpaired_conn)
        if(max(abs(nick(i).tour(1,:) - unpaired_conn(j).end3)) == 0 && ...
           max(abs(nick(i).tour(4,:) - unpaired_conn(j).end5)) == 0);
            bulge_curr = nick(i);
            bulge_curr.ssDNA_ID = unpaired_conn(j).ssDNA_ID;
            bulge_curr.len = numel(ssDNA(bulge_curr.ssDNA_ID).tour);
            bulge = cat(2, bulge, bulge_curr);
            
            tag = true;
            isBulge(j) = true;
            break;
        end
    end
    
    if(~tag)
        nick_real = cat(2, nick_real, nick(i));
    end
end

% Categorize the ssDNA connecting two dsDNA (excluding bulges) into
% (1) long ones (elastic spring) and (2) short ones (gap)
for i = 1 : numel(unpaired_conn)
    if(isBulge(i))
        continue;
    end
    
    conn_curr.end5 = unpaired_conn(i).end5;
    conn_curr.end3 = unpaired_conn(i).end3;
    conn_curr.len = numel(ssDNA(unpaired_conn(i).ssDNA_ID).tour);
    if(conn_curr.len >= N_ss)
        conn_es = cat(2, conn_es, conn_curr);
    else
        conn_gap = cat(2, conn_gap, conn_curr);
    end
end

end
