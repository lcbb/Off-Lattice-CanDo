function [junction, nick, conn_cluster] = find_junction_nick(conn_dsDNA_ssDNA)

conn_cluster = [];      % to be classified as either junction or nick
junction = [];
nick = [];

% Extract dsDNA-dsDNA connections
conn = conn_dsDNA_ssDNA(conn_dsDNA_ssDNA(:,1)==1 & conn_dsDNA_ssDNA(:,5)==1, :);
isVisited = false(size(conn,1),1);

i_cluster = 0;
while(~isempty(find(~isVisited,1)))
    i_cluster = i_cluster+1;
    conn_cluster(i_cluster).isCircular = false;
    i_conn = find(~isVisited,1);
    i_initConn = i_conn;
    
    % Find the starting point
    i_prevConn = find_prevConn(conn, i_conn);
    while(~isempty(i_prevConn))
        i_conn = i_prevConn;
        i_prevConn = find_prevConn(conn, i_conn);
        if(i_conn==i_initConn)
            conn_cluster(i_cluster).isCircular = true;
            break;
        end
    end
    i_initConn = i_conn;
    
    % Create the current cluster
    conn_cluster(i_cluster).tour = i_conn;
    isVisited(i_conn) = true;
    i_nextConn = find_nextConn(conn, i_conn);
    while(~isempty(i_nextConn))
        i_conn = i_nextConn;
        i_nextConn = find_nextConn(conn, i_conn);
        if(i_conn==i_initConn)
            break;
        end
        
        assert(~isVisited(i_conn));
        isVisited(i_conn) = true;
        
        conn_cluster(i_cluster).tour = cat(1,conn_cluster(i_cluster).tour,i_conn);
    end
end

% Create junctions & nicks
i_junction = 0;
i_nick = 0;
for i = 1:numel(conn_cluster)
    if(numel(conn_cluster(i).tour) == 1)
        % Is a nick
        i_nick = i_nick+1;
        assert(~conn_cluster(i).isCircular);
        nick(i_nick).tour = -ones(4,2);
        nick(i_nick).tour(2,:) = conn(conn_cluster(i).tour, 2:3);
        nick(i_nick).tour(3,:) = conn(conn_cluster(i).tour, 6:7);
        
        nick(i_nick).tour(1,:) = nick(i_nick).tour(2,:);
        nick(i_nick).tour(1,2) = find_comp_end_ID(nick(i_nick).tour(1,2));
        nick(i_nick).tour(4,:) = nick(i_nick).tour(3,:);
        nick(i_nick).tour(4,2) = find_comp_end_ID(nick(i_nick).tour(4,2));
        
    else
        % Is a N-way junction (N>=3)
        i_junction = i_junction+1;
        junction(i_junction).isCircular = conn_cluster(i).isCircular;
        if(junction(i_junction).isCircular)
            junction(i_junction).nWay = numel(conn_cluster(i).tour);
        else
            junction(i_junction).nWay = numel(conn_cluster(i).tour) + 1;
        end
        junction(i_junction).tour = -ones((numel(conn_cluster(i).tour)+1)*2, 2);
        
        for j = 1 : numel(conn_cluster(i).tour)
            junction(i_junction).tour(j*2,:) = conn(conn_cluster(i).tour(j), 2:3);
            junction(i_junction).tour(j*2+1,:) = conn(conn_cluster(i).tour(j), 6:7);
        end
        
        junction(i_junction).tour(1,:) = junction(i_junction).tour(2,:);
        junction(i_junction).tour(1,2) = find_comp_end_ID(junction(i_junction).tour(1,2));
        junction(i_junction).tour(end,:) = junction(i_junction).tour(end-1,:);
        junction(i_junction).tour(end,2) = find_comp_end_ID(junction(i_junction).tour(end,2));
        
        if(junction(i_junction).isCircular)
            junction(i_junction).tour(end-1:end, :) = [];
        end
    end
end

end


% Tool function: find the previous dsDNA-dsDNA connection
function [i_prevConn] = find_prevConn(conn, i_conn)

i_dsDNA = conn(i_conn,2);
i_end = conn(i_conn,3);

if(i_end == 2)
    i_prevend = 3;
elseif(i_end == 4)
    i_prevend = 1;
else
    error('Exception.');
end

i_prevConn = find(conn(:,6)==i_dsDNA & conn(:,7)==i_prevend);
assert(numel(i_prevConn)==0 || numel(i_prevConn)==1);

end


% Tool function: find the next dsDNA-dsDNA connection
function [i_nextConn] = find_nextConn(conn, i_conn)

i_dsDNA = conn(i_conn,6);
i_end = conn(i_conn,7);

if(i_end == 1)
    i_nextend = 4;
elseif(i_end == 3)
    i_nextend = 2;
else
    error('Exception.');
end

i_nextConn = find(conn(:,2)==i_dsDNA & conn(:,3)==i_nextend);
assert(numel(i_nextConn)==0 || numel(i_nextConn)==1);

end


% Tool function: find the complementary end ID in dsDNA
function [end_ID_comp] = find_comp_end_ID(end_ID)

if(end_ID == 1)
    end_ID_comp = 4;
elseif(end_ID == 2)
    end_ID_comp = 3;
elseif(end_ID == 3)
    end_ID_comp = 2;
elseif(end_ID == 4)
    end_ID_comp = 1;
else
    error('Exception.');
end

end