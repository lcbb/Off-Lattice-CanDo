% This function create a checktable for dsDNA-dsDNA connection
% Format of the checktable:
% Row i, Column 1-3 or 4-6 for nick/bulge:
%        Column 1-3: [dsDNA ID, 'left', type] or [dsDNA ID, 'right', type]
%            at the left end of the i-th dsDNA
%        Column 4-6: [dsDNA ID, 'left', type] or [dsDNA ID, 'right', type]
%            at the right end of the i-th dsDNA
% Row i, Column 1-3 or 4-6 for junction:
%        Column 1-3: [junction ID, arm ID, 'junction']
%            at the left end of the i-th dsDNA
%        Column 4-6: [junction ID, arm ID, 'junction']
%            at the right end of the i-th dsDNA

function checktable = find_checktable(nick_real, bulge, junction, n_dsDNA)

checktable = cell(n_dsDNA, 6);

% For real nicks & bulges
checktable = add_nick_to_checktable(checktable, nick_real, 'nick');
checktable = add_nick_to_checktable(checktable, bulge, 'bulge');
% For junctions
checktable = add_junction_to_checktable(checktable, junction);

end


function checktable = add_nick_to_checktable(checktable, nick, conn_type)

n_nick = numel(nick);
for i = 1 : n_nick
    t = nick(i).tour;
    assert(size(t, 1) == 4 && size(t, 2) == 2);
    id_1 = t(1, 1);
    id_2 = t(3, 1);
    assert(t(1, 1) == t(2, 1) && t(3, 1) == t(4, 1));

    if(t(1, 2) == 1)
        assert(t(2, 2) == 4);
        p_1 = 'left';
        i_1 = 1;
    elseif(t(1, 2) == 3)
        assert(t(2, 2) == 2);
        p_1 = 'right';
        i_1 = 4;
    else
        error('Exception.');
    end
    
    if(t(3, 2) == 1)
        assert(t(4, 2) == 4)
        p_2 = 'left';
        i_2 = 1;
    elseif(t(3, 2) == 3)
        assert(t(4, 2) == 2);
        p_2 = 'right';
        i_2 = 4;
    else
        error('Exception.');
    end
    
    assert(isempty(checktable{id_1, i_1}));
    checktable{id_1, i_1} = id_2;
    checktable{id_1, i_1 + 1} = p_2;
    checktable{id_1, i_1 + 2} = conn_type;
    assert(isempty(checktable{id_2, i_2}));
    checktable{id_2, i_2} = id_1;
    checktable{id_2, i_2 + 1} = p_1;
    checktable{id_2, i_2 + 2} = conn_type;
end

end


function checktable = add_junction_to_checktable(checktable, junction)

n_junction = numel(junction);
for i = 1 : n_junction
    t = junction(i).tour;
    assert(size(t, 1) == 8 && size(t, 2) == 2);
    id = t([1 3 5 7], 1);
    assert(t(1, 1) == t(2, 1) && t(3, 1) == t(4, 1) && t(5, 1) == t(6, 1) && t(7, 1) == t(8, 1));
    
    for j = 1 : numel(id)
        if(t(j * 2 - 1, 2) == 1)
            assert(t(j * 2, 2) == 4);
            ip = 1;
        elseif(t(j * 2 - 1, 2) == 3)
            assert(t(j * 2, 2) == 2);
            ip = 4;
        else
            error('Exception.');
        end
        
        assert(isempty(checktable{id(j), ip}));
        checktable{id(j), ip} = i;
        checktable{id(j), ip + 1} = j;
        checktable{id(j), ip + 2} = 'junction';
    end
end

end