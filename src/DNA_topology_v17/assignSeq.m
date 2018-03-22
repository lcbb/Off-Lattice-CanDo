function dnaTop = assignSeq(dnaTop, strand, seqPATH)

nBase = numel(dnaTop);
nStrand = numel(strand);
assert(nStrand == numel(seqPATH));

% Initialize the sequence information for each base
for i = 1:nBase
    dnaTop(i).seq = '';
end

% Read the sequence file and assign the sequence to the strands
for i = 1:nStrand
    if(~isempty(seqPATH{i}))
        seq = readSeq(seqPATH{i});
        for j = 1:numel(strand(i).tour)
            % For the current base
            baseID = strand(i).tour(j);
            dnaTop(baseID).seq = seq(j);
            
            % For the complementary base
            baseID_comp = dnaTop(baseID).across;
            seq_comp = wspair(seq(j));
            if(baseID_comp >= 0)
                assert(strcmp(dnaTop(baseID_comp).seq, '') || ...
                       strcmp(dnaTop(baseID_comp).seq, seq_comp));
                dnaTop(baseID_comp).seq = seq_comp;
            end
        end
    end
end

% Make sure that the sequence information are defined for the dsDNA
for i = 1:nBase
%    assert(~(~dnaTop(i).isSingleStrand & strcmp(dnaTop(i).seq,'')));
    if(dnaTop(i).across >= 0 && strcmp(dnaTop(i).seq,''))
        error('The first undefined dsDNA base is base %d (strand %d, residue %d)', i, dnaTop(i).strand, dnaTop(i).residue);
    end
end

end


% Read a single sequence file
function seq = readSeq(filename)

fid = fopen(filename);
seq = fscanf(fid,'%c'); 
fclose(fid);

for i = 1:numel(seq)
    assert(strcmp(seq(i),'A') || ...
           strcmp(seq(i),'G') || ...
           strcmp(seq(i),'C') || ...
           strcmp(seq(i),'T'));
end

end


% Get the complementary base
function y = wspair(x)

if('A'==x || 'a'==x)
    y = 'T';
elseif('G'==x || 'g'==x)
    y = 'C';
elseif('C'==x || 'c'==x)
    y = 'G';
elseif('T'==x || 't'==x)
    y = 'A';
else
    error('Illegal base.');
end

end