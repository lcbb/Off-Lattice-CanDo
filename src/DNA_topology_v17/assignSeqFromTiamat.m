function dnaTop = assignSeqFromTiamat(dnaTop, strand, seqPATH)

nBase = numel(dnaTop);
nStrand = numel(strand);

% Initialize the sequence information for each base
for i = 1:nBase
    dnaTop(i).seq = '';
end

% Read the sequence file and assign the sequence to the strands
seq = readSeqFromTiamat(seqPATH);
numseq = numel(seq);
assert(numel(seq) == nStrand);
for i = 1 : nStrand
    for j = 1:numel(strand(i).tour)
        % For the current base
        baseID = strand(i).tour(j);
        dnaTop(baseID).seq = seq{i}(j);

        % For the complementary base
        baseID_comp = dnaTop(baseID).across;
        seq_comp = wspair(seq{i}(j));
        if(baseID_comp >= 0)
            assert(strcmp(dnaTop(baseID_comp).seq, '') || ...
                   strcmp(dnaTop(baseID_comp).seq, seq_comp));
            dnaTop(baseID_comp).seq = seq_comp;
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


% Read all sequences from a Tiamat-exported file
function seq = readSeqFromTiamat(filename)

fid = fopen(filename);

ss = fgetl(fid);
lid = 0;

while(ischar(ss))
    lid = lid + 1;

    C = strsplit(ss, ': ');
    if(strcmp(C{1}(end), 'c'))
        % circular strand
        C{1} = C{1}(1 : end-1);
    end
    assert(numel(C) == 2 && str2double(C{1}) == lid);
    seq_curr = C{2};
    for i = 1 : numel(seq_curr)
        assert(strcmp(seq_curr(i),'A') || ...
               strcmp(seq_curr(i),'G') || ...
               strcmp(seq_curr(i),'C') || ...
               strcmp(seq_curr(i),'T'));
    end
    seq{lid} = seq_curr;

    ss = fgetl(fid);
end

fclose(fid);

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
