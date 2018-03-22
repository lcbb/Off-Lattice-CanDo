function strand = assignSeq(dnaTop, strand)

nStrand = numel(strand);

for i = 1:nStrand
    strand(i).seq = [];
    for j = 1:numel(strand(i).tour)
        currentSeq = dnaTop(strand(i).tour(j)).seq;
        strand(i).seq = cat(2, strand(i).seq, currentSeq);
    end
end

end