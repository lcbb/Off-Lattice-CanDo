function strand_NMA = transformMat_modify(strand,base2node,R_NMA,D_NMA)

strand_NMA = strand;
for i = 1:numel(strand_NMA)
    for j = 1:numel(strand_NMA(i).tour)
        iNode = base2node(strand(i).tour(j));
        if(iNode==0)
            continue;
        end
        strand_NMA(i).R{j} = R_NMA(:,:,iNode) * strand(i).R{j};
        strand_NMA(i).d{j} = strand(i).d{j} + D_NMA(:,iNode);
    end
end

end