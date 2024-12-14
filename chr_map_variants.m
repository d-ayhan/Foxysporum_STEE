window_size = 100000;
load('format_genome_file.mat')

for i = 1: length(new_f)
    F{i} = zeros(length(new_f(i).Sequence),2);
end

TEINS = readtable('Variants_202005.xlsx','Sheet','TEINS');
[~, idx] = unique(TEINS.mut_id);
TEINS = TEINS(idx, :);
pos = floor(mean([TEINS.start,TEINS.end_pos],2));

for i = 1: height(TEINS)
    Lia = strcmp(reordered_a.textdata, TEINS.CHROM{i});
    if any(Lia)
        ind = new_f_indexing(Lia,:);
        
        F{ind(1,2)}(ind(1,1) + pos(i),1) = F{ind(1,2)}(ind(1,1) + pos(i),1) + 1;
    end
end

SNP = readtable('Variants_202005.xlsx','Sheet','SNP_INDEL');
[~, idx] = unique(SNP.mut_id);
SNP = SNP(idx, :);
pos = SNP.POS;

for i = 1: height(SNP)
    Lia = strcmp(reordered_a.textdata, SNP.CHROM{i});
    if any(Lia)
        ind = new_f_indexing(Lia,:);
        
        F{ind(1,2)}(ind(1,1) + pos(i),2) = F{ind(1,2)}(ind(1,1) + pos(i),2) + 1;
    end
end
%get bins

chr = {};
position=[];
varsum = [];
for z = 1:2
    for i = 1: length(F)
        l = length(F{i});
        for j = 1:window_size:l
            k = window_size + j;
            if k > l
                k = l;
            end
            chr{end+1,1} = new_f(i).Header;
            position(end+1,1:2) = [j-1,k-1];
            varsum(end+1) = sum(F{i}(j:k,z));
        end
    end
    
    T = table(chr, position, varsum(:));
    if z == 1
        filen = 'TEins_for_circos_100k.txt';
    else
        filen = 'SNPs_for circos_100k.txt';
    end
    writetable(T, filen, 'Delimiter','\t', 'WriteVariableNames',0)
end