load('format_genome_file.mat')
files = {'H3K4me2_1.genomecov.txt','H3K4me2_2.genomecov.txt','H3K27me3_1.genomecov.txt','H3K27me3_2.genomecov.txt'};
window_size = 10000;
samples = {'H3K4me2','H3K27me3'};

readcounts = importdata('histone_methylation_total_read_counts.txt');

for i = 1: length(new_f)
    F{i} = zeros(length(new_f(i).Sequence),4);
end

for k=1:4
    T= readtable(files{k});
    for i = 1:length(reordered_a.textdata)
        Lia = strcmp(reordered_a.textdata{i},T.Var1);
        ind = new_f_indexing(i,:);
        F{ind(1,2)}(ind(1,1)+1:ind(1,1)+sum(Lia),k) = T.Var3(Lia);
    end
end
clearvars T


readcounts.data(1,1) = readcounts.data(1,1) + readcounts.data(2,1);
readcounts.data(2,1) = readcounts.data(3,1) + readcounts.data(4,1);
readcounts.data(3:4,:) = [];

for i = 1:length(F)
    F{i}(:,1) = F{i}(:,1) + F{i}(:,2);
    F{i}(:,2) = F{i}(:,3) + F{i}(:,4);
    F{i}(:,3:4) = [];
end
    
    

for z = 1:2
    chr = {};
position=[];
varsum = [];
    for i = 1:length(F)
        l = length(F{i});
        for j = 1:window_size:l
            k = window_size + j;
            if k > l
                k = l;
            end
            chr{end+1,1} = new_f(i).Header;
            position(end+1,1:2) = [j-1,k-1];
            varsum(end+1,1) = median(F{i}(j:k-1,z));
        end
    end
    if z == 1
    T1 = table(chr, position, varsum);
    else
        T2 = table(chr, position, varsum);
    end
    
    %writetable(T, [samples{z} '_for_circos.txt'], 'Delimiter','\t', 'WriteVariableNames',0)
end

scaling_factor =readcounts.data ./ readcounts.data(2);
temp = T2.varsum .* scaling_factor(1);
T2.varsum = temp;

writetable(T1, [samples{1} '_for_circos.txt'], 'Delimiter','\t', 'WriteVariableNames',0)
writetable(T2, [samples{2} '_for_circos.txt'], 'Delimiter','\t', 'WriteVariableNames',0)
