%LoadBedtoolsGenomecovFile('Fol4287_GCA_003315725_genomic.fna.repeatmasker202006_TEonly.bed','Fol4287_GCA_003315725_genomic.fna.fai',10000,'mat')

load('format_genome_file.mat')

window_size = 10000;

for i = 1: length(new_f)
    F{i} = zeros(length(new_f(i).Sequence),1);
end


T= readtable('Fol4287_GCA_003315725_genomic.fna.repeatmasker202006_TEonly.bed','FileType','text');
T.Var3(T.Var3 > 1) = 1;
for i = 1:length(reordered_a.textdata)
        Lia = strcmp(reordered_a.textdata{i},T.Var1);
        ind = new_f_indexing(i,:);
        F{ind(1,2)}(ind(1,1)+1:ind(1,1)+sum(Lia),1) = T.Var3(Lia);
    end

clearvars T

chr = {};
position=[];
varsum = [];
for z = 1:2
    for i = 1:length(F)
        l = length(F{i});
        for j = 1:window_size:l
            k = window_size + j;
            if k > l
                k = l;
            end
            chr{end+1,1} = new_f(i).Header;
            position(end+1,1:2) = [j-1,k-1];
            varsum(end+1,1) = sum(F{i}(j:k-1))./window_size;
        end
    end
    
    T = table(chr, position, varsum);
    
    writetable(T, 'TEdist_for_circos.txt', 'Delimiter','\t', 'WriteVariableNames',0)
end