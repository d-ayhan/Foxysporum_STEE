% generate genome file for circos

f = fastaread('Fol4287_GCA_003315725_genomic.fna');
a = importdata('Fol4287_GCA_003315725_genomic.fna.fai');

sections = [1,2,3,4,5,6,7,8,9,10,13,19,151,197,239]';
sections_headers = {'chr1', 'chr2', 'chr4', 'chr5', 'chr7', 'chr8', 'chr9', ...
    'chr10', 'chr11', 'chr12', 'chr13', 'chr2b', 'chr3', 'chr14', 'chr15'};

[coverage, order, ~] = xlsread('clustered_scaffolds_with_coverages.xlsx');
order = order(:,1);

[Lia,Loc] = ismember(order,a.textdata);

reordered_f = f(Loc(Lia));
reordered_a.textdata = a.textdata(Loc(Lia));
reordered_a.data = a.data(Loc(Lia),:);

for i = 1:length(sections)
    i1 = sections(i);
    if i == 1
         i0 = 0;
    else
        i0 = sections(i-1);
    end
    new_f(i).Header = sections_headers{i};
    new_f(i).Sequence = [reordered_f(i0+1:i1).Sequence];
    
    new_f_indexing(i0+1,1) = 0;
    new_f_indexing(i0+1,2) = i;
    temp = reordered_a.data(i0+1,1);
    if i1-i0 > 1
        for j = i0+2:i1
            new_f_indexing(j,1) = temp;
            temp = temp + reordered_a.data(j,1); % (j,1)
            new_f_indexing(j,2) = i;
        end
    end
end
lf = [];
for i = 1:length(new_f)
    lf(i,1) = length(new_f(i).Sequence);
end
        
save('format_genome_file.mat','reordered_a','new_f','new_f_indexing')
         