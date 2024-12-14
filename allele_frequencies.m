load('format_variants_file.mat')

conditions.short = {'P', 'Y' ,'M'};
conditions.long = {'Plant', 'Rich media', 'Minimal media'};

figname = 'allele_frequencies.pdf';
f = figure('Name',figname','Color','w','Position',[680, 558, 1188, 420]);
for i = 1:3
    ax(i) = subplot(2,3,i);
    ax(i+3) = subplot(2,3,i+3);

    indx = cellfun(@(S) strcmp(S(1),conditions.short{i}), inste.sample);
    histogram(ax(i+3), inste.AF(indx), 'BinWidth',0.1);
    
    indx = cellfun(@(S) strcmp(S(1),conditions.short{i}), snps.sample);
    histogram(ax(i), snps.AF(indx), 'BinWidth',0.1);
    ax(i).Title.String = ['TE insertions, ' conditions.long{i} '-passaged'];
    ax(i+3).Title.String = ['SNPs & INDELs, ' conditions.long{i} '-passaged'];
end

for i = 1:6
    ax(i).XLim = [0 1];
    ax(i).YLim = [0 40];
    ax(i).XLabel.String = 'AF';
    ax(i).YLabel.String = 'Counts';
end

print_pdf(f, figname)