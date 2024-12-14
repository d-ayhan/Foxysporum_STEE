% analyze the histone modifications on mutation positions
samplecolors = [54, 214, 52; ... P
    190, 92, 255; ...Y
    60, 72, 255, ...M
    ]./255;

%% get files
load('histone_methylation_covdata1000.mat')
a=importdata('chrs.txt');
load('format_variants_file_for_methylation_analysis.mat')

%clean the data
[~, idx] = unique(inste.mut_id);
inste = inste(idx,:);

[~, idx] = unique(snps.mut_id);
snps = snps(idx,:);

[Lia, Loc] = ismember(snps.CHROM, a.textdata);
snps.comp(Lia) = a.data(Loc(Lia));

[Lia, Loc] = ismember(inste.CHROM, a.textdata);
inste.comp(Lia) = a.data(Loc(Lia));

clear Lia Loc idx

%% find mutation loci in meth data
[~, Loc] = ismember(snps.CHROM, a.textdata);
meth.snp = nan([length(Loc),2]);

for j = 1: length(Loc)
    i = Loc(j);
    k = snps.POS(j);
    meth.snp(j,:) = dataCell{i,1}(k >= dataCell{i,1}(:,1) & k <= dataCell{i,1}(:,2),3:4);
end

[Lia, Loc] = ismember(inste.CHROM,a.textdata);
meth.te = nan([length(Loc),2]);

for j = 1: length(Loc)
    i = Loc(j);
    k = mean(inste.start(j),inste.end_pos(j));
    meth.te(j,:) = dataCell{i,1}(k >= dataCell{i,1}(:,1) & k <= dataCell{i,1}(:,2),3:4);
end
%% divide the genome by core (0), fast-core (1), and ls (2)
dC.core = cell2mat(dataCell(a.data(:,1)==0,1));
dC.core2 = cell2mat(dataCell(a.data(:,1)==1,1));
dC.ls = cell2mat(dataCell(a.data(:,1)==2,1));
dC.all = cell2mat(dataCell(:,1));

%% all methylation
figname1 = 'mutations_and_histone_methylation_background.pdf';
F1 = figure('Color','w', 'Name', figname1, 'Position', [301, 390, 1508, 420]);
subplot(1,3,1)
scatter(dC.core(:,3), dC.core(:,4), 'Marker','.','MarkerEdgeColor',  [0,0,0])
xlabel('H3K4me2')
ylabel('H3K27me3')
title('Core')
xlim([0,100])
ylim([0,100])

subplot(1,3,2)
scatter(dC.core2(:,3), dC.core2(:,4), 'Marker','.','MarkerEdgeColor',  [0,0,0])
xlabel('H3K4me2')
ylabel('H3K27me3')
title('Fast Core')
xlim([0,100])
ylim([0,100])


subplot(1,3,3)
scatter(dC.ls(:,3), dC.ls(:,4), 'Marker','.','MarkerEdgeColor', [0,0,0])
xlabel('H3K4me2')
ylabel('H3K27me3')
title('LS')
xlim([0,100])
ylim([0,100])

print_pdf(F1, figname1)

%% mutation data on background  + histogram in small window
figname2 = 'mutations_and_histone_methylation_mutation_loci.pdf';
F2 = figure('Color','w', 'Name', figname2, 'Position', [682, 281, 1013, 420]);
subplot(1,2,1)
scatter(dC.all(:,3), dC.all(:,4), 'Marker','.','MarkerEdgeColor', [.5,.5,.5],'DisplayName','genome')
hold on
x = meth.te(cellfun(@(S) strcmp(S(1), 'P'), inste.sample), :);
scatter(x(:,1), x(:,2), 'Marker', '^', 'MarkerEdgeColor', 'none', ...
    'MarkerFaceColor', samplecolors(1,:), 'DisplayName', 'Plant-TE')
x = meth.te(cellfun(@(S) strcmp(S(1), 'Y'), inste.sample), :);
scatter(x(:,1), x(:,2), 'Marker', 'o', 'MarkerEdgeColor', 'none', ...
    'MarkerFaceColor', samplecolors(2,:), 'DisplayName', 'Complete-TE')
x = meth.te(cellfun(@(S) strcmp(S(1), 'M'), inste.sample), :);
scatter(x(:,1), x(:,2), 'Marker', 's', 'MarkerEdgeColor', 'none', ...
    'MarkerFaceColor', samplecolors(3,:), 'DisplayName', 'Minimal-TE')
xlabel('H3K4me2')
ylabel('H3K27me3')
xlim([0,100])
ylim([0,100])
legend

subplot(1,2,2)
scatter(dC.all(:,3), dC.all(:,4), 'Marker','.','MarkerEdgeColor', [.5,.5,.5],'DisplayName','genome')
hold on
x = meth.snp(cellfun(@(S) strcmp(S(1), 'P'), snps.sample), :);
scatter(x(:,1), x(:,2), 'Marker', '^', 'MarkerEdgeColor', 'none', ...
    'MarkerFaceColor', samplecolors(1,:), 'DisplayName', 'Plant-SNP')
x = meth.snp(cellfun(@(S) strcmp(S(1), 'Y'), snps.sample), :);
scatter(x(:,1), x(:,2), 'Marker', 'o', 'MarkerEdgeColor', 'none', ...
    'MarkerFaceColor', samplecolors(2,:), 'DisplayName', 'Complete-SNP')
x = meth.snp(cellfun(@(S) strcmp(S(1), 'M'), snps.sample), :);
scatter(x(:,1), x(:,2), 'Marker', 's', 'MarkerEdgeColor', 'none', ...
    'MarkerFaceColor', samplecolors(3,:), 'DisplayName', 'Minimal-SNP')

xlabel('H3K4me2')
ylabel('H3K27me3')
xlim([0,100])
ylim([0,100])
legend

print_pdf(F2, figname2)

%% mutation data on background, combined
figname3 = 'mutations_and_histone_methylation_mutation_loci_combined.pdf';
F3 = figure('Color','w', 'Name', figname3,'Position', [897, 348, 542, 420]);

scatter(dC.all(:,3), dC.all(:,4), 'Marker','.','MarkerEdgeColor', [.5,.5,.5])
hold on

x = meth.te(inste.comp == 0,:);
scatter(x(:,1),x(:,2), 'Marker' ,'^','MarkerEdgeColor', 'none', 'MarkerFaceColor','r','DisplayName','Plant-TE')
x = meth.snp(snps.comp == 0,:);
scatter(x(:,1),x(:,2), 'Marker' ,'^','MarkerEdgeColor', 'none', 'MarkerFaceColor','b','DisplayName','Plant-SNP')

x = meth.te(inste.comp == 1,:);
scatter(x(:,1),x(:,2), 'Marker' ,'o','MarkerEdgeColor', 'none', 'MarkerFaceColor','r','DisplayName','Complete-TE')
x = meth.snp(snps.comp == 1,:);
scatter(x(:,1),x(:,2), 'Marker' ,'o','MarkerEdgeColor', 'none', 'MarkerFaceColor','b','DisplayName','Complete-SNP')

x = meth.te(inste.comp == 2,:);
scatter(x(:,1),x(:,2), 'Marker' ,'s','MarkerEdgeColor', 'none', 'MarkerFaceColor','r','DisplayName','Minimal-TE')
x = meth.snp(snps.comp == 2,:);
scatter(x(:,1),x(:,2), 'Marker' ,'s','MarkerEdgeColor', 'none', 'MarkerFaceColor','b','DisplayName','Minimal-SNP')

xlabel('H3K4me2')
ylabel('H3K27me3')
xlim([0,100])
ylim([0,100])

print_pdf(F3, figname3)

%%
isH3K27me3.te = meth.te(:,2) > meth.te(:,1);
isH3K27me3.snp = meth.snp(:,2) > meth.snp(:,1);
isH3K27me3.all = dC.all(:,4) > dC.all(:,3);
disp(['Total number of genomic locations with H3K27me2:    ', ...
    num2str(sum(isH3K27me3.all))]);
disp(['Total number of genomic locations without H3K27me2: ', ...
    num2str(sum(~isH3K27me3.all))]);
disp(['Number of TE insertion loci with H3K27me2:          ', ...
    num2str(sum(isH3K27me3.te))]);
disp(['Number of SNP&INDEL loci with H3K27me3:             ', ...
    num2str(sum(isH3K27me3.snp))]);

%% histogram of h3k27me3
figname4 = 'mutations_and_histone_methylation_histogram.pdf';
F4 = figure('Color','w', 'Name', figname4, 'Position', [666, 316, 972, 420]);
subplot(1,2,1)

histogram(dC.all(:,4), 'Normalization', 'probability', 'BinWidth', 2.5, ...
    'EdgeColor', 'none', 'DisplayName','All Genomic')
hold on
histogram(meth.te(:,2), 'Normalization', 'probability', 'BinWidth',2.5, ...
    'EdgeColor', 'none', 'DisplayName', 'TE insertion Sites')
xlabel('H3K27me3 levels')
ylabel('Probability')
legend

subplot(1,2,2)
histogram(dC.all(:,4), 'Normalization', 'probability', 'BinWidth', 2.5, ...
    'EdgeColor', 'none', 'DisplayName', 'All Genomic')
hold on
histogram(meth.snp(:,2), 'Normalization', 'probability', 'BinWidth',2.5, ...
    'EdgeColor', 'none', 'DisplayName', 'SNP&INDEL sites')
legend
xlabel('H3K27me3 levels')
ylabel('Probability')

print_pdf(F4, figname4);

%% stats
g = length( dC.all(:,4));
g0 = sum( dC.all(:,4) == 0);
t = length( meth.te(:,2));
t0 = sum(meth.te(:,2) == 0);
s = length( meth.snp(:,2));
s0 = sum(meth.snp(:,2) == 0);
disp('One-tail Fisher Exact Test')
testt = [g0, t0; g-g0, t-t0];
[h,p,stats] = fishertest(testt,'Tail','right');
disp('H3K27me3 = 0 and TE insertion happens:')
disp(['H: ' num2str(h) ', p-value: ' num2str(p)])
tests = [g0, s0; g-g0, s-s0];
[h,p,stats] = fishertest(tests,'Tail','right');
disp('H3K27me3 = 0 and SNP or INDEL happens:')
disp(['H: ' num2str(h) ', p-value: ' num2str(p)])