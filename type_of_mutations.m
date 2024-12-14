load('format_variants_file.mat')

conditions.short = {'P', 'Y' ,'M'};
conditions.long = {'Plant', 'Rich media', 'Minimal media'};

% filter table
AFcutOff = 0.1;
snps(snps.AF < AFcutOff, :) = [];
inste(inste.AF < AFcutOff, :) = [];

% sort into categories
vn = {'sample','mut_id','AF','chrom','loc','dist','type'};
for i = 1:3
    for j = 1:5
        if i == 1 && j == 3 %ignore P3
            continue 
        end
        sn = [conditions.short{i}, num2str(j)];
        idx = cellfun(@(S) strcmp(S(1:2),sn), snps.sample);
        T1 = table(snps.sample(idx), snps.mut_id(idx), snps.AF(idx), snps.CHROM(idx), ...
            snps.loc(idx), snps.dist(idx), snps.TYPE(idx), 'VariableNames', vn);
        
        idx = cellfun(@(S) strcmp(S(1:2),sn), inste.sample);
        T2 = table(inste.sample(idx), inste.mut_id(idx), inste.AF(idx), inste.CHROM(idx), ...
            inste.loc(idx), inste.dist(idx), inste.TE(idx), 'VariableNames', vn);
        
        data.(sn) = [T1;T2];
    end
end

%% summarize
samples = cell(15,1);
finalpassages = nan(15,3);
allpassages = nan(15,3);
muttype = {'T','S','I'};
for i = 1:3 %muttype
    for j = 1:3 %conditions
        for k = 1:5 %replicates
            if j == 1 && k == 3 %ignore P3
                continue
            end
            sn = [conditions.short{j}, num2str(k)];
            finalsample = [sn '-10'];
            Ta = data.(sn);
            Tf = Ta(ismember(Ta.sample, finalsample), :);
            Tam = unique(Ta.mut_id);
            Tfm = unique(Tf.mut_id);
            allpassages((j-1).*5+k,i) = sum(cellfun(@(S) strcmp(S(1),muttype{i}), Tam));
            finalpassages((j-1)*5+k,i) = sum(cellfun(@(S) strcmp(S(1),muttype{i}), Tfm));
            samples{(j-1)*5+k} = finalsample;
        end
    end
end

T = table(samples, finalpassages(:,1), finalpassages(:,2), finalpassages(:,3), 'VariableNames', {'Samples', 'TEins', 'SNPs', 'INDELs'});
T(3,:) = [];

writetable(T,'mutation_summary_final_passage.txt')
%% arrange the data for bar graph
AP = [allpassages(1:5,:); nan(1,3); allpassages(6:10,:); nan(1,3); allpassages(11:15,:)];
AP(3,:) = [];
FP = [finalpassages(1:5,:); nan(1,3); finalpassages(6:10,:); nan(1,3); finalpassages(11:15,:)];
FP(3,:) = [];

%%
figname = 'type_of_mutations.pdf';
h = figure('Name',figname,'Color','w','Position',[680, 349, 1074, 629]);

ax(1) = subplot(2,2,1); % raw counts, final passage
bar(ax(1),FP,1,'stacked')
title('Raw Counts - Final Passage')
legend('TE insertion','SNPs', 'Indels');

ax(2) = subplot(2,2,2); % raw counts, final passage
bar(ax(2),AP,1,'stacked')
title('Raw Counts - All Passages')


ax(3) = subplot(2,2,3); % raw counts, final passage
bar(ax(3),FP./sum(FP,2),1,'stacked')
title('Ratios - Final Passage')


ax(4) = subplot(2,2,4); % raw counts, final passage
bar(ax(4),AP./sum(AP,2),1 ,'stacked')
title('Ratios - All Passages')

for i = 1:4
    ax(i).XLim = [0 16];
    ax(i).XTick = [2.5, 8, 14];
    ax(i).XTickLabel = conditions.long;
end
print_pdf(h, figname);