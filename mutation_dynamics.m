create_color_arrays
load('format_variants_file.mat')
type = {'TE', 'inste'; 'SNP', 'snps'};
samp = {'P2','Y3','M4'};


for j = 1:2 %snp or te
    if j == 1
        t = inste;
    else
        t = snps;
    end
    for i = 1:3 %samples
        T.(type{j,1}).(samp{i}) = t(contains(t.sample,samp{i}),:);
        muts = unique(T.(type{j,1}).(samp{i}).mut_id);
        mut.(samp{i}).(type{j,1}) = muts;
        D.(samp{i}).(type{j,1}) = zeros(length(muts),11);
        for z = 1:length(muts)
            mutid = muts{z};
            temp = T.(type{j,1}).(samp{i})(ismember(T.(type{j,1}).(samp{i}).mut_id, mutid),:);
            for k = 1:10
                lia = endsWith(temp.sample, ['-' num2str(k)]);
                if any(lia)
                    D.(samp{i}).(type{j,1})(z,k+1) = temp.AF(lia);
                end
            end
        end
    end
end

clearvars j t i muts z mutid temp k lia

for i = 1:3
    DD.(samp{i}).AF = [D.(samp{i}).TE; D.(samp{i}).SNP];
    DD.(samp{i}).mutid = [mut.(samp{i}).TE; mut.(samp{i}).SNP];
end
save('mutation_dynamics.mat','DD')

% filter out the mutations with max AF < 0.1
for i = 1:3
    if i ~= 1
        DD.(samp{i}).AF(:,[3:5 7:9]) = [];
    end
    f = max(DD.(samp{i}).AF,[],2) < 0.1;
    DD.(samp{i}).AF(f,:) = [];
    DD.(samp{i}).mutid(f) = [];
end

%% plots
figname = 'mutation_dynamics.pdf';
f = figure('Color','w', 'Name', figname, 'Position', [353, 42, 863, 954]);

for i = 1:3
    ax(i) = subplot(3,1,i);
    ax(i).Title.String = samp{i};
    
    if i == 1
        x = 0:10;
    else
        x = [0 1 5 9 10];
    end
    mut = DD.(samp{i}).mutid;
    for j = 1:length(mut)
        lh(j) = line(ax(i), x, DD.(samp{i}).AF(j,:),'Color',color_qual(j,:),'DisplayName',mut{j});
        hold on
        if startsWith(mut{j}, 'T')
            lh(j).Marker = 'v';
            t = inste;
        else
            lh(j).Marker = 'o';
            t = snps;
        end
        ind = find(ismember(t.mut_id, mut{j}),1);
        if startsWith(t.CHROM(ind), 'Chr') || strcmp(t.CHROM(ind), 'U_2') || strcmp(t.CHROM(ind), 'U_3')
            lh(j).LineStyle = '-';
        else
            lh(j).LineStyle = '--';
        end
        
        if strcmp(t.loc(ind),'coding')
            lh(j).LineWidth = 2;
        else
            lh(j).LineWidth = 1;
        end
    end
    
    ax(i).XLim = [0 10];
    ax(i).YLim = [0, 1.1];
    ax(i).XLabel.String = 'Passages';
    ax(i).YLabel.String = 'AF';
    
    %   if i == 2
    %       legend('Location','westoutside')
    %   else
    legend('Location','eastoutside')
    %   end
end
clearvars i ax x mut j lh t ind

print_pdf(f, figname)
