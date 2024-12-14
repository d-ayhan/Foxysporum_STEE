load('format_variants_file.mat')
load('mutation_dynamics.mat')
samp = {'P2','Y3','M4'};

for z = 1:2
    if z == 1
        typ = 'fixed';
        zz = 1;
    else
        typ = 'extinct';
        zz = 0;
    end
    for i = 1:3 %samp
        M.(samp{i}).(typ) = DD.(samp{i}).mutid(DD.(samp{i}).AF(:,end) == zz);
        
        for j = 1:length(M.(samp{i}).(typ))
            mut = M.(samp{i}).(typ){j,1};
            if startsWith(mut, 'T')
                t = inste;
            else
                t = snps;
            end
            
            ind = find(ismember(t.mut_id, mut), 1);
            
            if startsWith(t.CHROM(ind), 'Chr') || strcmp(t.CHROM(ind), 'U_2') || strcmp(t.CHROM(ind), 'U_3')
                M.(samp{i}).(typ){j,2} = 'core';
            else
                M.(samp{i}).(typ){j,2} = 'ls';
            end
            
            if strcmp(t.loc(ind),'coding')
                M.(samp{i}).(typ){j,3} = 'coding';
            else
                M.(samp{i}).(typ){j,3} = 'non-coding';
            end
            
            if startsWith(mut, 'T')
                M.(samp{i}).(typ){j,4} = 'T';
            else
                M.(samp{i}).(typ){j,4} = 'S';
            end
        end
    end
end

for k = 1:3 %groups
    switch k
        case 1
            a = 'core';
            b = 'ls';
        case 2
            a = 'coding';
            b = 'non-coding';
        case 3
            a = 'T';
            b = 'S';
    end
    for i = 1:3 %sample
        
        for j = 1:2 %fixed extinct
            switch j
                case 1
                    typ = 'fixed';
                case 2
                    typ = 'extinct';
            end
            
            D(k).A(i,j) = sum(ismember(M.(samp{i}).(typ)(:,k+1),a));
            D(k).B(i,j) = sum(ismember(M.(samp{i}).(typ)(:,k+1),b));
            
        end
    end
end

figname = 'extinct_vs_fixed.pdf';
f = figure('Name',figname,'Color','w','Position',[576, 125, 664, 853]);
for i = 1:3
    switch i
        case 1
            a = 'Core';
            b = 'LS';
        case 2
            a = 'Coding';
            b = 'Non-coding';
        case 3
            a = 'TE insertion';
            b = 'SNPs and INDELs';
    end
    ax(i,1) = subplot(3,2,2*i-1);
    bar(ax(i,1), D(i).A);
    ax(i,1).XTick = 1:3;
    ax(i,1).XTickLabels = samp;
    ax(i,1).YLim = [0 15];
    ax(i,1).Title.String = a;
    
    
    ax(i,2) = subplot(3,2,2*i);
    bar(ax(i,2), D(i).B);
    ax(i,2).XTick = 1:3;
    ax(i,2).XTickLabels = samp;
    ax(i,2).YLim = [0 15];
    ax(i,2).Title.String = b;
end
legend(ax(3,2),{'Fixed', 'Extinct'})

print_pdf(f, figname)