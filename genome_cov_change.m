% compare genomic covarage of evolved popualtions
%faifile = 'Fol4287_GCA_003315725_genomic.fna.fai';
%a=importdata(faifile);

%contigCount= length(a.data);
%[~,loc] = ismember(samples,wgscovs.textdata);



load('format_genome_file.mat')
window_size = 10000;
samples = {'WT','P1','P2','P3','P4','P5',...
    'Y1','Y2','Y3','Y4','Y5','M1','M2','M3','M4','M5'};

% wgscovs=importdata('wgs_cov_stats.txt');
% [~,loc] = ismember(samples,wgscovs.textdata);
% MedCovs = wgscovs.data(loc);

for i = 1: length(new_f)
    F{i} = zeros(length(new_f(i).Sequence),16);
end

for k=1:16
    T= readtable(['genomecov\' samples{k} '.genomecov.txt']);
    for i = 1:length(reordered_a.textdata)
        Lia = strcmp(reordered_a.textdata{i},T.Var1);
        ind = new_f_indexing(i,:);
        F{ind(1,2)}(ind(1,1)+1:ind(1,1)+sum(Lia),k) = T.Var3(Lia);
    end
end
clearvars T
%
%
% readcounts.data(1,1) = readcounts.data(1,1) + readcounts.data(2,1);
% readcounts.data(2,1) = readcounts.data(3,1) + readcounts.data(4,1);
% readcounts.data(3:4,:) = [];
%

save('coverage_for_circos.mat','F', '-v7.3')

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
        varsum(end+1,:) = median(F{i}(j:k-1,:));
    end
end

Z = varsum ./ median(varsum,1);
Y = Z ./ Z(:,1);
W = mean(Y(:,2:16),2);
V = std(Y(:,2:16),[],2);
R = var(Y(:,2:16),[],2);
subplot(3,1,1)
plot(W)
subplot(3,1,2)
plot(V)
subplot(3,1,3)
plot(R)
T = table(chr, position, W);
writetable(T,  'cov_change_for_circos.txt', 'Delimiter','\t', 'WriteVariableNames',0)
