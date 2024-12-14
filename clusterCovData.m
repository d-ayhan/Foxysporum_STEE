[data,textdata,~]=xlsread('CovData20191001.xlsx','filtered');

%columns: observations; rows: samples
D = data(:,1:15)';
%correlation matrix
rho = corr(D);
%dissimilarities
diss = 1- rho;
%hierarchical clustering
Z= linkage(diss,'complete');


figure
mysubplot(2,1,1)
%plot hierarchical clustering
[H,T,outperm] = dendrogram(Z,0,'ColorThreshold',16);
xticks([])

mysubplot(2,1,2)
%reorder by clusters
dd=D(:,outperm);
%show by heatmap
heatmap(dd)
colormap('Jet')
grid off

%generate clusters
groups = cluster(Z,'Cutoff',16,'criterion','distance');


% print bins
td=textdata(2:end,1);
Tab=table(groups(outperm),td(outperm),D(:,outperm)');
writetable(Tab,'Filtered_clustered.txt')
% % longitudinal 
% figure
% heatmap(dd') 
% colormap('Jet')
% grid off