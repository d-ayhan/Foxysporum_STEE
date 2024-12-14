load('ChrCov_dataCell_ws10000.mat')
a=importdata('Fol4287_GCA_003315725_genomic.fna.fai');

samples = {'WT','P1','P2','P3','P4','P5',...
    'Y1','Y2','Y3','Y4','Y5','M1','M2','M3','M4','M5'};

sections = [1,2,3,4,5,6,7,8,9,10,13,19,151,197,239]';

[coverage, order, ~] = xlsread('clustered_scaffolds_with_coverages.xlsx');
order = order(:,1);

[Lia,Loc] = ismember(order,a.textdata);

%reshape dataCell according to order
D = dataCell(Loc(Lia),2);

contigbinlengths = cellfun('size',D(:,1),1);
contigbinlengths(:,2) = cumsum(contigbinlengths(:,1));
xaxisticks = [contigbinlengths(sections,2)]+.5;


for i = 1:length(D)
    D{i}(:,3:end) = D{i}(:,3:end).*coverage(i);
    D{i}(:,2) = ones(size(D{i},1),1).*coverage(i);
end

%fix chr 13 end 
D{11,1}(137:end,3:17) = D{11,1}(137:end,3:17).*2;
D{11}(137:end,2) = ones(45,1).*2;

DD = cell2mat(D);
DD(:,1) = [];

%%
CO = [.2,.2,.2; repmat([47,130,69],5,1);repmat([163,50,191],5,1);repmat([66,78,179],5,1)]./255; %black, 5green, 5purple, 5blue

x = 1:contigbinlengths(end,2);
y2 = zeros(size(x));
x2 = [x, fliplr(x)];
%% all
figname1 = 'coverage_plot_all.pdf';
f1 = figure('Color', 'w', 'Name', figname1,  'Renderer', 'Painters', 'Position', [1, 41, 1920, 963]);
ax = axes;
for i = 1:16
    y = DD(:,i);
    y1p = y'+((16-i).*5);
    y2p = y2+((16-i).*5);
    h = line(ax,x, y1p, 'Color', 'none','LineWidth',1);
    hold on
    h2 = line(ax,x, y2p, 'Color', 'none','LineWidth',1);
    inbetween = [y1p, y2p];
    fill(x2, inbetween,CO(i,:),'EdgeColor','none')
    hold on
end
ax.Color = 'none';
temp = reshape(1:80,[5,16]);
temp = temp(1:3,:);
ax.YTick = reshape(temp,1,[]);
ax.YGrid='on';
ax.XGrid ='on';
ax.Box='on';
ax.XLim=[1 x(end)];
ax.YLim=[0 80];
ax.YTickLabel=[];
ax.XTick=xaxisticks;
ax.XTickLabel = {'1','2','4','5','7','8','9','10','11','12','13','2a','3','14','15'};
ax.YLabel.String='Normalized Coverage';
ax.XLabel.String='Position (Mb)';

print_pdf(f1, figname1)

%% subset
%remove P3 
DD(:,4) = [];
CO(4,:) = [];

%remove core chr
remv = xaxisticks(1)+.5 : xaxisticks(10)-0.5;
DD(remv,:) = [];
x = 1:size(DD,1);
y2 = zeros(size(x));
x2 = [x, fliplr(x)];
xaxisticks = [xaxisticks(1);  xaxisticks(11:end) - xaxisticks(10) + xaxisticks(1)];

%
figname2 = 'coverage_plot_subset.pdf';
f2 = figure('Color', 'w', 'Name', figname1, 'Renderer', 'Painters', 'Position', [680, 70, 560, 908]);
ax = axes;

for i = 1:15
    y = DD(:,i);
    y1p = y'+((15-i).*5);
    y2p = y2+((15-i).*5);
    h = line(ax,x, y1p, 'Color', 'none','LineWidth',1);
    hold on
    h2 = line(ax,x, y2p, 'Color', 'none','LineWidth',1);
    inbetween = [y1p, y2p];
    fill(x2, inbetween,CO(i,:),'EdgeColor','none')
    hold on
end
ax.Color = 'none';
temp = reshape(1:75,[5,15]);
temp = temp(1:3,:);
ax.YTick = reshape(temp,1,[]);
ax.YGrid='on';
ax.XGrid ='on';
ax.Box='on';
ax.XLim=[1 x(end)];
ax.YLim=[0 75];
ax.YTickLabel=[];
ax.XTick=xaxisticks;
ax.XTickLabel = {'1','13','2a','3','14','15'};
ax.YLabel.String='Normalized Coverage';
ax.XLabel.String='Position (Mb)';

print_pdf(f2, figname2)
