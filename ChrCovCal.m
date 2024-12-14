% compare genomic covarage of evolved popualtions
faifile = 'Fol4287_GCA_003315725_genomic.fna.fai';
a=importdata(faifile);
wgscovs=importdata('wgs_cov_stats.txt');
files = dir('genomecov\*.txt');
windowSize = 10000;
fileCount = length(files);
samples = {'WT','P1','P2','P3','P4','P5',...
    'Y1','Y2','Y3','Y4','Y5','M1','M2','M3','M4','M5'};
contigCount= length(a.data);
[~,loc] = ismember(samples,wgscovs.textdata);
MedCovs = wgscovs.data(loc);

% generate final dataCell with positions of the windows inside
dataCell = cell(contigCount,2); % (:,1) -> for raw; (:,2) -> for normalized
for i = 1:contigCount
    t = a.data(i,1);
    tt = 1:t;
    m = tt(1:windowSize:end)'; % goes like: 1,1001,2001...
    n = m+windowSize-1; % range of the windows (m:n)
    n(end) = t; % last point is lat pos.
        
    dataCell{i,1} = NaN(length(m),fileCount + 2); % first 2 are the range, then the wt&samples
    dataCell{i,2} = NaN(length(m),fileCount + 1); % first 2 are the range, then the samples
    dataCell{i,1}(:,1:2) = [m,n];
    dataCell{i,2}(:,1:2) = [m,n];
end
pos = cell(fileCount,contigCount);
for k = 1:fileCount
    kk{k} = files(k).name(1:2);
end
[~,l] = ismember(samples,kk);
for j = 1:fileCount
    cov = LoadBedtoolsGenomecovFile([files(l(j)).folder '\' files(l(j)).name], ...
        faifile, windowSize, 'median');
    for i=1:contigCount
        A = dataCell{i,1};
        pos{j,i} = floor(A(:,1) + (A(:,2) - A(:,1) +1)/2);
        dataCell{i,1}(:,j+2) = cov{i,1}(pos{j,i})./MedCovs(j);
        if j ~= 1 % for the samples, do normalization
            dataCell{i,2}(:,j+1) = dataCell{i,1}(:,j+2)./dataCell{i,1}(:,1+2);
        end
    end
end

save(['ChrCov_dataCell_ws' num2str(windowSize) '.mat'],'dataCell')