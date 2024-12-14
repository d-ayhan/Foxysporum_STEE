% format genomecoverage file

a=importdata('Fol4287_GCA_003315725_genomic.fna.fai');

files = {'H3K4me2_1.genomecov.txt','H3K4me2_2.genomecov.txt','H3K27me3_1.genomecov.txt','H3K27me3_2.genomecov.txt'};
windowSize = 1000;
samples = {'H3K4me2_1','H3K4me2_2','H3K27me3_1','H3K27me3_2'};
contigCount= length(a.data);

readcounts = importdata('histone_methylation_total_read_counts.txt');

dataCell = cell(contigCount,2); % (:,1) -> for raw; (:,2) -> for normalized
for i = 1:contigCount
    t = a.data(i,1);
    tt = 1:t;
    m = tt(1:windowSize:end)'; % goes like: 1,1001,2001...
    n = m+windowSize-1; % range of the windows (m:n)
    n(end) = t; % last point is lat pos.
    
    dataCell{i,1} = NaN(length(m),1 + 2); % first 2 are the range, then the wt&samples
    dataCell{i,2} = NaN(length(m),1 + 1); % first 2 are the range, then the samples
    dataCell{i,1}(:,1:2) = [m,n];
    dataCell{i,2}(:,1:2) = [m,n];
end
pos = cell(1,contigCount);

for j = 1:length(files)
    cov = LoadBedtoolsGenomecovFile(files{j},...
        'Fol4287_GCA_003315725_genomic.fna.fai',windowSize,'median');
    for i=1:contigCount
        A = dataCell{i,1};
        pos{j,i} = floor(A(:,1) + (A(:,2) - A(:,1) +1)/2);
        dataCell{i,1}(:,j+2) = cov{i,1}(pos{j,i});
    end
end

for i=1:length(dataCell)
    dataCell{i,1}(:,3) =  dataCell{i,1}(:,3) + dataCell{i,1}(:,4);
    dataCell{i,1}(:,4) =  dataCell{i,1}(:,5) + dataCell{i,1}(:,6);
    dataCell{i,1}(:,5:6) = [];
end

save(['histone_methylation_covdata' num2str(windowSize) '.mat'],'dataCell')