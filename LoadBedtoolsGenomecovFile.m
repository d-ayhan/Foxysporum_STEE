function [out] = LoadBedtoolsGenomecovFile(genomecov,genome,windowsize,outputType)
% genomecov: str : bedtools genomecov output file
% genome: str : includes chromosome names and lengths (bedtools .genome file)
% windowsize:int : to generate moving window
% outType: str :{'mat','median',''}

a=importdata(genome);
b=importdata(genomecov);
C=cell(length(a.data),1);
covMovMed=cell(length(a.data),1);
c=1;
for i=1:length(a.data)
    C{i}=b.data(c:(c+a.data(i)-1),2);
    covMovMed{i}=movmedian(C{i},windowsize);
    c=c+a.data(i);
end

if strcmp(outputType,'mat')
    save(strcat(genomecov, '.mat'),'C','covMovMed')
    out = [];
elseif strcmp(outputType,'median')
    out = covMovMed;
else
    out = C;
end

end