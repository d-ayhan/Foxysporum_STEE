% print data sheets from dataCell mat files of ChrCovCal.m
load('ChrCov_dataCell_ws10000.mat')
dataCell;

a=importdata('Fol4287_GCA_003315725_genomic.fna.fai');
samples = {'WT','P1','P2','P3','P4','P5',...
    'Y1','Y2','Y3','Y4','Y5','M1','M2','M3','M4','M5'};


for k = 1:2
    if k == 1
        s = samples;
    else
        s =samples(2:end);
    end
    
    C = [{'pos_range'}, s];
    
    for i = 1:length(a.textdata)
        [t,~] = size(dataCell{i,k});
        tempC = [cellstr(...
            [char(repmat(a.textdata(i,1),t,1)),...
            repmat(char({':'}),t,1), num2str(dataCell{i,k}(:,1)),...
            repmat(char({'-'}),t,1), num2str(dataCell{i,k}(:,2))]), ....
            num2cell(dataCell{i,k}(:,3:end))];
        C = [C; tempC];
    end
    if k == 1
        xlswrite('CovData.xlsx',C,'raw_10k')
    else
        xlswrite('CovData.xlsx',C,'normalized_10k')
    end
end


