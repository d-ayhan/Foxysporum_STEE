load('format_variants_file.mat')
Class = readtable('TE_classification.csv');

for k = 1:2
    if k == 1
        inste_e = inste(contains(inste.sample, '-10'),:);
        tablename = 'active_TEs_finalonly.txt';
    else
        inste_e = inste;
        tablename = 'active_TEs_all.txt';
    end
    
    [~, idx] = unique(inste_e.mut_id);
    TE = inste_e(idx,:);
    
    s.P = TE(contains(TE.sample, 'P'),:);
    s.Y = TE(contains(TE.sample, 'Y'),:);
    s.M = TE(contains(TE.sample, 'M'),:);
    
    family = unique(Class.Family);
    
    for i = 1:length(family)
        temp = Class.name(ismember(Class.Family,family{i}));
        
        c.P(i,1) = sum(ismember(s.P.TE, temp));
        c.Y(i,1) = sum(ismember(s.Y.TE, temp));
        c.M(i,1) = sum(ismember(s.M.TE, temp));
    end
    
    T = table(family,c.P,c.Y,c.M,'VariableNames',{'family','P','Y','M'});
    writetable(T,tablename)
end