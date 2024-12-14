% GC content
load('format_genome_file.mat')
window_size = 10000;
chr = {};
position=[];
GCcont = [];
for i = 1: length(new_f)
    seq = new_f(i).Sequence;
    l = length(seq);
    for j = 1:window_size:l
        k = window_size + j;
        if k > l
            k = l;
        end
        
        base = basecount(seq(j:k));
        gc = (base.G + base.C)./(base.A + base.T + base.G + base.C);
        
        chr{end+1,1} = new_f(i).Header;
        position(end+1,1:2) = [j-1,k-1];
        GCcont(end+1,1) = gc;
    end
end
        
T = table(chr, position, GCcont);
writetable(T, 'GC_densitiy.txt', 'Delimiter','\t', 'WriteVariableNames',0)