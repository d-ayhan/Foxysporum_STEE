% Main scripts for Lopez-Diaz et al., 2020
% Written by Ayhan, DH 2020/06

%% Section 1 experimental evo method & phenotype

%% Section 2 genotype
    % Table: all mutations
    format_variants_file
    % allele frequencies histograms
    allele_frequencies
    % type of mutations in all vs type of mutations in 10th
    type_of_mutations
    % TE table
    active_TEs
    % methylation bias
    histone_methylation_data
    mutations_and_histone_methylation
    % circos plot
    %%% Layer 0: karyotype
    format_genome_file
    %%% Layer 1: GC content
    genome_GC_content
    %%% Layer 2: TE content
    genome_TE_content
    %%% Layer 3-4: methylation
    chr_map_hismeth
    %%% Layer 5-6: SNPs and TEinsertions
    chr_map_variants
    %%% Layer 7: Coverage change
    genome_cov_change
    % GO term enrichment for mutated genes

%% Section 3 intermediates
    % AF change in passages
    mutation_dynamics
    % extinct vs fixed
    extinct_VS_fixed
    % muller plot
    
%% Section 4 velvet case

%% Section 5 CNV
    % change over all    
    ChrCovCal
    ChrCovPrint
    clusterCovData
    coverage_plot
    