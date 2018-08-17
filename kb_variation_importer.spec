/*
A KBase module: kb_variation_importer
*/

module kb_variation_importer {

    /* A boolean - 0 for false, 1 for true.
        @range (0, 1)
    */
    typedef int boolean;

    /* An X/Y/Z style reference
    */
    typedef string obj_ref;

    /*
        required params:
        genome_ref: KBaseGenomes.Genome object reference
        staging_file_subdir_path: path to VCF in staging area
        location_file_subdir_path: path to location file in staging area.

        optional params:


        *** Filtering ***
        maf_threshold: percent threshold for filtering by maf
        geno_missingness: percent threshold to exclude SNPs based on missing calls
        indiv_missingness: percent threshold to remove individuals with missing calls
        hwe_threshold: p-value threshold to remove samples not in Hardy-Weinberg Equilibrium

        *** Visualization ***
        plot_maf: generate histogram of minor allele frequencies
        plot_hwe: generate histogram of Hardy-Weinberg Equilibrium p-values

        will_perform_gwas: groom data output for EMMAX association.


    */
    
    typedef structure {
        string workspace_name;
        obj_ref genome_ref;
        string staging_file_subdir_path;
        string location_file_subdir_path;
        string command_line_args;
        boolean will_perform_gwas;
    } import_variation_params;

    typedef structure {
        string report_name;
        string report_ref;
        obj_ref variation_ref;
    } import_variation_results;


    funcdef import_variation(import_variation_params) 
        returns (import_variation_results) authentication required;

};
