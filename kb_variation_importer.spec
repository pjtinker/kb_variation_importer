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
        genome_ref: KBaseGenomes.Genome-8.3 object reference
        staging_file_subdir_path: path to VCF in staging area

        optional params:
        will_perform_gwas: groom data output for EMMAX association.

    */
    
    typedef structure {
        string workspace_name;
        obj_ref genome_ref;
        string staging_file_subdir_path;
        string command_line_args;
        boolean will_perform_gwas;
    } import_snp_params;

    typedef structure {
        string report_name;
        string report_ref;
        obj_ref variation_ref;
    } snp_import_results;


    funcdef import_snp_data(import_snp_params) 
        returns (snp_import_results) authentication required;

};
