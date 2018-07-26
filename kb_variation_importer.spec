/*
A KBase module: kb_variation_importer
*/

module kb_variation_importer {
    /*
        Insert your typespec information here.
    */
    
    typedef structure {
        string workspace_name;
        string staging_file_subdir_path;
        int will_perform_gwas;
    } import_snp_params;

    typedef structure {
        string report_name;
        string report_ref;
        string vcf_version;
    } snp_import_results;


    funcdef import_snp_data(import_snp_params) 
        returns (snp_import_results) authentication required;

};
