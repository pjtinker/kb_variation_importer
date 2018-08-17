# -*- coding: utf-8 -*-
#BEGIN_HEADER
import os 
import subprocess
import uuid
import errno
import json
from kb_variation_importer.Utils import variation_importer_utils

# from DataFileUtil.DataFileUtilClient import DataFileUtil
# from KBaseReport.KBaseReportClient import KBaseReport
# from GenomeAnnotationAPI.GenomeAnnotationAPIServiceClient import GenomeAnnotationAPI


#END_HEADER


class kb_variation_importer:
    '''
    Module Name:
    kb_variation_importer

    Module Description:
    A KBase module: kb_variation_importer
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "0.0.1"
    GIT_URL = "https://github.com/pjtinker/kb_variation_importer.git"
    GIT_COMMIT_HASH = "058c2667a5ea5420ebfc44611faf3cc4318cde3d"

    #BEGIN_CLASS_HEADER
    def _mkdir_p(self, path):
        """
        _mkdir_p: make directory for given path
        """
        if not path:
            return
        try:
            os.mkdir(path)
        except OSError as exc:
            if exc.errno == errno.EEXIST and os.path.isdir(path):
                pass
            else:
                raise
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.config = config
        self.scratch = config['scratch']
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        #END_CONSTRUCTOR
        pass


    def import_variation(self, ctx, import_variation_params):
        """
        :param import_variation_params: instance of type "import_variation_params"
           (Insert your typespec information here.) -> structure: parameter
           "workspace_name" of String, parameter "staging_file_subdir_path"
           of String, parameter "will_perform_gwas" of Long
        :returns: instance of type "snp_import_results" -> structure:
           parameter "report_name" of String, parameter "report_ref" of
           String, parameter "vcf_version" of String
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN import_variation

        returnVal = {}
        # TODO: Validate params

        utility_params = self.config
        utility_params['token'] = ctx['token']
        utility_params['callback_url'] = self.callback_url

        self.vu = variation_importer_utils.variation_importer_utils(utility_params)


        # This is the process if staging file rights have been granted
        # vcf_staging_area_path = self.dfu.download_staging_file(
        #     {'staging_file_subdir_path': import_variation_params['staging_file_subdir_path']
        # }).get('copy_file_path')
        # vcf_staging_area_path = self.vu.pretend_download_staging_file(
        #     import_variation_params['staging_file_subdir_path'], self.scratch).get('copy_file_path')
        
        # print("Scratch file path produced by DFU: {}".format(vcf_staging_area_path))

        try:
            returnVal = self.vu.validate_vcf(import_variation_params)
        except Exception as e:
            print("Error importing variation data!")
            raise ValueError(e)
        
        #END import_variation

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method import_variation return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]

    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
