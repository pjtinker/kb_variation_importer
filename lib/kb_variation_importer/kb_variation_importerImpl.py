# -*- coding: utf-8 -*-
#BEGIN_HEADER
import os 
import subprocess
import uuid
import errno

from kb_variation_importer.Utils import variation_importer_utils

from DataFileUtil.DataFileUtilClient import DataFileUtil
from KBaseReport.KBaseReportClient import KBaseReport
# from GenomeAnnotationAPI.GenomeAnnotationAPIServiceClient import GenomeAnnotationAPI
class InvalidVCFError(Exception):
    def __init__(self, file_path, message):
        self.file_path = file_path
        self.message = message
STORAGE_DIR = "/kb/module/work/tmp/html/"

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
    GIT_URL = ""
    GIT_COMMIT_HASH = ""

    #BEGIN_CLASS_HEADER
    def _mkdir_p(self, path):
        """
        _mkdir_p: make directory for given path
        """
        if not path:
            return
        try:
            os.makedirs(path)
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
        self.scratch_file_path = os.path.join(config['scratch'], 'import_variation_' + str(uuid.uuid4()) + '/')
        self._mkdir_p(self.scratch_file_path)
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.dfu = DataFileUtil(self.callback_url)
        self.vu = variation_importer_utils.variation_importer_utils(self.config['srv-wiz-url'], self.dfu, self.scratch_file_path)

        #END_CONSTRUCTOR
        pass

    def import_snp_data(self, ctx, import_snp_params):
        """
        :param workspace_name: instance of String
        :param import_snp_params: instance of type "import_snp_params" (Input
        parameters) -> structure: parameter "staging_file_subdir_path" of String,
        parameter "will_perform_gwas" of Int
        :returns: instance of type "import_snp_results" (Output results) ->
        structure:
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN import_snp_data

        print("Params passed to import_snp_data: {}".format(import_snp_params))
        vcf_version = None

        # This is the process if staging file rights have been granted
        # vcf_scratch_path = self.dfu.download_staging_file(
        #     {'staging_file_subdir_path': import_snp_params['staging_file_subdir_path']
        # }).get('copy_file_path')
        vcf_scratch_path = self.vu.pretend_download_staging_file(
            import_snp_params['staging_file_subdir_path'], self.scratch_file_path).get('copy_file_path')
        
        print("Scratch file path produced by DFU: {}".format(vcf_scratch_path))
        genome_ref='18590/2/8'

        try:
            vcf_version = self.vu.validate_vcf(vcf_scratch_path, genome_ref)
        except Exception as e:
            print("Error importing variation data!")
            print(e)
            raise ValueError(e)

        try:
            self.vu.generate_vcf_stats(import_snp_params['command_line_args'], vcf_scratch_path, genome_ref)
        except Exception as e:
            print("Error generating summary statistics!")
            print(e)
            raise

        file_extensions = ['frq', 'log']
        indexHTML = "<head><body> "
        #TODO: Create links for all relevant files generated. 
        indexHTML += "<a href='./frequencies.frq'>Frequencies</a> "
        indexHTML += "<a href='./frequencies.log'>Log</a>"
        indexHTML += "</body></head>"
        
        with open(os.path.join(self.scratch_file_path, 'index.html'), 'w') as html:
            html.write(str(indexHTML))

        try:
            html_upload_ret = self.dfu.file_to_shock({'file_path': self.scratch_file_path, 'make_handle': 0, 'pack': 'zip'})
            print("File to shock info: {}".format(html_upload_ret))
        except:
            raise ValueError('Error uploading HTML to shock')

        reportObj = {'objects_created': [],
                     'message': '',
                     'direct_html': None, #Is this where to include html displayed in narrative?
                     'direct_html_index': 0,
                     'file_links': [],
                     'html_links': [{
                                    'shock_id': html_upload_ret['shock_id'],
                                    'name': 'index.html',
                                    'label': 'View generated files'
                                    }
                                   ],
                     'html_window_height': 220,
                     'workspace_name': import_snp_params['workspace_name'],
                     'report_object_name': 'vcf_summary_stats_' + str(uuid.uuid4())
                     }

        report = KBaseReport(self.callback_url, token=ctx['token'])
        report_info = report.create_extended_report(reportObj)
        returnVal = { 'report_name': report_info['name'], 
                      'report_ref': report_info['ref'],
                      'vcf_version' : vcf_version
                     }  
        
        if not isinstance(returnVal, dict):
            raise ValueError('Method import_snp_data return value ' +
                             'returnVal is not type dict as required.')
        #END import_snp_data

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
