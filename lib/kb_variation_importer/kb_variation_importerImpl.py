# -*- coding: utf-8 -*-
#BEGIN_HEADER
import os 
import subprocess
import uuid

from kb_variation_importer.Utils import variation_importer_utils
from DataFileUtil.DataFileUtilClient import DataFileUtil
from KBaseReport.KBaseReportClient import KBaseReport

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

    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.config = config
        self.shared_folder = config['scratch']
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        #END_CONSTRUCTOR
        pass

    def import_snp_data(self, ctx, import_snp_params):
        """
        :param workspace_name: instance of String
        :param import_snp_params: instance of type "import_snp_params" (Input
        parameters) -> structure: parameter "input_file_path" of String,
        parameter "will_perform_gwas" of Int
        :returns: instance of type "import_snp_results" (Output results) ->
        structure:
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN import_snp_data
        variation_utils = variation_importer_utils.variation_importer_utils(STORAGE_DIR)
        print("Params passed to import_snp_data: {}".format(import_snp_params))
        vcf_version = None
        try:
            vcf_version = variation_utils.generate_vcf_stats(import_snp_params)
        except Exception as e:
            print("Error importing variation data!")
            raise ValueError(e)

        file_extensions = ['frq', 'log']
        indexHTML = "<head><body> "
        #for ext in file_extensions:
        indexHTML += "<a href='./frequencies.frq'>Frequencies</a>"
        indexHTML += "<a href='./frequencies.log'>Log</a>"
        indexHTML += "</body></head>"
        # indexHTML += """
        #         <a href="./frequencies.frq">Frequencies</a>
        #         <a href="./frequencies.log">Log</a>
        #     </body>
        # </head>
        # """
        with open(STORAGE_DIR + '/index.html', 'w') as html:
            html.write(str(indexHTML))
        dfu = DataFileUtil(self.callback_url)

        try:
            html_upload_ret = dfu.file_to_shock({'file_path': STORAGE_DIR, 'make_handle': 0, 'pack': 'zip'})
        except:
            raise ValueError('Error uploading HTML to shock')

        report_name = 'vcf_summary_stats_' + str(uuid.uuid4())
        reportObj = {'objects_created': [],
                     'message': '',
                     'direct_html': None,
                     'direct_html_index': 0,
                     'file_links': [],
                     'html_links': [{'shock_id': html_upload_ret['shock_id'],
                                    'name': 'index.html',
                                    'label': 'View generated files'
                                    }
                                   ],
                     'html_window_height': 220,
                     'workspace_name': import_snp_params['workspace_name'],
                     'report_object_name': report_name
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
