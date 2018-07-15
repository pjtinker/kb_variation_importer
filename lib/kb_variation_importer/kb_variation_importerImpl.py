# -*- coding: utf-8 -*-
#BEGIN_HEADER
import os 
import subprocess
import uuid
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
    vcf_version = None

    def validate_vcf(self, vcf_path):
        """
            :param vcf_path: string defining directory of VCF file
        """
        #BEGIN validate_vcf
        # TODO determine size of file.  May want to use HDF5 for
        print('\nValidating VCF...')
        try:
            f = open(vcf_path, "r")
        except Exception as e:
            print("Error opening file: {}".format(e))
            raise InvalidVCFError(vcf_path, e)
        
        line = f.readline()
        tokens = line.split('=')
        f.close()
        if(tokens[0] != "##fileformat" or int(tokens[1][4]) != 4):
            # TODO: Add messages based on incorrect VCF version or basic formatting error
            # TODO: add additional validation procedures
            print("{} format is invalid!".format(vcf_path.split('/')[-1]))
            raise InvalidVCFError(vcf_path, "{} format is invalid!".format(vcf_path.split('/')[-1]))

        vcf_version = tokens[1][4:7]
        print("Valid VCF file, version: {}".format(vcf_version))
        self.vcf_version = vcf_version
        #END validate_vcf

        return 

    def generate_vcf_stats(self, params, vcf_path):
        """
            :param commments go here
        """
        #BEGIN generate_vcf_stats
        
        if(not os.path.isdir(STORAGE_DIR)):
            print("\nCreating storage directory {}...".format(STORAGE_DIR))
            try:
                os.mkdir(STORAGE_DIR)
            except Exception as e:
                raise ValueError(e)

        print("Results will be written to {}".format(STORAGE_DIR))
        try:
            self.validate_vcf(vcf_path)
        except InvalidVCFError as ive:
            raise ValueError(ive.message)

        ## TODO: Use params to build linux command
        plink_cmd = ["plink"]
        plink_cmd.append('--vcf')
        plink_cmd.append(vcf_path)
        cmds = params['command_line_args'].split(';')
        for cmd in cmds:
            plink_cmd.append(cmd)
        plink_cmd.append('--freq')

        plink_cmd.append('--out')
        plink_cmd.append('frequencies')
        print("PLINK arguments: {}".format(plink_cmd))

        plink_output = []
        p = subprocess.Popen(plink_cmd, \
                            cwd = STORAGE_DIR, \
                            stdout = subprocess.PIPE, \
                            stderr = subprocess.STDOUT, \
                            shell = False)
        while True:
            line = p.stdout.readline()
            if not line: break
            #print(line.replace('\n', ''))
            tokens = line.split(':')
            if(tokens[0] == 'Error'):
                raise ValueError('PLINK 1.9 error: ' + line)
            elif(tokens[0] == 'Warning'):
                plink_output.append(line)
                print(line)
            elif(tokens[0] == 'Note'):
                plink_output.append(line)
                print(line)
        
        p.stdout.close()
        p.wait()

        if p.returncode != 0:
            raise ValueError("Error running PLINK, return code: " + str(p.returncode))
        # TODO: correct for the user supplied output file name.  Should I allow them to name?
        if not os.path.isfile(STORAGE_DIR+"frequencies.frq"):
            raise ValueError("PLINK failed to create frequency file {} frequencies.frq".format(STORAGE_DIR))
        
        return
        
        #END generate_vcf_stats
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
        print("Params passed to import_snp_data: {}".format(import_snp_params))
        try:
            self.generate_vcf_stats(import_snp_params, import_snp_params['input_file_path'])
        except Exception as e:
            raise ValueError(e)

        file_extensions = ['frq', 'log']
        indexHTML = "<head><body> "
        #for ext in file_extensions:
        indexHTML += """
                <a href="./frequencies.frq">Frequencies</a>
                <a href="./frequencies.log">Log</a>
            </body>
        </head>
        """
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
                                    'label': 'Save promoter_download.zip'
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
                      'vcf_version' : self.vcf_version
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
