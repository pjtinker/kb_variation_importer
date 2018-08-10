import gzip
import glob
import json
import os
import random
import shutil
import string
import subprocess
import time
import uuid
import zipfile

from DataFileUtil.DataFileUtilClient import DataFileUtil
from GenomeAnnotationAPI.GenomeAnnotationAPIServiceClient import \
    GenomeAnnotationAPI
from KBaseReport.KBaseReportClient import KBaseReport


class InvalidVCFException(Exception):
    def __init__(self, file_path, message):
        self.file_path = file_path
        self.message = message
    def __str__(self):
        return repr(self.message + self.file_path)

def log(message, prefix_newline=False):
    """Logging function, provides a hook to suppress or redirect log messages."""
    print(('\n' if prefix_newline else '') + '{0:.2f}'.format(time.time()) + ': ' + str(message))

template_dir = "/kb/module/lib/kb_variation_importer/Utils/invalid_report_template.html"
#TODO: All manner of input validation checks.

class variation_importer_utils:

    def __init__(self, utility_params):
        self.params = utility_params
        # self.scratch = utility_params['scratch']
        self.scratch = os.path.join(utility_params['scratch'], 'variation_importer_'+str(uuid.uuid4()))
        os.mkdir(self.scratch)
        self.service_wiz_url = utility_params['srv-wiz-url']
        self.callback_url = utility_params['callback_url']

        self.dfu = DataFileUtil(self.callback_url)
        self.kbr = KBaseReport(self.callback_url, token=utility_params['token'])

    def _create_fake_location_data(self):
        location = {
            'lat' : random.uniform(-90, 90),
            'lon' : random.uniform(-180, 180),
            'elevation' : random.uniform(0, 100),
            'description' : "".join([random.choice(string.ascii_letters) for n in xrange(20)])
        }
        return location

    def _create_fake_straininfo(self, genotype_id):
        straininfo = {
            'source_id' : genotype_id,
            'location_info' : self._create_fake_location_data()
        }
        return straininfo

    def _create_fake_population(self, genotypes):
        population = {
            'description' : 'Faker population data.',
            'strains' : []
        }
        for genome in genotypes:
            population['strains'].append(self._create_fake_straininfo(genome))
        return population

    def create_fake_kinship_matrix(self):
        kinship = {
            'row_ids' : ['one', 'two'],
            'col_ids' : ['one', 'two'],
            'kinship_coefficients' : [
                [0.1, 0.1],
                [0.1, 0.1]
            ]
        }
        return kinship

    def pretend_download_staging_file(self, vcf_filename, scratch):
        print("Value passed to pretend_download: {}".format(vcf_filename))
        vcf_filepath = os.path.join(scratch, vcf_filename)
        shutil.copy('/kb/module/data/' + vcf_filename, vcf_filepath)
        print("Copy file path: {}".format(vcf_filepath))
        return { 'copy_file_path': vcf_filepath }

    # Retrieve contigs from assembly file.  
    def _get_contigs_from_assembly(self, assembly_ref, type='Assembly'):
        try:
            assembly_data = self.dfu.get_objects(
                {'object_refs': [assembly_ref]}
            )['data'][0]['data']
        except Exception as e:
            print("Unable to retrieve Assembly reference: {}".format(assembly_ref))
        raw_contigs = assembly_data['contigs']
        contigs = {}

        # Contigs returns just a dict with key and contig_id
        for key, value in raw_contigs.iteritems():
            contigs[str(key)] = value['contig_id']
        return raw_contigs
        
    def _get_contigs_genotypes_vcf(self, vcf_file_path):
        contigs = []
        genotypes = []
        # with open(vcf_file_path) as vcf:
        with(gzip.open if vcf_file_path.endswith('.gz') else open)(vcf_file_path, 'rt') as vcf:
            for line in vcf:
                if line.startswith("#CHROM"):
                    print("#CHROM encountered, exiting loop.")
                    genotypes = line.split()[9:]
                    print("Number Genotypes in vcf: {}".format(len(genotypes)))
                    break
                tokens = line.split("=")
                if tokens[0] == "##contig":
                    contigs.append(tokens[2][:-2])
        return contigs, genotypes

    # Arabidopsis ref: 18590/2/8
    def _get_assembly_ref_from_genome(self, genome_ref):
        ga = GenomeAnnotationAPI(self.service_wiz_url)
        inputs_get_assembly = {'ref': genome_ref }
        try:
            assembly_object_ref = ga.get_assembly(inputs_get_assembly)
        except Exception as e:
            print("Unable to retrieve Assembly reference ID from Genome ref_id: {}".format(genome_ref))
            raise Exception(e)
        
        return assembly_object_ref

    def _generate_output_file_list(self):
        log('Start packing result files')
        output_files = list()

        # output_directory = os.path.join(self.scratch, 'importer_output')
        # os.mkdir(output_directory)
        result_file = os.path.join(self.scratch, 'variation_importer_results.zip')


        with zipfile.ZipFile(result_file, 'w',
                             zipfile.ZIP_DEFLATED,
                             allowZip64=True) as zip_file:
            for root, dirs, files in os.walk(self.scratch):
                for file in files:
                    if not (file.endswith('.zip') or
                            file.endswith('.vcf') or
                            file.endswith('.vcf.gz') or
                            file.endswith('.html') or
                            file.endswith('.DS_Store')):
                        zip_file.write(os.path.join(root, file), file)


        output_files.append({'path': result_file,
                             'name': os.path.basename(result_file),
                             'label': os.path.basename(result_file),
                             'description': 'File(s) generated by Variation Importer'})
        log("Importer output generated: {}".format(output_files))
        
        return output_files


    def _generate_html_report(self, variation_results, stats_output=None):
        """
            _generate_html_report: generate html report from output files
        """
        html_report = list()
        try:
            with open(template_dir, 'r') as html, open(variation_results['validation_output_filepath'], 'r') as validation:
                
                validation_content = '<p><h4>{} '.format(variation_results['variation_filename'])
                if variation_results.get('valid_variation_file'):
                    validation_content += 'is a valid variation file.</h4>'
                else:
                    validation_content += 'is not a valid variation file. Details below.</h4>'
                validation_content += '</p>'

                report = html.read()

                # Discard the first line of the validation file.  It is irrelevant.
                validation.readline()     
                validation_content += '<p>Errors and warning generated by VCF validator:</p>'
                validation_content += '<ul>'    
                for line in validation.readlines():
                    validation_content += '<li>{}</li>'.format(line)
                validation_content += '</ul>'
                
                if variation_results.get('invalid_contigs'):
                    validation_content += '<h4>The following Contigs were not found in Genome {}.  The possible contigs have been written to the file valid_contigs.txt.  Please see the associated links to download.</h4>'.format(variation_results.get('genome_ref'))
                    validation_content += '<ul>'
                    for contig in variation_results.get('invalid_contigs'):
                        validation_content += '<li>{}</li>'.format(contig)
                    validation_content += '</ul>'

                    invalid_contig_filepath = os.path.join(self.scratch, 'valid_contigs.txt') 
                    with open(invalid_contig_filepath, 'w') as icf:
                        for contig in validation_content.get('contigs'):
                            icf.write(contig)

                if not variation_results.get('contigs'):
                    validation_content += '<h4>No contig information was included in the VCF file header!  Please recreate the VCF file with each contig described in the meta description </h4>'
                report = report.replace('Validation_Results', validation_content)

                # report = report.replace("GENOME_REF", variation_results.get('reference_genome'))
        except Exception as e:
            print("Error generating HTML report.")
            raise ValueError(e)

        report_file_path = os.path.join(self.scratch, 'html')
        os.mkdir(report_file_path)
        with open(os.path.join(report_file_path, 'index.html'), 'w') as output:
            output.write(report)

        try:
            html_upload_ret = self.dfu.file_to_shock({
                'file_path' : report_file_path,
                'make_handle' : 0,
                'pack' : 'zip'
            })
            log("Variation HTML report to shock ref: {}".format(html_upload_ret))
        except:
            raise ValueError('Error uploading HTML to shock')
        
        html_report.append({
            'shock_id' : html_upload_ret['shock_id'],
            'name' : os.path.basename(report_file_path),
            'label' : os.path.basename(report_file_path),
            'description' : 'HTML report for Variation Importer'
        })

        return html_report

  

    def _generate_variation_stats(self, cmd_line_args, variation_filepath):
        """
            :param commments go here
        """
        file_output_directory = os.path.join(self.scratch, 'stats_' + str(uuid.uuid4()))
        os.mkdir(file_output_directory)

        image_output_directory = os.path.join(self.scratch, 'stats_images_' + str(uuid.uuid4())) 
        os.mkdir(image_output_directory)
        ## TODO: Validate user supplied params and build PLINK command
        plink_cmd = ["plink"]
        plink_cmd.append('--vcf')
        plink_cmd.append(variation_filepath)
        if(cmd_line_args is not None):
            cmds = cmd_line_args.split(';')
            for cmd in cmds:
                plink_cmd.append(cmd)
        # plink_cmd.append('--recode12')
        # plink_cmd.append('transpose')
        # plink_cmd.append('--output-missing-genotype')
        # plink_cmd.append("0")
        plink_cmd.append('--freq')
        # plink_cmd.append('gz')
        plink_cmd.append('--out')
        plink_cmd.append(variation_filepath)

        print("PLINK arguments: {}".format(plink_cmd))

        plink_output = {
            "errors" : [],
            "warnings" : []
            # "notes" : []
        }
        p = subprocess.Popen(plink_cmd, \
                            cwd = file_output_directory, \
                            stdout = subprocess.PIPE, \
                            stderr = subprocess.STDOUT, \
                            shell = False)
        while True:
            line = p.stdout.readline()
            if not line: break
            # log(line)
            tokens = line.split(':')
            if(tokens[0] == 'Error'):
                plink_output['errors'].append(line)
                raise ValueError('PLINK 1.9 error: ' + line)
            elif(tokens[0] == 'Warning'):
                plink_output['warnings'].append(line)
                print(line)
            # elif(tokens[0] == 'Note'):
            #     plink_output['notes'].append(line)
            #     print(line)
        
        p.stdout.close()
        p.wait()
        plink_output_filepath = os.path.join(file_output_directory, 'plink_commandline_output.txt')
        with open(plink_output_filepath, 'w') as plink:
            for data in plink_output:
                plink.write("{}: {}".format(data, plink_output[data]))

        plink_output_files = [f for f in os.listdir(self.scratch) if f.startswith(os.path.basename(variation_filepath) + '.')]

        for file in plink_output_files:
            shutil.move(os.path.join(self.scratch, file), file_output_directory)

        if p.returncode != 0:
            # TODO: Handle PLINK errors.
            log("PLINK encountered an error during runtime.  Please see log file.")    
            

        # TODO: generate visualizations and store in directory
        
        # TODO: Check for whatever file type has been generated by user commands.
        # if not os.path.isfile():
        #     raise ValueError("PLINK failed to create frequency file frequencies.frq")

        return { 'stats_file_dir' : file_output_directory, \
                 'stats_img_dir' : image_output_directory }


    def _save_variation_to_ws(self, workspace_name, variation_results, vcf_staging_area_path, kinship_matrix):
        ws_id = self.dfu.ws_name_to_id(workspace_name)
        print("workspace id: {}".format(ws_id))
        ws_id = '18590' # TODO: remove hard-coded test case!
        try:
            # TODO: Change this from static dir to staging area upload?
            vcf_shock_return = self.dfu.file_to_shock({
                                        'file_path': vcf_staging_area_path,
                                        'make_handle': 1,
                                        'pack' : 'gzip'})
        except Exception as e:
            print("Error uploading file to shock!")
            raise ValueError(e)

        variation = {
            'population' : variation_results.get('population'),
            'comment' : 'Variation comment here',
            'assay' : 'Assay data here',
            'originator' : 'PI / LAB',
            'genome' : variation_results.get('genome_ref'),
            'pubmed_id' : 'Pub Med ID',
            'contigs' : variation_results.get('contigs'),
            'variation_file_reference' : vcf_shock_return.get('shock_id'),
            'kinship_info' : kinship_matrix
        }
        info = self.dfu.save_objects(
            {
                'id' : ws_id,
                'objects': [{
                    'type': 'KBaseGwasData.Variations',
                    'data': variation,
                    'name': workspace_name
                }]
            })[0]

        variation_ref = "%s/%s/%s" % (info[6], info[0], info[4])
        print("Variation reference created: {}".format(variation_ref))
        return variation_ref



    def _generate_report(self, params, variation_results, \
                         variation_file_path, validation_output_dir = None ):

        stats_results = self._generate_variation_stats(params['command_line_args'], \
                        variation_file_path)

        html_report = self._generate_html_report(variation_results, stats_results)

        file_links = self._generate_output_file_list()

        report_params = {
                'objects_created': [], # Do I put the Variation object here?
                'message': '',
                'direct_html': None, # Is this where to include html displayed in narrative?
                'direct_html_index': 0,
                'file_links': file_links,
                'html_links': html_report,
                'html_window_height': 220,
                'workspace_name': params['workspace_name'],
                'report_object_name': 'variation_importer_report_' + str(uuid.uuid4())
        }
        kbr_output = self.kbr.create_extended_report(report_params)
        report_output = {
            'report_name' : kbr_output['name'],
            'report_ref' : kbr_output['ref']
        }

        return report_output

    def validate_vcf(self, params):
        """
            :param params: dict containing all input parameters.
        """
        returnVal = {}

        variation_results = {
            "valid_variation_file" : True,
            "genome_ref" : params['genome_ref'],
            "invalid_contigs" : None
        }

        try:
            vcf_file_path = self.pretend_download_staging_file(
                params['staging_file_subdir_path'], self.scratch).get('copy_file_path')
        except Exception as e:
            raise Exception("Unable to download {} from staging area.".format(params['staging_file_subdir_path']))
        
        variation_results['variation_filename'] = os.path.basename(vcf_file_path)
        # Check file size 
        print("{} file size: {}".format(vcf_file_path, os.path.getsize(vcf_file_path)))
        print('\nValidating {}...'.format(vcf_file_path))

        try:
            with (gzip.open if vcf_file_path.endswith('.gz') else open)(vcf_file_path, 'rt') as f:
                line = f.readline()
                tokens = line.split('=')
        except Exception as e:
            print("Error opening file: {}".format(e))
            raise InvalidVCFException(vcf_file_path, e)

        # TODO: Handle the version string more thoroughly.  Make it less brittle.
        if(tokens[0] != "##fileformat" or tokens[1][4] != '4'):
            # TODO: Add messages based on incorrect VCF version or basic formatting error
            # TODO: add additional validation procedures
            print("{} version is invalid.".format(vcf_file_path.split('/')[-1]))
            raise InvalidVCFException(vcf_file_path, "{} version is not 4.1 or higher!".format(vcf_file_path.split('/')[-1]))
        else:
            vcf_version = tokens[1][4:7]

        variation_results["vcf_version"] = vcf_version

        #TODO: If version is below 4.1, attempt to convert using VCF tools?
        #TODO: Decide which validator to use.  

        validation_output_dir = os.path.join(self.scratch, 'validation_' + str(uuid.uuid4()))
        os.mkdir(validation_output_dir)

        if float(vcf_version) >= 4.1:
            print("Using vcf_validator_linux...")
            validator_cmd = ["vcf_validator_linux"] 
            validator_cmd.append("-i")
            validator_cmd.append(vcf_file_path)
            validator_cmd.append("-o")            
            validator_cmd.append(validation_output_dir)
        else:
            print("Using vcftools to validate...")
            validator_cmd = ["vcf-validator"]
            validator_cmd.append(vcf_file_path)
            print("VCF version below 4.1.  No validation logging.")

        print("Validator command: {}".format(validator_cmd))
        p = subprocess.Popen(validator_cmd, \
                            cwd = self.scratch, \
                            stdout = subprocess.PIPE, \
                            stderr = subprocess.STDOUT, \
                            shell = False)
        validator_output = []
        while True:
            line = p.stdout.readline()
            if not line: break
            validator_output.append(line)

        p.wait()

        validation_output_filename = [f for f in os.listdir(validation_output_dir) if f.endswith('.txt')][0]
        validation_output_filepath = os.path.join(validation_output_dir, validation_output_filename)
        variation_results['validation_output_filepath'] = validation_output_filepath
        if not validation_output_filename:
            print('Validator did not generate log file!')
            raise Exception("Validator did not generate a log file.")

        print("Validator output filepath: {}".format(validation_output_filepath))    

        print("Return code from Popen {}".format(p.returncode))
        if p.returncode != 0:
            variation_results['valid_variation_file'] = False
            returnVal = self._generate_report(params['workspace_name'], variation_results)
            return returnVal

        vcf_contigs, vcf_genotypes = self._get_contigs_genotypes_vcf(vcf_file_path)
        # TODO: If no contigs, use file provided by user?  Or use file and if none found, use header?
        if not vcf_contigs:
            print("No contig data in {}!".format(vcf_file_path))
            variation_results['valid_variation_file'] = False        

        variation_results['contigs'] = vcf_contigs
        variation_results['num_genotypes'] = len(vcf_genotypes)

        # Retrieve Assembly object reference associated with genome.  
        try:
            assembly_ref = self._get_assembly_ref_from_genome(params['genome_ref'])
        except Exception as e:
            print("Unable to retrieve {}".format(params['genome_ref']))
            raise ValueError(e)

        # Retrieve contig list from Assembly object.
        try:
            assembly_contigs = self._get_contigs_from_assembly(assembly_ref)
        except Exception as e:
            print("Unable to retrieve contigs from Assembly ref: {}".format(assembly_ref))
            raise ValueError(e)

        # Compare contig IDs from VCF to those in the Assembly object
        invalid_contigs = []
        for contig in vcf_contigs:
            if contig not in assembly_contigs.keys():
                invalid_contigs.append(contig)
        
        if invalid_contigs:
            print("Invalid contig IDs found in {}".format(vcf_file_path))
            variation_results['invalid_contigs'] = invalid_contigs
            variation_results['valid_variation_file'] = False

        variation_results['population'] = self._create_fake_population(vcf_genotypes)

        # validation_files = self._generate_output_file_list()
        # stats_results = self._generate_variation_stats(params['command_line_args'], vcf_file_path)
        # validation_files.append(stats_results.get('stats_files'))
        # print("Validation files after append: {}".format(validation_files))

        returnVal = self._generate_report(params, variation_results, vcf_file_path,validation_output_dir)
        
        return returnVal 
