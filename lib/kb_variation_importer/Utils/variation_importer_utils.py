import os
import subprocess
import uuid
import time
import shutil
import json
import numpy as np
import random
import string

from DataFileUtil.DataFileUtilClient import DataFileUtil
from KBaseReport.KBaseReportClient import KBaseReport
from GenomeAnnotationAPI.GenomeAnnotationAPIServiceClient import GenomeAnnotationAPI

class InvalidVCFException(Exception):
    def __init__(self, file_path, message):
        self.file_path = file_path
        self.message = message
    def __str__(self):
        return repr(self.message + self.file_path)

def log(message, prefix_newline=False):
    """Logging function, provides a hook to suppress or redirect log messages."""
    print(('\n' if prefix_newline else '') + '{0:.2f}'.format(time.time()) + ': ' + str(message))

#TODO: All manner of input validation checks.

class variation_importer_utils:

    def __init__(self, ctx, service_wiz_url, callback_url, scratch_file_path):
        self.scratch_file_path = scratch_file_path
        self.service_wiz_url = service_wiz_url
        self.dfu = DataFileUtil(callback_url)
        self.kbr = KBaseReport(callback_url, token=ctx['token'])

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

    def pretend_download_staging_file(self, vcf_file_path, scratch_path):
        print("Value passed to pretend_download: {}".format(vcf_file_path))
        shutil.copy('/kb/module/data/' + vcf_file_path, scratch_path + vcf_file_path)
        return { 'copy_file_path': os.path.join(scratch_path, vcf_file_path) }

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
        with open(vcf_file_path) as vcf:
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

    def validate_vcf(self, vcf_file_path, genome_ref):
        """
            :param vcf_file_path: string defining directory of VCF file
        """
        valid_vcf = False
        # TODO determine size of file.  May want to use HDF5 
        print("{} file size: {}".format(vcf_file_path, os.path.getsize(vcf_file_path)))
        print('\nValidating {}...'.format(vcf_file_path))
        try:
            f = open(vcf_file_path, "r")
        except Exception as e:
            print("Error opening file: {}".format(e))
            raise InvalidVCFException(vcf_file_path, e)
        
        line = f.readline()
        tokens = line.split('=')
        f.close()
        
        if(tokens[0] != "##fileformat" or int(tokens[1][4]) != 4):
            # TODO: Add messages based on incorrect VCF version or basic formatting error
            # TODO: add additional validation procedures
            print("{} version is invalid.".format(vcf_file_path.split('/')[-1]))
            vcf_version = -1
            # raise InvalidVCFException(vcf_file_path, "{} invalid version.".format(vcf_file_path.split('/')[-1]))
        else:
            vcf_version = tokens[1][4:7]

        vcf_contigs, vcf_genotypes = self._get_contigs_genotypes_vcf(vcf_file_path)
        # TODO: If no contigs, use file provided by user?  Or use file and if none found, use header?
        if not vcf_contigs:
            print("No contig data in {}!".format(vcf_file_path))
            # raise InvalidVCFException(vcf_file_path, "No contig data in VCF header!")
        population_data = self._create_fake_population(vcf_genotypes)
        #TODO: If version is below 4.1, attempt to convert using VCF tools?
        #TODO: Decide which validator to use.  
        if float(vcf_version) >= 4.1:
            print("Using vcf_validator_linux...")
            validator_cmd = ["vcf_validator_linux"] 
            validator_cmd.append("-i")
            validator_cmd.append(vcf_file_path)
            validator_cmd.append("-o")            
            # validation_output_dir = os.path.join(self.scratch_file_path, 'validation_report')
            # os.mkdir(validation_output_dir)
            validator_cmd.append(self.scratch_file_path)
        else:
            print("Using vcftools to validate...")
            validator_cmd = ["vcf-validator"]
            validator_cmd.append(vcf_file_path)
            print("VCF version below 4.1.  No validation logging.")

        print("Validator command: {}".format(validator_cmd))
        p = subprocess.Popen(validator_cmd, \
                            cwd = self.scratch_file_path, \
                            stdout = subprocess.PIPE, \
                            stderr = subprocess.STDOUT, \
                            shell = False)
        validator_output = []
        while True:
            line = p.stdout.readline()
            if not line: break
            # TODO: Output from vcf-validator is striclty warnings or errors.  
            # Warnings are readable, but errors are not (stack-trace).  I need to make them readable.  
            # Looks like the most readable error will begin with 'Vcf::throw...'
            # vcf_validator_linux is more readable by far, but it outputs exclusively to file!
            if line.startswith('[info] Summary'):
                tokens = line.split(':')
                print("Tokens: {}".format(tokens))


        p.wait()
       
        print("Return code from Popen {}".format(p.returncode))
        if p.returncode != 0:
            #TODO: Handle invalid VCF file type.  Output the file generated by validation?
            # It may be better to handle this outside of the function.  
            pass
        else:
            valid_vcf = True
        # p = subprocess.Popen(validator_cmd, \
        #                     cwd=self.scratch_file_path, \
        #                     stdout=subprocess.PIPE, \
        #                     stderr=subprocess.PIPE, \
        #                     shell=False)
        # out, err = p.communicate()
        # if p.returncode != 0:
        #     print("vcf_validator_linux return code: {}".format(p.returncode))
        #     # TODO: If I get a return code, parse and return the message to the user.  
        #     # TODO: Generate a report to return to the user?  
        path = tokens[1].split('/')
        validator_output_filename = path[-1]
        validator_output_filename.replace(' ', '')

        print("Validator output filename: {}".format(validator_output_filename))    
        print("Scratch file path in validate_vcf: {}".format(self.scratch_file_path))
        validator_output_path = os.path.join(self.scratch_file_path, validator_output_filename)
        # Manually removed newline character for some reason...
        validator_output_path = validator_output_path[:-1]
        if not os.path.isfile(validator_output_path):
            raise ValueError("Validator failed to create summary file!")

        # Retrieve Assembly object reference associated with genome.  
        try:
            assembly_ref = self._get_assembly_ref_from_genome(genome_ref)
        except Exception as e:
            print("Unable to retrieve {}".format(genome_ref))
            raise ValueError(e)

        # print("Assembly ref returned from genome: {}".format(assembly_ref))
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
            valid_vcf = False
        
        # TODO: Handle invalid contigs.  Number of invalid contigs could potentially be in the thousands. 
        # How can this be handled?  Write to file?

        returnVal = {
            "valid_vcf" : valid_vcf,
            "vcf_version" : vcf_version,
            "invalid_contigs" : invalid_contigs,
            "num_genotypes" : len(vcf_genotypes),
            "validation_output_path" : validator_output_path,
            "population" : population_data,
            "reference_genome" : genome_ref,
            "contigs" : vcf_contigs
        }

        return returnVal   

    def generate_vcf_stats(self, cmd_line_args, scratch_subdir_path, genome_ref):
        """
            :param commments go here
        """
        returnVal = {
            'plink_output' : []
        }
        ## TODO: Validate user supplied params and build PLINK command
        plink_cmd = ["plink"]
        plink_cmd.append('--vcf')
        plink_cmd.append(scratch_subdir_path)
        if(cmd_line_args is not None):
            cmds = cmd_line_args.split(';')
            for cmd in cmds:
                plink_cmd.append(cmd)
        # plink_cmd.append('--recode12')
        # plink_cmd.append('transpose')
        # plink_cmd.append('--output-missing-genotype')
        # plink_cmd.append("0")
        plink_cmd.append('--freq')
        plink_cmd.append('--out')
        plink_cmd.append('frequencies')

        print("PLINK arguments: {}".format(plink_cmd))

        plink_output = {
            "errors" : [],
            "warnings" : [],
            "notes" : []
        }
        p = subprocess.Popen(plink_cmd, \
                            cwd = self.scratch_file_path, \
                            stdout = subprocess.PIPE, \
                            stderr = subprocess.STDOUT, \
                            shell = False)
        while True:
            line = p.stdout.readline()
            if not line: break
            tokens = line.split(':')
            if(tokens[0] == 'Error'):
                plink_output['errors'].append(line)
                raise ValueError('PLINK 1.9 error: ' + line)
            elif(tokens[0] == 'Warning'):
                plink_output['warnings'].append(line)
                print(line)
            elif(tokens[0] == 'Note'):
                plink_output['notes'].append(line)
                print(line)
        
        p.stdout.close()
        p.wait()

        if p.returncode != 0:
            pass
            # TODO: Do something here?
            # raise ValueError("Error running PLINK, return code: " + str(p.returncode))


        # TODO: Check for whatever file type has been generated by user commands.
        if not os.path.isfile(self.scratch_file_path + "frequencies.frq"):
            raise ValueError("PLINK failed to create frequency file {} frequencies.frq".format(self.scratch_file_path))
        
        returnVal['plink_output'] = plink_output
        return returnVal

    def save_variation_to_ws(self, workspace_name, variation_results, vcf_staging_area_path, kinship_matrix):
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

    def generate_invalid_vcf_html_report(self, workspace_name, variation_results, scratch_path):
        template_dir = "/kb/module/lib/kb_variation_importer/Utils/invalid_report_template.html"
        try:
            with open(template_dir, 'r') as html, open(variation_results['validation_output_path'], 'r') as variation:
                template = html.read()
                variation_html = '<ul>'
                for line in variation.readlines():
                    variation_html += '<li>' + line + '</li>'
                variation_html += '</ul>'
                template = template.replace("###VCF_VERSION###", str(variation_results.get('vcf_version')))
                if(variation_results.get('invalid_contigs')):
                    contig_html = '<ul>'
                    for contig in variation_results.get('invalid_contigs'):
                        contig_html += '<li>' + contig + '</li>'
                    contig_html += '</ul>'
                    template = template.replace("###INVALID_CONTIG_TXT###", contig_html)
                template = template.replace("###VALIDATION_REPORT_TXT###", variation_html)
                template = template.replace("###GENOME_REF###", variation_results.get('reference_genome'))
        except Exception as e:
            print("Error generating invalid VCF Report.")
            raise ValueError(e)
        print("HTML Template: {}".format(template))

        report_file_path = os.path.join(scratch_path, 'html')
        os.mkdir(report_file_path)
        with open(os.path.join(report_file_path, 'index.html'), 'w') as output:
            output.write(template)

        try:
            html_upload_ret = self.dfu.file_to_shock({
                'file_path' : report_file_path,
                'make_handle' : 0,
                'pack' : 'zip'
            })
            print("invalid html report to shock ref: {}".format(html_upload_ret))
        except:
            raise ValueError('Error uploading HTML to shock')
        
        returnVal = {
            'shock_id' : html_upload_ret['shock_id'],
            'name' : os.path.basename(report_file_path),
            'label' : os.path.basename(report_file_path),
            'description' : 'Invalid_VCF_report'
        }
        return [returnVal]

    def build_report(self, workspace_name, html_links):
        report_params = {
                'objects_created': [], # Do I put the Variation object here?
                'message': '',
                'direct_html': None, # Is this where to include html displayed in narrative?
                'direct_html_index': 0,
                'file_links': [],
                'html_links': html_links,
                'html_window_height': 220,
                'workspace_name': workspace_name,
                'report_object_name': html_links[0].get('description') + str(uuid.uuid4())
        }
        kbr_output = self.kbr.create_extended_report(report_params)
        report_output = {
            'report_name' : kbr_output['name'],
            'report_ref' : kbr_output['ref']
        }

        return report_output
