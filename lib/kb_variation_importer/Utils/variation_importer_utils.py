import csv
import glob
import gzip
import json
import os
import random
import shutil
import string
import subprocess
import time
import uuid
import zipfile
from collections import Counter

import pandas as pd

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
# TODO: All manner of input validation checks.


class variation_importer_utils:

    def __init__(self, utility_params):
        self.params = utility_params
        # self.scratch = utility_params['scratch']
        self.scratch = os.path.join(
            utility_params['scratch'], 'variation_importer_'+str(uuid.uuid4()))
        os.mkdir(self.scratch)
        self.service_wiz_url = utility_params['srv-wiz-url']
        self.callback_url = utility_params['callback_url']

        self.dfu = DataFileUtil(self.callback_url)
        self.kbr = KBaseReport(self.callback_url, token=utility_params['token'])

    def _create_fake_location_data(self):
        location = {
            'lat': random.uniform(-90, 90),
            'lon': random.uniform(-180, 180),
            'elevation': random.uniform(0, 100),
            'description': "".join([random.choice(string.ascii_letters) for n in xrange(20)])
        }
        return location

    def _create_fake_straininfo(self, genotype_id):
        straininfo = {
            'source_id': genotype_id,
            'location_info': self._create_fake_location_data()
        }
        return straininfo

    def _create_fake_population(self, genotypes):
        population = {
            'description': 'Faker population data.',
            'strains': []
        }
        for genome in genotypes:
            population['strains'].append(self._create_fake_straininfo(genome))
        return population

    def _create_fake_kinship_matrix(self):
        kinship = {
            'row_ids': ['one', 'two'],
            'col_ids': ['one', 'two'],
            'kinship_coefficients': [
                [0.1, 0.1],
                [0.1, 0.1]
            ]
        }
        return kinship

    def _compare(self, s, t):
        return Counter(s) == Counter(t)

    def pretend_download_staging_file(self, vcf_filename, scratch):
        vcf_filepath = os.path.join(scratch, vcf_filename)
        shutil.copy('/kb/module/data/' + vcf_filename, vcf_filepath)
        return {'copy_file_path': vcf_filepath}

    def _generate_population(self, location_filepath, genotypes, population_description="None Provided"):
        locations = pd.read_csv(location_filepath, delimiter='\t')

        # Drop any missing data from id, latitude, or longitude.
        locations.dropna(subset=['id', 'latitude', 'longitude'], inplace=True)

        # Compare the location IDs with the genotype IDs
        if not(self._compare(locations.iloc[:, 0].astype(str).tolist(), genotypes)):
            log("Location IDs do not match Sample IDs in Variation file!")
            raise ValueError("Location IDs do not match Sample IDs in Variation file!")

        col_names = [x.lower() for x in locations.columns.values]
        expected_columns = ['id', 'latitude', 'longitude']
        optional_columns = ['elevation', 'description']

        # CHeck that first three columns match the expected columns.
        if not(self._compare(col_names[0:3], expected_columns)):
            raise ValueError("Missing or unexpected column names in {}".format(location_filepath))

        # If optional columns are not present, give default value for each.
        for col in optional_columns:
            if col not in col_names:
                if col == 'elevation':
                    locations[col] = 0.0
                else:
                    locations[col] = "None provided."

        population = {
            'description': population_description,
            'strains': []
        }
        for idx, row in locations.iterrows():
            population['strains'].append(
                {
                    'source_id': str(row['id']),
                    'location_info': {
                        'lat': row['latitude'],
                        'lon': row['longitude'],
                        'elevation': row['elevation'],
                        'description': row['description']
                    }
                }
            )

        return population

    def _validate_vcf(self, vcf_filepath, vcf_version):
        validation_output_dir = os.path.join(self.scratch, 'validation_' + str(uuid.uuid4()))
        os.mkdir(validation_output_dir)
        ## TODO: Make this choice more robust.  
        ## Attempt conversion to 4.1?
        if vcf_version >= 4.1:
            print("Using vcf_validator_linux...")
            validator_cmd = ["vcf_validator_linux"]
            validator_cmd.append("-i")
            validator_cmd.append(vcf_filepath)
            validator_cmd.append("-o")
            validator_cmd.append(validation_output_dir)
        else:
            print("Using vcftools to validate...")
            validator_cmd = ["vcf-validator"]
            validator_cmd.append(vcf_filepath)
            print("VCF version below 4.1.  No validation logging.")

        print("Validator command: {}".format(validator_cmd))
        p = subprocess.Popen(validator_cmd,
                             cwd=self.scratch,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT,
                             shell=False)
        validator_output = []
        while True:
            line = p.stdout.readline()
            if not line:
                break
            validator_output.append(line)

        p.wait()

        validation_output_filename = [f for f in os.listdir(
            validation_output_dir) if f.endswith('.txt')][0]
        validation_output_filepath = os.path.join(validation_output_dir, validation_output_filename)

        if not validation_output_filename:
            print('Validator did not generate log file!')
            raise Exception("Validator did not generate a log file.")

        log("Validator output filepath: {}".format(validation_output_filepath))

        log("Return code from validator {}".format(p.returncode))

        return validation_output_filepath, p.returncode

    # Retrieve contigs from assembly file.
    def _get_contigs_from_assembly(self, assembly_ref, type='Assembly'):
        try:
            assembly_data = self.dfu.get_objects(
                {'object_refs': [assembly_ref]}
            )['data'][0]['data']
        except Exception as e:
            print("Unable to retrieve Assembly reference: {}".format(assembly_ref))
            raise ValueError(e)
        raw_contigs = assembly_data['contigs']
        contigs = {}

        # Contigs returns just a dict with key and contig_id
        for key, value in raw_contigs.iteritems():
            contigs[str(key)] = value['contig_id']
        return raw_contigs

    def _get_version_contigs_genotypes(self, vcf_filepath):
        contigs = []
        genotypes = []
        version = ''
        with(gzip.open if vcf_filepath.endswith('.gz') else open)(vcf_filepath, 'rt') as vcf:
            line = vcf.readline()
            tokens = line.split('=')

            if not(tokens[0].startswith('##fileformat')):
                log("Invalid VCF.  ##fileformat line in meta is improperly formatted.")
                raise ValueError("Invalid VCF.  ##fileformat line in meta is improperly formatted.")
            version = float(tokens[1][-4:].rstrip())
            log("VCF version: {}".format(version))
            for line in vcf:
                if line.startswith("#CHROM"):
                    log("#CHROM encountered, exiting loop.")
                    genotypes = line.split()[9:]
                    log("Number Genotypes in vcf: {}".format(len(genotypes)))
                    break
                tokens = line.split("=")

                if tokens[0].startswith('##contig'):
                    contigs.append(tokens[2][:-2])
        return version, contigs, genotypes

    # Arabidopsis ref: 18590/2/8
    def _get_assembly_ref_from_genome(self, genome_ref):
        ga = GenomeAnnotationAPI(self.service_wiz_url)
        inputs_get_assembly = {'ref': genome_ref}
        try:
            assembly_object_ref = ga.get_assembly(inputs_get_assembly)
        except Exception as e:
            print("Unable to retrieve Assembly reference ID from Genome ref_id: {}".format(genome_ref))
            raise Exception(e)

        return assembly_object_ref

    def _generate_output_file_list(self):
        log('Start packing result files')
        output_files = list()

        result_file = os.path.join(self.scratch, 'variation_importer_results.zip')
        excluded_extensions = ['.zip', '.vcf', '.vcf.gz', '.html', '.DS_Store']
        with zipfile.ZipFile(result_file, 'w',
                             zipfile.ZIP_DEFLATED,
                             allowZip64=True) as zip_file:
            for root, dirs, files in os.walk(self.scratch):
                for file in files:
                    if not (file.endswith(tuple(excluded_extensions))
                            # file.endswith('.zip') or
                            # file.endswith('.vcf') or
                            # file.endswith('.vcf.gz') or
                            # file.endswith('.html') or
                            # file.endswith('.DS_Store')
                            ):
                        zip_file.write(os.path.join(root, file), file)

        output_files.append({'path': result_file,
                             'name': os.path.basename(result_file),
                             'label': os.path.basename(result_file),
                             'description': 'File(s) generated by Variation Importer'})
        log("Importer output generated: {}".format(output_files))

        return output_files

    def _generate_report(self, params, variation_results,
                         variation_file_path):

        stats_results = self._generate_variation_stats(None,
                                                       variation_file_path)

        html_report = self._generate_html_report(variation_results, stats_results)

        file_links = self._generate_output_file_list()
        objects = []
        if (variation_results['valid_variation_file']):
            objects = [{
                'ref': variation_results['variation_obj_ref'],
                'description': 'Variation Object created by VCF Importer'
            }]
        
        report_params = {
            'objects_created': objects,
            'message': '',
            'direct_html_link_index': 0,
            'file_links': file_links,
            'html_links': html_report,
            'html_window_height': 330,
            'workspace_name': params['workspace_name'],
            'report_object_name': 'variation_importer_report_' + str(uuid.uuid4())
        }
        kbr_output = self.kbr.create_extended_report(report_params)
        report_output = {
            'report_name': kbr_output['name'],
            'report_ref': kbr_output['ref'],
            'variation_ref': variation_results['variation_obj_ref']
        }
        log("Returning from _generate_report!")
        return report_output

    def _generate_html_report(self, variation_results, stats_output=None):
        """
            _generate_html_report: generate html report from output files
        """
        html_report = list()
        print("Validation output filepath passed to html report: {}".format(
            variation_results['validation_output_filepath']))
        try:
            report_dir = os.path.join(self.scratch, 'html')
            os.mkdir(report_dir)

            with open(template_dir, 'r') as html, open(variation_results['validation_output_filepath'], 'r') as validation:

                validation_content = '<p><h4>{} '.format(variation_results['variation_filename'])
                if variation_results.get('valid_variation_file'):
                    validation_content += '<em><i>is</i> a valid </em> variation file.'
                else:
                    validation_content += '<em><i>is not</i> a valid </em>variation file. Details below.'
                validation_content += '</h4></p>'

                report = html.read()

                # Discard the first line of the validation file.  It is irrelevant.
                validation.readline()

                validation_content += '<p><h4>Errors and warning generated by VCF validator:</h4></p>'
                validation_content += '<ul>'
                for line in validation.readlines():
                    validation_content += '<li>{}</li>'.format(line)
                validation_content += '</ul>'

                if variation_results.get('invalid_contigs'):
                    validation_content += '<h4>The following Contigs were not found in the reference genome.  The possible contigs have been written to the file {}.  Please see the associated links to download.</h4>'.format(
                        variation_results.get('genome_ref'), 'valid_contigs.txt')
                    validation_content += '<ul>'
                    for contig in variation_results.get('invalid_contigs'):
                        validation_content += '<li>{}</li>'.format(contig)
                    validation_content += '</ul>'

                # if not variation_results.get('contigs'):
                #     validation_content += '<h4>No contig information was included in the VCF file header!  Please recreate the VCF file with each contig described in the meta description </h4>'
                report = report.replace('Validation_Results', validation_content)

                if(stats_output.get('stats_file_dir')):
                    summary_results = '<p><h4>Summary Statistics</h4></p>'
                    summary_results += '''
                                        <table>
                                            <tr>
                                                <th>Number of SNPs</th>
                                                <th>Number of Genotypes </th>
                                            </tr>
                                        '''
                    summary_results += '<tr>'
                    summary_results += '<td>{}</td><td>{}</td>'.format(
                        'To be added later', variation_results['num_genotypes'])
                    summary_results += '</tr></table>'
                    report = report.replace('Variation_Statistics', summary_results)

                # visualization
                image_content = ''
                if(stats_output.get('stats_img_dir')):
                    image_dir = stats_output.get('stats_img_dir')

                    for file in glob.glob(os.path.join(image_dir, '*.png')):
                        shutil.move(file, report_dir)

                    for image in glob.glob(report_dir + "/*.png"):
                        image = image.replace(report_dir + '/', '')
                        caption = image.replace(report_dir + '/', '').replace('.png', '')
                        image_content += '<p style="text-align:center"><img align="center" src="{}" ' \
                            '></a><a target="_blank"><br>' \
                            '<p align="center">{}</p></p>'.format(image, caption)

                else:
                    image_content += 'No visualizations generated.'

                report = report.replace("Visualization_Results", image_content)
        except Exception as e:
            print("Error generating HTML report.")
            raise

        report_file_path = os.path.join(report_dir, 'index.html')
        with open(report_file_path, 'w') as output:
            output.write(report)
        try:
            html_upload_ret = self.dfu.file_to_shock({
                'file_path': report_file_path,
                'make_handle': 0,
                'pack': 'zip'
            })
            log("Variation HTML report to shock ref: {}".format(html_upload_ret))
        except:
            raise ValueError('Error uploading HTML to shock')

        html_report.append({
            'shock_id': html_upload_ret['shock_id'],
            'name': os.path.basename(report_file_path),
            'label': os.path.basename(report_file_path),
            'description': 'HTML report for Variation Importer'
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

        # TODO: Validate user supplied params and build PLINK command
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
        plink_cmd.append('--hardy')
        # plink_cmd.append('gz')
        plink_cmd.append('--out')
        plink_cmd.append(variation_filepath)

        print("PLINK arguments: {}".format(plink_cmd))

        plink_output = {
            "errors": [],
            "warnings": []
            # "notes" : []
        }
        p = subprocess.Popen(plink_cmd,
                             cwd=file_output_directory,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT,
                             shell=False)
        while True:
            line = p.stdout.readline()
            if not line:
                break
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
        plink_output_filepath = os.path.join(file_output_directory, 'plink_cli_output.txt')
        with open(plink_output_filepath, 'w') as plink:
            for data in plink_output:
                plink.write("{}: {}\n".format(data, plink_output[data]))

        plink_output_files = [f for f in os.listdir(self.scratch) if f.startswith(
            os.path.basename(variation_filepath) + '.')]

        for file in plink_output_files:
            shutil.move(os.path.join(self.scratch, file), file_output_directory)

        if p.returncode != 0:
            log("PLINK encountered an error during runtime.  Please see log file.")

        variation_filename = os.path.basename(variation_filepath)
        base_filepath = os.path.join(file_output_directory, variation_filename)
        freq_filepath = base_filepath + '.frq'

        maf_script_filepath = '/kb/module/lib/kb_variation_importer/Utils/MAF_check.R'
        hwe_script_filepath = '/kb/module/lib/kb_variation_importer/Utils/HWE.R'
        log("Frequency filepath: {}".format(freq_filepath))
        # TODO: make function to do Rscript calls.
        # generate visualizations and store in directory
        maf_command = ['Rscript']
        maf_command.append('--no-save')
        maf_command.append('--vanilla')
        maf_command.append(maf_script_filepath)
        maf_command.append(freq_filepath)
        maf_command.append("Minor Allele Frequencies.png")
        print("MAF command: {}".format(maf_command))
        r = subprocess.Popen(maf_command,
                             cwd=image_output_directory,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT,
                             shell=False)
        r.wait()

        if r.returncode != 0:
            log("Error creating MAF histogram in R")

        hwe_filepath = base_filepath + '.hwe'
        zoom_filepath = hwe_filepath + '.zoom'
        log("HWE filepath: {}".format(hwe_filepath))
        zoom_command = '''awk '{{ if ($9 < 0.00001) print $0 }}' {} > {}'''.format(
            hwe_filepath, zoom_filepath)
        log("Zoom cmd: {}".format(zoom_command))

        try:
            z = subprocess.Popen(zoom_command,
                                 cwd=file_output_directory,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT,
                                 shell=True)
            z.wait()

            if z.returncode != 0:
                log("Error creating HWE zoom file.")

        except Exception as e:
            log("Error creating zoom HWE file: {}".format(e))

        hwe_command = ['Rscript']
        hwe_command.append('--no-save')
        hwe_command.append('--vanilla')
        hwe_command.append(hwe_script_filepath)
        hwe_command.append(hwe_filepath)
        hwe_command.append("Hardy-Weinberg Equilibrium.png")
        hwe_command.append(zoom_filepath)
        hwe_command.append("Hardy-Weinberg Equilibrium Zoom.png")
        print("MAF command: {}".format(hwe_command))
        h = subprocess.Popen(hwe_command,
                             cwd=image_output_directory,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT,
                             shell=False)
        h.wait()

        if h.returncode != 0:
            log("Error generating HWE Zoom plot")

        return {'stats_file_dir': file_output_directory,
                'stats_img_dir': image_output_directory}

    def _save_variation_to_ws(self, workspace_name, variation_object_name, variation_obj, variation_filepath, kinship_matrix):
        ws_id = self.dfu.ws_name_to_id(workspace_name)
        try:
            vcf_shock_return = self.dfu.file_to_shock({
                'file_path': variation_filepath,
                'make_handle': 1,
                'pack': 'gzip'})
        except Exception as e:
            print("Error uploading file to shock!")
            raise ValueError(e)

        variation_obj['variation_file_reference'] = vcf_shock_return.get('shock_id')

        info = self.dfu.save_objects(
            {
                'id': ws_id,
                'objects': [{
                    'type': 'KBaseGwasData.Variations',
                    'data': variation_obj,
                    'name': variation_object_name
                }]
            })[0]

        variation_ref = "%s/%s/%s" % (info[6], info[0], info[4])
        log("Variation reference created: {}".format(variation_ref))
        return variation_ref

    def validate_vcf(self, params):
        """
            :param params: dict containing all input parameters.
        """

        returnVal = {}
        valid_vcf_file = True

        try:
            vcf_filepath = self.pretend_download_staging_file(
                params['variation_file_subdir_path'], self.scratch).get('copy_file_path')

            location_filepath = self.pretend_download_staging_file(
                params['variation_attributes_subdir_path'], self.scratch).get('copy_file_path')

        except Exception as e:
            raise Exception("Unable to download {} from staging area.".format(
                params['variation_file_subdir_path']))

        try:
            location_filepath = self.pretend_download_staging_file(
                params['variation_attributes_subdir_path'], self.scratch).get('copy_file_path')

        except Exception as e:
            raise Exception("Unable to download {} from staging area.".format(
                params['variation_attributes_subdir_path']))

        # Check file size
        log("{} file size: {}".format(vcf_filepath, os.path.getsize(vcf_filepath)))
        log('\nValidating {}...'.format(vcf_filepath))

        vcf_version, vcf_contigs, vcf_genotypes = self._get_version_contigs_genotypes(vcf_filepath)

        if not vcf_contigs:
            log("No contig data in {} header.".format(vcf_filepath))
            raise ValueError("No contig data in {} header.".format(vcf_filepath))

        if (vcf_version < 4.1):
            log("VCF file is version {}.  Must be at least version 4.1".format(vcf_version))
            raise ValueError(
                "VCF file is version {}.  Must be at least version 4.1".format(vcf_version))

        # Generate population object
        population = self._generate_population(location_filepath, vcf_genotypes)

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

        log("Length of assembly contigs: {}".format(len(assembly_contigs)))
        # Compare contig IDs from VCF to those in the Assembly object
        invalid_contigs = []
        for contig in vcf_contigs:
            if contig not in assembly_contigs.keys():
                invalid_contigs.append(contig)

        if invalid_contigs:
            log("Invalid contig IDs found in {}".format(vcf_filepath))
            valid_contig_filepath = os.path.join(
                self.scratch, 'valid_contigs.txt')
            log("Writing valid contigs to file: {}".format(valid_contig_filepath))
            with open(valid_contig_filepath, 'w') as icf:
                for contig in assembly_contigs:
                    icf.write(contig + '\n')
            valid_vcf_file = False

        validation_output_filepath, returncode = self._validate_vcf(vcf_filepath, vcf_version)

        if returncode != 0:
            valid_vcf_file = False

        kinship_matrix = self._create_fake_kinship_matrix()

        variation_obj_ref = ''
        if valid_vcf_file:
            variation_object = {
                "genome": params['genome_ref'],
                "population": population,
                "contigs": vcf_contigs,
                "comment": "Comments go here",
                "assay": "Assay data goes gere.",
                "originator": "PI/Lab info goes here",
                "pubmed_id": "PubMed ID goes here",
                "kinship_info": kinship_matrix
            }

            variation_obj_ref = self._save_variation_to_ws(params['workspace_name'],
                                                           params['variation_object_name'],
                                                           variation_object,
                                                           vcf_filepath,
                                                           kinship_matrix)

        log("Variation object reference: {}".format(variation_obj_ref))
        variation_report_metadata = {
            'valid_variation_file': valid_vcf_file,
            'variation_obj_ref': variation_obj_ref,
            'variation_filename': os.path.basename(vcf_filepath),
            'validation_output_filepath': validation_output_filepath,
            'vcf_version': vcf_version,
            'num_genotypes': len(vcf_genotypes),
            'num_contigs': len(vcf_contigs),
            'invalid_contigs': invalid_contigs
        }

        returnVal = self._generate_report(params, variation_report_metadata, vcf_filepath)

        return returnVal
