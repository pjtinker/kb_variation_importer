import os
import subprocess
import uuid
import time
import shutil
import json

from DataFileUtil.DataFileUtilClient import DataFileUtil
from KBaseReport.KBaseReportClient import KBaseReport
from GenomeAnnotationAPI.GenomeAnnotationAPIServiceClient import GenomeAnnotationAPI

class InvalidVCFError(Exception):
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

    def __init__(self, service_wiz_url, dfu, storage_dir):
        self.storage_dir = storage_dir
        self.service_wiz_url = service_wiz_url
        self.dfu = dfu

    def pretend_download_staging_file(self, vcf_file_path, scratch_path):
        print("Value passed to pretend_download: {}".format(vcf_file_path))
        shutil.copy('/kb/module/data/' + vcf_file_path, scratch_path + vcf_file_path)
        return { 'copy_file_path': scratch_path + vcf_file_path }

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

        for key, value in raw_contigs.iteritems():
            contigs[str(key)] = 0
        return contigs
        
    def _get_contigs_from_vcf(self, vcf_file_path):
        contigs = []
        with open(vcf_file_path) as vcf:
            for line in vcf:
                if line.startswith("#CHROM"):
                    print("#CHROM encountered, exiting loop.")
                    break
                tokens = line.split("=")
                if tokens[0] == "##contig":
                    contigs.append(tokens[2][:-2])
        return contigs
    # Arabidopsis ref: 18590/2/8
    def _get_assembly_ref_from_genome(self, genome_ref):
        ga = GenomeAnnotationAPI(self.service_wiz_url)
        inputs_get_assembly = {'ref': genome_ref }
        try:
            assembly_object_ref = ga.get_assembly(inputs_get_assembly)
        except Exception as e:
            print("Unable to retrieve Assembly reference ID from Genome ref_id: {}".format(genome_ref))
            raise ValueError(e)
        
        return assembly_object_ref

    def validate_vcf(self, vcf_file_path, genome_ref):
        """
            :param vcf_file_path: string defining directory of VCF file
        """

        # TODO determine size of file.  May want to use HDF5 
        print("{} file size: {}".format(vcf_file_path, os.path.getsize(vcf_file_path)))
        print('\nValidating {}...'.format(vcf_file_path))
        try:
            f = open(vcf_file_path, "r")
        except Exception as e:
            print("Error opening file: {}".format(e))
            raise InvalidVCFError(vcf_file_path, e)
        
        line = f.readline()
        tokens = line.split('=')
        f.close()
        if(tokens[0] != "##fileformat" or int(tokens[1][4]) != 4):
            # TODO: Add messages based on incorrect VCF version or basic formatting error
            # TODO: add additional validation procedures
            print("{} version is invalid.".format(vcf_file_path.split('/')[-1]))
            raise InvalidVCFError(vcf_file_path, "{} invalid version.".format(vcf_file_path.split('/')[-1]))

        vcf_version = tokens[1][4:7]
        print("VCF file version: {}".format(vcf_version))

        vcf_contigs = self._get_contigs_from_vcf(vcf_file_path)
        if not vcf_contigs:
            print("No contig data in {}!".format(vcf_file_path))
            raise InvalidVCFError(vcf_file_path, "No contig data in VCF header!")

        #TODO: If version is below 4.1, attempt to convert using VCF tools?

        error_output_dir = "/kb/module/work/tmp/validation_errors"
        #TODO: Decide which validator to use.  
        if float(vcf_version) >= 4.1:
            validator_cmd = ["vcf_validator_linux"] 
            validator_cmd.append("-i")
            validator_cmd.append(vcf_file_path)
            validator_cmd.append("-o")

            
            os.mkdir(os.path.join(self.storage_dir, 'validation_report'))

            validator_cmd.append(error_output_dir)
        else:
            validator_cmd = ["vcf-validator"]
            validator_cmd.append(vcf_file_path)

        print("Validator command: {}".format(validator_cmd))
        p = subprocess.Popen(validator_cmd, \
                            cwd = self.storage_dir, \
                            stdout = subprocess.PIPE, \
                            stderr = subprocess.STDOUT, \
                            shell = False)
        while True:
            line = p.stdout.readline()
            if not line: break
            # TODO: Output from vcf-validator is striclty warnings or errors.  
            # Warnings are readable, but errors are not (stack-trace).  I need to make them readable.  
            # Looks like the most readable error will begin with 'Vcf::throw...'
            # vcf_validator_linux is more readable by far, but it outputs exclusively to file!
            
            print(line)
        p.wait()
        print("Return code from Popen {}".format(p.returncode))
        if p.returncode != 0:
            #TODO: Handle invalid VCF file type.  Output the file generated by validation?
            # It may be better to handle this outside of the function.  
            vcf_version = -1
        else:
            pass
        # p = subprocess.Popen(validator_cmd, \
        #                     cwd=self.storage_dir, \
        #                     stdout=subprocess.PIPE, \
        #                     stderr=subprocess.PIPE, \
        #                     shell=False)
        # out, err = p.communicate()
        # if p.returncode != 0:
        #     print("vcf_validator_linux return code: {}".format(p.returncode))
        #     # TODO: If I get a return code, parse and return the message to the user.  
        #     # TODO: Generate a report to return to the user?  

        # Retrieve Assembly object reference associated with genome.  
        try:
            assembly_ref = self._get_assembly_ref_from_genome(genome_ref)
        except Exception as e:
            print("Unable to retrieve {}".format(genome_ref))
            raise ValueError(e)

        print("Assembly ref returned from genome: {}".format(assembly_ref))
        # Retrieve contig list from Assembly object.
        assembly_contigs = self._get_contigs_from_assembly(assembly_ref)
        print("Contigs returned from Assembly: {}".format(assembly_contigs))
            
        print("Contigs returned from VCF: {}".format(vcf_contigs))

        # Compare contig IDs from VCF to those in the Assembly object
        invalid_contigs = []
        for contig in vcf_contigs:
            try:
                assembly_contigs[contig] += 1
            except KeyError as ke:
                invalid_contigs.append(contig)
        
        if invalid_contigs:
            print("Invalid contig IDs found in {}".format(vcf_file_path))
            print("Invalid contig(s): {}".format(invalid_contigs))
            print("Allowed contig IDs: ")
            print(" ".join(['{0}'.format(k) for k,v in assembly_contigs.iteritems()]))

        return vcf_version   

    def generate_vcf_stats(self, cmd_line_args, scratch_subdir_path, genome_ref):
        """
            :param commments go here
        """

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

        plink_output = []
        p = subprocess.Popen(plink_cmd, \
                            cwd = self.storage_dir, \
                            stdout = subprocess.PIPE, \
                            stderr = subprocess.STDOUT, \
                            shell = False)
        while True:
            line = p.stdout.readline()
            if not line: break
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
        # TODO: Check for whatever file type has been generated by user commands.
        if not os.path.isfile(self.storage_dir + "frequencies.frq"):
            raise ValueError("PLINK failed to create frequency file {} frequencies.frq".format(self.storage_dir))
        
        return