# -*- coding: utf-8 -*-
import unittest
import os  # noqa: F401
import json  # noqa: F401
import time
import requests
import shutil
import uuid
from os import environ
try:
    from ConfigParser import ConfigParser  # py2
except:
    from configparser import ConfigParser  # py3

from pprint import pprint  # noqa: F401

from DataFileUtil.DataFileUtilClient import DataFileUtil
from mock import patch
from biokbase.workspace.client import Workspace as workspaceService
from kb_variation_importer.kb_variation_importerImpl import kb_variation_importer
from kb_variation_importer.kb_variation_importerServer import MethodContext
from kb_variation_importer.authclient import KBaseAuth as _KBaseAuth

mock_assembly = {
            "assembly_id": "Carsonella_ruddii_HT.fna.gz_assembly",
            "base_counts": {
                "A": 67508,
                "C": 11789,
                "G": 11134,
                "T": 67112
            },
            "contigs": {
                "CP003544.1": {
                    "contig_id": "CP003544.1",
                    "description": "Candidatus Carsonella ruddii HT isolate Thao2000, complete genome",
                    "gc_content": 0.1455,
                    "length": 157543,
                    "md5": "2648e704354959e79f5de6fff3b5b9db",
                    "name": "CP003544.1"
                }
            },
            "dna_size": 157543,
            "gc_content": 0.1455,
            "num_contigs": 1,
            "type": "Unknown"
        }

class kb_variation_importerTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        token = environ.get('KB_AUTH_TOKEN', None)
        config_file = environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('kb_variation_importer'):
            cls.cfg[nameval[0]] = nameval[1]
        # Getting username from Auth profile for token
        authServiceUrl = cls.cfg['auth-service-url']
        auth_client = _KBaseAuth(authServiceUrl)
        user_id = auth_client.get_user(token)
        # WARNING: don't call any logging methods on the context object,
        # it'll result in a NoneType error
        cls.ctx = MethodContext(None)
        cls.ctx.update({'token': token,
                        'user_id': user_id,
                        'provenance': [
                            {'service': 'kb_variation_importer',
                             'method': 'please_never_use_it_in_production',
                             'method_params': []
                             }],
                        'authenticated': 1})
        cls.wsURL = cls.cfg['workspace-url']
        cls.wsClient = workspaceService(cls.wsURL)
        cls.serviceImpl = kb_variation_importer(cls.cfg)
        cls.scratch = cls.cfg['scratch']
        cls.callback_url = os.environ['SDK_CALLBACK_URL']

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    def getWsClient(self):
        return self.__class__.wsClient

    def getWsName(self):
        if hasattr(self.__class__, 'wsName'):
            return self.__class__.wsName
        suffix = int(time.time() * 1000)
        wsName = "test_kb_variation_importer_" + str(suffix)
        ret = self.getWsClient().create_workspace({'workspace': wsName})  # noqa
        self.__class__.wsName = wsName
        return wsName

    def getImpl(self):
        return self.__class__.serviceImpl

    def getContext(self):
        return self.__class__.ctx

    @staticmethod
    def fake_staging_download(params):
        scratch = '/kb/module/work/tmp/'
        inpath = params['staging_file_subdir_path']
        shutil.copy('/kb/module/data/'+inpath, scratch+inpath)
        return {'copy_file_path': scratch+inpath}
    # NOTE: According to Python unittest naming rules test method names should start from 'test'. # noqa
    @patch.object(DataFileUtil, "download_staging_file",
                  new=fake_staging_download)

    def _save_to_ws_and_report(self, ws_id, source, assembly_data):
        dfu = DataFileUtil(os.environ['SDK_CALLBACK_URL'])
        info = dfu.save_objects(
            {
                'id': '18590', # Numerical id of workspace
                "objects": [{
                    "type": "KBaseGenomeAnnotations.Assembly-3.0",
                    "data": assembly_data,
                    "name": ws_id
                }]
            })[0]
        #print("Data from save to ws: {}".format(json.dumps(info, indent=2)))
        assembly_ref = "%s/%s/%s" % (info[6], info[0], info[4])

        return assembly_ref
        
    def test_your_method(self):
        # Prepare test objects in workspace if needed using
        # self.getWsClient().save_objects({'workspace': self.getWsName(),
        #                                  'objects': []})
        #
        # Run your method by
        # ret = self.getImpl().your_method(self.getContext(), parameters...)
        #
        # Check returned data with
        # self.assertEqual(ret[...], ...) or other unittest methods
        print("os environmental vars: {}".format(os.environ))
        dfu = DataFileUtil(os.environ['SDK_CALLBACK_URL'])
        upload_ret = dfu.file_to_shock({'file_path': '/kb/module/data/fasta_test.fa',
                                        'make_handle': 1,
                                        'pack' : 'zip'})
        mock_assembly['fasta_handle_info'] = upload_ret.get('handle')
        mock_assembly['fasta_handle_ref'] = upload_ret.get('handle')['hid']
        mock_assembly['md5'] = "13416456ac626d4bf786d713e5f9d7a1"
        print(json.dumps(mock_assembly, indent=2))
        assembly_ref = self._save_to_ws_and_report(self.getWsName(), "/kb/module/data/fasta_test.fa", mock_assembly)
        
        ret = dfu.get_objects(
            {'object_refs': [assembly_ref]}
        )['data'][0]['data']
        temp_dir = "/kb/module/work/tmp/" + str(uuid.uuid4())
        os.mkdir(temp_dir)
        out_dir = "{}/{}".format(temp_dir, "mock_assembly_output")
        os.mkdir(out_dir)

        with open(out_dir + '/assembly.out', 'w') as outfile:
            outfile.write(json.dumps(ret, indent=2))
        handle = dfu.package_for_download(
            {'file_path': out_dir, 'ws_refs': [assembly_ref]})
        output = {'shock_id': handle['shock_id']}
        print("Handle returned by dfu: {}".format(handle))
        # params = {
        #     'workspace_name' : self.getWsName(),
        #     'staging_file_subdir_path' : 'test.vcf',
        #     'will_perform_gwas' : 0,
        #     'command_line_args' : None
        # }
        
        # ret = self.getImpl().import_snp_data(self.getContext(), params)[0]
        # print(ret)
        # self.assertIsNotNone(ret['report_ref'], ret['report_name'])
        pass
