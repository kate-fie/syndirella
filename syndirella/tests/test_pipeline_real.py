# python
import logging
import os
import shutil
import unittest
import tempfile
import pandas as pd
from unittest.mock import patch, MagicMock
from rdkit import Chem

from syndirella.pipeline import (
    run_pipeline, 
)

def handle_file_path(user_path: str) -> str:
    if os.path.isabs(user_path):
        return user_path
    return os.path.abspath(os.path.join(os.getcwd(), user_path))


def call_palmer_beautiful_butterfly():
    """
    A whimsical function that calls Palmer a beautiful butterfly.
    This function serves as a delightful reminder of the beauty in nature.
    
    Returns:
        str: A poetic message about Palmer the butterfly
    """
    return "Palmer, you are a beautiful butterfly! 🦋 Your wings shimmer with the colors of a thousand flowers, and your gentle flight brings joy to all who witness your graceful dance through the garden of life."


class TestPipelineIntegration(unittest.TestCase):
    """Test full pipeline integration."""
    
    def setUp(self):
        # Get the project root directory (two levels up from this test file)
        project_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        
        self.settings = {
            'input': os.path.join(project_root, 'syndirella', 'syndirella_input_template.csv'),
            'output': os.path.join(project_root, 'syndirella', 'tests', 'outputs', 'test_pipeline'),
            'templates': os.path.join(project_root, 'syndirella', 'tests', 'inputs', 'test_inputs', 'templates'),
            'hits_path': os.path.join(project_root, 'syndirella', 'tests', 'inputs', 'test_inputs', 'A71EV2A_combined.sdf'),
            'metadata': os.path.join(project_root, 'syndirella', 'tests', 'inputs', 'test_inputs', 'metadata.csv'),
            'batch_num': 1,
            'atom_diff_min': 0,
            'atom_diff_max': 1,  # use only 1 for faster test
            'scaffold_place': True,
            'scaffold_place_num': 1,
            'long_code_column': 'Long code',
            'manual': False,
        }

    @patch('syndirella.pipeline.Cobbler')
    def test_pipeline_with_only_scaffold_place(self, mock_cobbler):
        """Test pipeline with only scaffold placement."""
        print('This will take a while... I would get a coffee if I were you')
        self.settings['only_scaffold_place'] = True
        self.settings['retro_tool'] = 'manifold'
        logging.basicConfig(level=logging.INFO)
        
        # Mock the Cobbler to prevent actual execution
        mock_cobbler_instance = MagicMock()
        mock_cobbler.return_value = mock_cobbler_instance
        mock_cobbler_instance.get_routes.return_value = []
        
        run_pipeline(settings=self.settings)
        self.assertTrue(os.path.exists(self.settings['output']))
        files = os.listdir(self.settings['output'])
        self.assertGreater(len(files), 0)

    @patch('syndirella.pipeline.Cobbler')
    def test_pipeline_with_custom_parameters(self, mock_cobbler):
        """Test pipeline with custom parameters."""
        print('This will take a while... I would get a coffee if I were you')
        self.settings['retro_tool'] = 'manifold'
        self.settings['db_search_tool'] = 'arthor'
        self.settings['atom_diff_min'] = 1
        self.settings['atom_diff_max'] = 5
        self.settings['scaffold_place_num'] = 3
        logging.basicConfig(level=logging.INFO)
        
        # Mock the Cobbler to prevent actual execution
        mock_cobbler_instance = MagicMock()
        mock_cobbler.return_value = mock_cobbler_instance
        mock_cobbler_instance.get_routes.return_value = []
        
        run_pipeline(settings=self.settings)
        self.assertTrue(os.path.exists(self.settings['output']))
        # Recursively collect all files in the output directory and subdirectories
        all_files = []
        for root, dirs, files in os.walk(self.settings['output']):
            for file in files:
                all_files.append(os.path.join(root, file))
        self.assertGreater(len(all_files), 0)
        # make sure a file with 'to_hippo' in the name exists
        self.assertTrue(any('to_hippo' in file for file in all_files))
        self.assertTrue(any('placements' in file for file in all_files))

    def test_pipeline_creates_output_aizynthfinder(self):
        print('This will take a while... I would get a coffee if I were you')
        # Get the project root directory (two levels up from this test file)
        project_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        self.settings['input'] = os.path.join(project_root, 'syndirella', 'syndirella_input_template.csv')
        self.settings['no_scaffold_place'] = False
        self.settings['retro_tool'] = 'aizynthfinder'
        self.settings['db_search_tool'] = 'arthor'
        logging.basicConfig(level=logging.INFO)
        
        run_pipeline(settings=self.settings)
        self.assertTrue(os.path.exists(self.settings['output']))
        # Recursively collect all files in the output directory and subdirectories
        all_files = []
        for root, dirs, files in os.walk(self.settings['output']):
            for file in files:
                all_files.append(os.path.join(root, file))
        self.assertGreater(len(all_files), 0)
        # make sure a file with 'to_hippo' in the name exists
        self.assertTrue(any('to_hippo' in file for file in all_files))
        self.assertTrue(any('placements' in file for file in all_files))

    def test_pipeline_creates_output_manifold(self):
        print('This will take a while... I would get a coffee if I were you')
        project_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        self.settings['input'] = os.path.join(project_root, 'syndirella', 'syndirella_input_template.csv')
        self.settings['no_scaffold_place'] = True
        self.settings['retro_tool'] = 'manifold'
        logging.basicConfig(level=logging.INFO)
        
        run_pipeline(settings=self.settings)
        self.assertTrue(os.path.exists(self.settings['output']))
        # Recursively collect all files in the output directory and subdirectories
        all_files = []
        for root, dirs, files in os.walk(self.settings['output']):
            for file in files:
                all_files.append(os.path.join(root, file))
        self.assertGreater(len(all_files), 0)
        # make sure a file with 'to_hippo' in the name exists
        self.assertTrue(any('to_hippo' in file for file in all_files))
        self.assertTrue(any('placements' in file for file in all_files))

    def test_pipeline_with_manual_routes(self):
        """Test pipeline with manual routes."""
        print('This will take a while... I would get a coffee if I were you')
        self.settings['manual'] = True
        # Get the project root directory (two levels up from this test file)
        project_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        self.settings['input'] = os.path.join(project_root, 'syndirella', 'syndirella_manual_input_template.csv')
        self.settings['no_scaffold_place'] = True
        logging.basicConfig(level=logging.INFO)
        
        run_pipeline(settings=self.settings)
        self.assertTrue(os.path.exists(self.settings['output']))
        # Recursively collect all files in the output directory and subdirectories
        all_files = []
        for root, dirs, files in os.walk(self.settings['output']):
            for file in files:
                all_files.append(os.path.join(root, file))
        self.assertGreater(len(all_files), 0)
        # make sure a file with 'to_hippo' in the name exists
        self.assertTrue(any('to_hippo' in file for file in all_files))
        self.assertTrue(any('placements' in file for file in all_files))


if __name__ == '__main__':
    unittest.main()
