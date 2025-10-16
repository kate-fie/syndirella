#!/usr/bin/env python3
"""
Minimal comprehensive test for pipeline functionality with mocking.
"""
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
    PipelineConfig, 
    assert_scaffold_placement,
    start_elaboration
)
from syndirella.constants import RetrosynthesisTool, DatabaseSearchTool


def abs_path(*parts):
    """Helper to get absolute path from test file location."""
    return os.path.abspath(os.path.join(os.path.dirname(__file__), *parts))


class TestPipelineConfig(unittest.TestCase):
    """Test the PipelineConfig dataclass."""
    
    def setUp(self):
        self.settings = {
            'input': 'test_input.csv',
            'output': 'tests/outputs/test_pipeline/test_output/',
            'templates': 'tests/inputs/test_inputs/templates',
            'hits_path': 'tests/inputs/test_inputs/A71EV2A_combined.sdf',
            'batch_num': 1,
            'atom_diff_min': 0,
            'atom_diff_max': 1,
            'scaffold_place_num': 1,
            'retro_tool': 'manifold',
            'db_search_tool': 'arthor'
        }

    def test_pipeline_config_creation(self):
        """Test creating PipelineConfig from settings."""
        config = PipelineConfig.from_settings(self.settings)
        
        self.assertEqual(config.csv_path, 'test_input.csv')
        self.assertEqual(config.output_dir, 'tests/outputs/test_pipeline/test_output/')
        self.assertEqual(config.retro_tool, RetrosynthesisTool.MANIFOLD)
        self.assertEqual(config.db_search_tool, DatabaseSearchTool.ARTHOR)
        self.assertFalse(config.manual_routes)
        self.assertTrue(config.scaffold_place)

    def test_pipeline_config_with_manual_routes(self):
        """Test creating PipelineConfig with manual routes."""
        self.settings['manual'] = True
        config = PipelineConfig.from_settings(self.settings)
        
        self.assertTrue(config.manual_routes)

    def test_pipeline_config_with_only_scaffold_place(self):
        """Test creating PipelineConfig with only scaffold placement."""
        self.settings['only_scaffold_place'] = True
        config = PipelineConfig.from_settings(self.settings)
        
        self.assertTrue(config.only_scaffold_place)

    def test_pipeline_config_missing_required_field(self):
        """Test PipelineConfig creation with missing required field."""
        del self.settings['input']
        
        with self.assertRaises(KeyError):
            PipelineConfig.from_settings(self.settings)

    def test_pipeline_config_default_values(self):
        """Test PipelineConfig default values."""
        config = PipelineConfig.from_settings(self.settings)
        
        self.assertEqual(config.additional_columns, ['compound_set'])
        self.assertFalse(config.elab_single_reactant)


class TestPipelineFunctions(unittest.TestCase):
    """Test individual pipeline functions with minimal mocking."""
    
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        self.test_scaffold = "O=C(NC1CC(C(F)(F)F)C1)c1cc2ccsc2[nH]1"
        self.test_output_dir = os.path.join(self.temp_dir, "test_output")

    def tearDown(self):
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    @patch('syndirella.pipeline.SlipperFitter')
    def test_assert_scaffold_placement_valid(self, mock_slipper_fitter):
        """Test scaffold placement assertion with valid SMILES."""
        mock_fitter = MagicMock()
        mock_slipper_fitter.return_value = mock_fitter
        mock_fitter.check_scaffold.return_value = "test_placement_path"
        
        result = assert_scaffold_placement(
            scaffold=self.test_scaffold,
            template_path="test_template.pdb",
            hits_path="test_hits.sdf",
            hits_names=['Ax0310a', 'Ax0528a'],
            output_dir=self.test_output_dir,
            scaffold_place_num=3
        )
        
        self.assertIsNotNone(result)
        mock_slipper_fitter.assert_called_once()
        self.assertGreaterEqual(mock_fitter.check_scaffold.call_count, 1)

    def test_assert_scaffold_placement_invalid_smiles(self):
        """Test scaffold placement assertion with invalid SMILES."""
        with self.assertRaises(Exception):
            assert_scaffold_placement(
                scaffold="invalid_smiles",
                template_path="test_template.pdb",
                hits_path="test_hits.sdf",
                hits_names=['Ax0310a', 'Ax0528a'],
                output_dir=self.test_output_dir,
                scaffold_place_num=3
            )

    @patch('syndirella.pipeline.assert_scaffold_placement')
    def test_start_elaboration(self, mock_assert_scaffold):
        """Test start elaboration function."""
        mock_assert_scaffold.return_value = {"test": "placement"}
        
        result = start_elaboration(
            product=self.test_scaffold,
            template_path="test_template.pdb",
            hits_path="test_hits.sdf",
            hits=['Ax0310a', 'Ax0528a'],
            output_dir=self.test_output_dir,
            scaffold_place_num=3,
            scaffold_place=True
        )
        
        self.assertIsNotNone(result)
        mock_assert_scaffold.assert_called_once()


class TestPipelineErrorHandling(unittest.TestCase):
    """Test pipeline error handling."""
    
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        self.test_output_dir = os.path.join(self.temp_dir, "test_output")

    def tearDown(self):
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    def test_pipeline_with_invalid_input_file(self):
        """Test pipeline with invalid input file."""
        settings = {
            'input': 'nonexistent_file.csv',
            'output': self.test_output_dir,
            'templates': 'test/templates',
            'hits_path': 'test/hits.sdf',
            'metadata': 'test/metadata.csv',
            'batch_num': 1,
            'atom_diff_min': 0,
            'atom_diff_max': 5,
            'scaffold_place_num': 3,
        }
        
        with self.assertRaises(FileNotFoundError):
            run_pipeline(settings)

    def test_pipeline_with_invalid_settings(self):
        """Test pipeline with invalid settings."""
        settings = {
            'output': self.test_output_dir,
            # Missing required 'input' field
        }
        
        with self.assertRaises(KeyError):
            PipelineConfig.from_settings(settings)


if __name__ == '__main__':
    unittest.main() 