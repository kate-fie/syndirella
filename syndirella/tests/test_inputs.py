#!/usr/bin/env python3
"""
Minimal comprehensive test for input validation functionality.
"""
import logging
import os
import shutil
import sys
import unittest
from unittest.mock import patch, MagicMock
import tempfile
import pandas as pd

from syndirella.utils import check_inputs


def abs_path(*parts):
    """Helper to get absolute path from test file location."""
    return os.path.abspath(os.path.join(os.path.dirname(__file__), *parts))


class TestInputValidations(unittest.TestCase):
    """Minimal comprehensive test for input validation functionality."""
    
    def setUp(self):
        self.csv_path = abs_path('inputs', 'test_inputs', 'syndirella_input_template.csv')
        self.manual_csv_path = abs_path('inputs', 'test_inputs', 'syndirella_manual_input_template.csv')
        self.template_dir = abs_path('inputs', 'test_inputs', 'templates')
        self.hits_path = abs_path('inputs', 'test_inputs', 'A71EV2A_combined.sdf')
        self.metadata = abs_path('inputs', 'test_inputs', 'metadata.csv')
        
        # Create temporary test files
        self.temp_dir = tempfile.mkdtemp()
        self.temp_csv = os.path.join(self.temp_dir, 'test.csv')
        self.temp_metadata = os.path.join(self.temp_dir, 'metadata.csv')
        self.temp_template_dir = os.path.join(self.temp_dir, 'templates')
        self.temp_hits_path = os.path.join(self.temp_dir, 'hits.sdf')
        
        # Create test template directory
        os.makedirs(self.temp_template_dir, exist_ok=True)

    def tearDown(self):
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    def test_check_csv_valid(self):
        """Test valid CSV file validation."""
        check_inputs.check_csv(self.csv_path)
        self.assertTrue(os.path.exists(self.csv_path))

    def test_check_csv_manual_valid(self):
        """Test valid manual CSV file validation."""
        check_inputs.check_csv(self.manual_csv_path)
        self.assertTrue(os.path.exists(self.manual_csv_path))

    def test_check_csv_file_not_found(self):
        """Test CSV file not found error."""
        with self.assertRaises(FileNotFoundError):
            check_inputs.check_csv('nonexistent_file.csv')

    def test_check_csv_missing_required_columns(self):
        """Test CSV with missing required columns."""
        df = pd.DataFrame({'smiles': ['CCO'], 'compound_set': ['test']})
        df.to_csv(self.temp_csv, index=False)
        
        with self.assertRaises(ValueError):
            check_inputs.check_csv(self.temp_csv)


    def test_check_template_paths_valid(self):
        """Test template paths validation."""
        result = check_inputs.check_template_paths(self.template_dir, self.csv_path)
        self.assertIsInstance(result, set)
        self.assertGreater(len(result), 0)

    def test_check_template_paths_invalid_directory(self):
        """Test template directory not found error."""
        with self.assertRaises(NotADirectoryError):
            check_inputs.check_template_paths('nonexistent_dir', self.csv_path)

    def test_check_manual_valid(self):
        """Test manual route validation."""
        check_inputs.check_manual(self.manual_csv_path)

    def test_check_manual_missing_columns(self):
        """Test manual CSV with missing required columns."""
        df = pd.DataFrame({'smiles': ['CCO'], 'compound_set': ['test']})
        df.to_csv(self.temp_csv, index=False)
        
        with self.assertRaises(ValueError):
            check_inputs.check_manual(self.temp_csv)

    def test_check_apo_template_valid(self):
        """Test apo template validation."""
        pdb_content = """ATOM      1  N   ALA A   1      27.044  14.477   5.456  1.00 20.00           N
ATOM      2  CA  ALA A   1      26.206  13.246   5.456  1.00 20.00           C
ATOM      3  C   ALA A   1      26.206  12.456   6.756  1.00 20.00           C
ATOM      4  O   ALA A   1      27.044  12.456   7.756  1.00 20.00           O
END"""
        pdb_file = os.path.join(self.temp_template_dir, 'test.pdb')
        with open(pdb_file, 'w') as f:
            f.write(pdb_content)
        
        check_inputs.check_apo_template(pdb_file)

    def test_check_apo_template_invalid_file(self):
        """Test apo template file not found error."""
        with self.assertRaises(FileNotFoundError):
            check_inputs.check_apo_template('nonexistent.pdb')

    def test_format_additional_info(self):
        """Test additional info formatting."""
        row = pd.Series({
            'compound_set': 'test_set',
            'smiles': 'CCO',
            'template': 'test_template',
            'hit1': 'test_hit'
        })
        additional_columns = ['compound_set']
        
        result = check_inputs.format_additional_info(row, additional_columns)
        self.assertIsInstance(result, dict)
        self.assertIn('compound_set', result)
        self.assertEqual(result['compound_set'], 'test_set')

    def test_get_exact_hit_names(self):
        """Test exact hit names retrieval."""
        row = pd.Series({
            'hit1': 'Ax0556a',
            'hit2': 'Ax0450a'
        })
        
        result = check_inputs.get_exact_hit_names(row, self.hits_path)
        self.assertIsInstance(result, list)
        self.assertGreater(len(result), 0)

    def test_get_exact_hit_names_invalid(self):
        """Test exact hit names retrieval with invalid hit name."""
        row = pd.Series({
            'hit1': 'InvalidHitName',
            'hit2': 'Ax0450a'
        })
        
        with self.assertRaises(ValueError) as context:
            check_inputs.get_exact_hit_names(row, self.hits_path)
        
        self.assertIn("not found in SDF file", str(context.exception))

    def test_get_template_path(self):
        """Test template path retrieval."""
        template = 'Ax0310a'
        result = check_inputs.get_template_path(self.template_dir, template)
        self.assertIsInstance(result, str)
        self.assertTrue(os.path.exists(result))

    def test_get_template_path_not_found(self):
        """Test template path not found error."""
        template = 'nonexistent_template'
        with self.assertRaises(FileNotFoundError):
            check_inputs.get_template_path(self.template_dir, template)

    def test_check_additional_columns(self):
        """Test additional columns validation."""
        check_inputs.check_additional_columns(self.csv_path, ['compound_set'])

    def test_check_additional_columns_missing(self):
        """Test missing additional columns error."""
        with self.assertRaises(ValueError):
            check_inputs.check_additional_columns(self.csv_path, ['nonexistent_column'])

    def test_format_manual_route(self):
        """Test manual route formatting."""
        row = pd.Series({
            'reactant_step1': 'CCO',
            'reactant2_step1': 'CCCO',
            'reaction_name_step1': 'test_reaction_1'
        })
        
        reactants, reactions, num_steps = check_inputs.format_manual_route(row)
        self.assertIsInstance(reactants, list)
        self.assertIsInstance(reactions, list)
        self.assertIsInstance(num_steps, int)
        self.assertEqual(len(reactants), 1)
        self.assertEqual(len(reactions), 1)
        self.assertEqual(num_steps, 1)

    def test_check_pipeline_inputs_valid(self):
        """Test valid pipeline inputs validation."""
        check_inputs.check_pipeline_inputs(
            csv_path=self.csv_path,
            template_dir=self.template_dir,
            hits_path=self.hits_path,
            additional_columns=['compound_set'],
            manual_routes=False
        )

    def test_check_pipeline_inputs_manual_routes(self):
        """Test manual pipeline inputs validation."""
        check_inputs.check_pipeline_inputs(
            csv_path=self.manual_csv_path,
            template_dir=self.template_dir,
            hits_path=self.hits_path,
            additional_columns=['compound_set'],
            manual_routes=True
        )



if __name__ == '__main__':
    unittest.main()
