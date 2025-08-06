# python
import logging
import os
import shutil
import sys
import unittest
from unittest.mock import patch, MagicMock
import tempfile
import pandas as pd

import syndirella.check_inputs as check_inputs
from syndirella.cli import main


class TestInputValidations(unittest.TestCase):
    def setUp(self):
        self.csv_path = 'syndirella_input_template.csv'
        self.manual_csv_path = 'syndirella_manual_input_template.csv'
        self.template_dir = 'tests/inputs/test_inputs/templates'
        self.hits_path = 'tests/inputs/test_inputs/A71EV2A_combined.sdf'
        self.metadata = 'tests/inputs/test_inputs/metadata.csv'
        
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
        logging.basicConfig(level=logging.INFO)
        check_inputs.check_csv(self.csv_path)
        self.assertTrue(os.path.exists(self.csv_path))

    def test_check_csv_manual_valid(self):
        logging.basicConfig(level=logging.INFO)
        check_inputs.check_csv(self.manual_csv_path)
        self.assertTrue(os.path.exists(self.manual_csv_path))

    def test_check_csv_file_not_found(self):
        with self.assertRaises(FileNotFoundError):
            check_inputs.check_csv('nonexistent_file.csv')

    def test_check_csv_missing_required_columns(self):
        # Create CSV with missing required columns
        df = pd.DataFrame({'smiles': ['CCO'], 'compound_set': ['test']})  # Missing 'template' and 'hit1'
        df.to_csv(self.temp_csv, index=False)
        
        with self.assertRaises(ValueError):
            check_inputs.check_csv(self.temp_csv)

    def test_check_hit_names(self):
        logging.basicConfig(level=logging.INFO)
        check_inputs.check_hit_names(self.csv_path, self.hits_path, self.metadata, 'Long code')

    def test_check_hit_names_invalid_metadata(self):
        with self.assertRaises(FileNotFoundError):
            check_inputs.check_hit_names(self.csv_path, self.hits_path, 'nonexistent_metadata.csv', 'Long code')

    def test_metadata_dict_valid(self):
        result = check_inputs.metadata_dict(self.metadata, 'Long code')
        self.assertIsInstance(result, dict)
        self.assertGreater(len(result), 0)

    def test_metadata_dict_file_not_found(self):
        with self.assertRaises(FileNotFoundError):
            check_inputs.metadata_dict('nonexistent_metadata.csv', 'Long code')

    def test_metadata_dict_missing_columns(self):
        # Create metadata CSV with missing required columns
        df = pd.DataFrame({'other_column': ['test']})
        df.to_csv(self.temp_metadata, index=False)
        
        with self.assertRaises(ValueError):
            check_inputs.metadata_dict(self.temp_metadata, 'Long code')

    def test_metadata_dict_with_crystal_name_fallback(self):
        # Create metadata CSV with crystal_name column
        df = pd.DataFrame({'crystal_name': ['test1', 'test2']})
        df.to_csv(self.temp_metadata, index=False)
        
        result = check_inputs.metadata_dict(self.temp_metadata, 'Long code')
        self.assertEqual(result, {'test1': 'test1', 'test2': 'test2'})

    def test_check_template_paths_valid(self):
        result = check_inputs.check_template_paths(self.template_dir, self.csv_path, self.metadata)
        self.assertIsInstance(result, set)
        self.assertGreater(len(result), 0)

    def test_check_template_paths_invalid_directory(self):
        with self.assertRaises(NotADirectoryError):
            check_inputs.check_template_paths('nonexistent_dir', self.csv_path, self.metadata)

    def test_check_template_paths_no_templates_in_csv(self):
        # Create CSV without template column
        df = pd.DataFrame({'smiles': ['CCO'], 'compound_set': ['test'], 'hit1': ['test']})
        df.to_csv(self.temp_csv, index=False)
        
        with self.assertRaises(ValueError):
            check_inputs.check_template_paths(self.template_dir, self.temp_csv, self.metadata)

    def test_check_manual_valid(self):
        check_inputs.check_manual(self.manual_csv_path)

    def test_check_manual_missing_columns(self):
        # Create manual CSV with missing required columns
        df = pd.DataFrame({'smiles': ['CCO'], 'compound_set': ['test']})  # Missing manual route columns
        df.to_csv(self.temp_csv, index=False)
        
        with self.assertRaises(ValueError):
            check_inputs.check_manual(self.temp_csv)

    def test_check_apo_template_valid(self):
        # Create a mock PDB file
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
        with self.assertRaises(FileNotFoundError):
            check_inputs.check_apo_template('nonexistent.pdb')

    def test_format_additional_info(self):
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
        row = pd.Series({
            'hit1': 'Ax0556a',
            'hit2': 'Ax0450a'
        })
        
        result = check_inputs.get_exact_hit_names(row, self.metadata, self.hits_path)
        self.assertIsInstance(result, list)
        self.assertGreater(len(result), 0)

    def test_get_template_path(self):
        template = 'Ax0310a'
        result = check_inputs.get_template_path(self.template_dir, template, self.metadata)
        self.assertIsInstance(result, str)
        self.assertTrue(os.path.exists(result))

    def test_get_template_path_not_found(self):
        template = 'nonexistent_template'
        with self.assertRaises(FileNotFoundError):
            check_inputs.get_template_path(self.template_dir, template, self.metadata)

    def test_check_additional_columns(self):
        check_inputs.check_additional_columns(self.csv_path, ['compound_set'])

    def test_check_additional_columns_missing(self):
        with self.assertRaises(ValueError):
            check_inputs.check_additional_columns(self.csv_path, ['nonexistent_column'])

    def test_format_manual_route(self):
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
        check_inputs.check_pipeline_inputs(
            csv_path=self.csv_path,
            template_dir=self.template_dir,
            hits_path=self.hits_path,
            metadata_path=self.metadata,
            additional_columns=['compound_set'],
            manual_routes=False,
            long_code_column='Long code'
        )

    def test_check_pipeline_inputs_manual_routes(self):
        check_inputs.check_pipeline_inputs(
            csv_path=self.manual_csv_path,
            template_dir=self.template_dir,
            hits_path=self.hits_path,
            metadata_path=self.metadata,
            additional_columns=['compound_set'],
            manual_routes=True,
            long_code_column='Long code'
        )


class TestCLI(unittest.TestCase):
    def setUp(self):
        self.csv_path = 'syndirella_input_template.csv'
        self.output_dir = 'tests/outputs/test_inputs/cli'
        self.manual_csv_path = 'syndirella_manual_input_template.csv'
        self.template_dir = 'tests/inputs/test_inputs/templates'
        self.hits_path = 'tests/inputs/test_inputs/A71EV2A_combined.sdf'
        self.metadata = 'tests/inputs/test_inputs/metadata.csv'
        if os.path.exists(self.output_dir):
            shutil.rmtree(self.output_dir)

    def run_main_with_args(self, args):
        with patch.object(sys, 'argv', args):
            # Mock the pipeline execution to avoid running the full pipeline
            with patch('syndirella.cli.run_pipeline') as mock_run_pipeline, \
                 patch('syndirella.cli.run_justretroquery') as mock_run_justretroquery:
                mock_run_pipeline.return_value = None
                mock_run_justretroquery.return_value = None
                try:
                    main()
                except SystemExit as e:
                    return e.code
                return 0

    def test_valid_pipeline_arguments(self):
        args = ['syndirella', '--input', self.csv_path, '--output', self.output_dir, 
                '--hits_path', self.hits_path, '--metadata', self.metadata, '--templates', self.template_dir]
        exit_code = self.run_main_with_args(args)
        self.assertEqual(exit_code, 0)

    def test_aizynth_pipeline_argument(self):
        args = ['syndirella', '--input', self.csv_path, '--output', self.output_dir, '--retro_tool', 'aizynthfinder', 
                '--hits_path', self.hits_path, '--metadata', self.metadata, '--templates', self.template_dir]
        exit_code = self.run_main_with_args(args)
        self.assertEqual(exit_code, 0)

    def test_manifold_pipeline_argument(self):
        args = ['syndirella', '--input', self.csv_path, '--output', self.output_dir, '--retro_tool', 'manifold',
                '--hits_path', self.hits_path, '--metadata', self.metadata, '--templates', self.template_dir]
        exit_code = self.run_main_with_args(args)
        self.assertEqual(exit_code, 0)

    def test_valid_just_retro_arguments(self):
        args = ['syndirella', '--input', self.csv_path, '--output', self.output_dir, '--just_retro']
        exit_code = self.run_main_with_args(args)
        self.assertEqual(exit_code, 0)

    def test_manual_routes_argument(self):
        args = ['syndirella', '--input', self.manual_csv_path, '--output', self.output_dir, '--manual',
                '--hits_path', self.hits_path, '--metadata', self.metadata, '--templates', self.template_dir]
        exit_code = self.run_main_with_args(args)
        self.assertEqual(exit_code, 0)

    def test_only_scaffold_place_argument(self):
        args = ['syndirella', '--input', self.csv_path, '--output', self.output_dir, '--only_scaffold_place',
                '--hits_path', self.hits_path, '--metadata', self.metadata, '--templates', self.template_dir]
        exit_code = self.run_main_with_args(args)
        self.assertEqual(exit_code, 0)

    def test_no_scaffold_place_argument(self):
        args = ['syndirella', '--input', self.csv_path, '--output', self.output_dir, '--no_scaffold_place',
                '--hits_path', self.hits_path, '--metadata', self.metadata, '--templates', self.template_dir]
        exit_code = self.run_main_with_args(args)
        self.assertEqual(exit_code, 0)

    def test_custom_parameters(self):
        args = [
            'syndirella', 
            '--input', self.csv_path, 
            '--output', self.output_dir,
            '--hits_path', self.hits_path, 
            '--metadata', self.metadata, 
            '--templates', self.template_dir,
            '--batch_num', '5',
            '--atom_diff_min', '1',
            '--atom_diff_max', '5',
            '--scaffold_place_num', '3'
        ]
        exit_code = self.run_main_with_args(args)
        self.assertEqual(exit_code, 0)

    def test_missing_required_arguments(self):
        args = ['syndirella', '--input', self.csv_path]
        with patch.object(sys, 'argv', args):
            with self.assertRaises(SystemExit):
                main()

    def test_missing_env_variables(self):
        args = ['syndirella', '--input', self.csv_path, '--output', 'outputs',
                '--hits_path', self.hits_path, '--metadata', self.metadata, '--templates', self.template_dir]
        with patch.object(sys, 'argv', args):
            with self.assertRaises(SystemExit):
                main()

    def test_invalid_retro_tool(self):
        args = ['syndirella', '--input', self.csv_path, '--output', self.output_dir, '--retro_tool', 'invalid_tool',
                '--hits_path', self.hits_path, '--metadata', self.metadata, '--templates', self.template_dir]
        with patch.object(sys, 'argv', args):
            with self.assertRaises(SystemExit):
                main()

    def test_invalid_db_search_tool(self):
        args = ['syndirella', '--input', self.csv_path, '--output', self.output_dir, '--db_search_tool', 'invalid_tool',
                '--hits_path', self.hits_path, '--metadata', self.metadata, '--templates', self.template_dir]
        with patch.object(sys, 'argv', args):
            with self.assertRaises(SystemExit):
                main()

    def test_cli_argument_parsing(self):
        """Test that CLI arguments are parsed correctly without running the pipeline."""
        args = ['syndirella', '--input', self.csv_path, '--output', self.output_dir, 
                '--hits_path', self.hits_path, '--metadata', self.metadata, '--templates', self.template_dir,
                '--retro_tool', 'manifold', '--db_search_tool', 'arthor', '--batch_num', '10']
        
        with patch.object(sys, 'argv', args):
            with patch('syndirella.cli.run_pipeline') as mock_run_pipeline:
                mock_run_pipeline.return_value = None
                main()
                
                # Verify that run_pipeline was called with the correct settings
                mock_run_pipeline.assert_called_once()
                call_args = mock_run_pipeline.call_args[0][0]  # First argument is settings dict
                self.assertEqual(call_args['input'], self.csv_path)
                self.assertEqual(call_args['output'], self.output_dir)
                self.assertEqual(call_args['hits_path'], self.hits_path)
                self.assertEqual(call_args['metadata'], self.metadata)
                self.assertEqual(call_args['templates'], self.template_dir)
                self.assertEqual(call_args['retro_tool'], 'manifold')
                self.assertEqual(call_args['db_search_tool'], 'arthor')
                self.assertEqual(call_args['batch_num'], 10)

    def test_just_retro_argument_parsing(self):
        """Test that just_retro arguments are parsed correctly."""
        args = ['syndirella', '--input', self.csv_path, '--output', self.output_dir, '--just_retro',
                '--retro_tool', 'aizynthfinder']
        
        with patch.object(sys, 'argv', args):
            with patch('syndirella.cli.run_justretroquery') as mock_run_justretroquery:
                mock_run_justretroquery.return_value = None
                try:
                    main()
                    exit_code = 0
                except SystemExit as e:
                    exit_code = e.code
                
                # Verify that run_justretroquery was called with the correct settings
                mock_run_justretroquery.assert_called_once()
                call_args = mock_run_justretroquery.call_args[0][0]  # First argument is settings dict
                self.assertEqual(call_args['input'], self.csv_path)
                self.assertEqual(call_args['output'], self.output_dir)
                self.assertEqual(call_args['retro_tool'], 'aizynthfinder')
                self.assertEqual(exit_code, 0)


if __name__ == '__main__':
    unittest.main()
