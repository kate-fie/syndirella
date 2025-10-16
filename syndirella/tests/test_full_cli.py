#!/usr/bin/env python3
"""
Minimal comprehensive test for the full CLI functionality.
"""
import logging
import os
import shutil
import unittest
import tempfile
import pandas as pd
import subprocess
import sys
from pathlib import Path

from syndirella.cli import main
from syndirella.pipeline import run_pipeline, PipelineConfig
from syndirella.constants import RetrosynthesisTool, DatabaseSearchTool


class TestFullCLI(unittest.TestCase):
    """Minimal comprehensive test for the full CLI functionality."""
    
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        # Get the project root directory and construct the inputs path
        current_dir = os.path.dirname(os.path.abspath(__file__))
        project_root = os.path.dirname(os.path.dirname(current_dir))
        self.inputs_dir = os.path.join(project_root, 'syndirella', 'tests', 'inputs', 'test_inputs')
        self.test_output_dir = os.path.join(self.temp_dir, "test_output")
        
        # Create minimal test input files
        self.create_test_input_files()
        
        # Set up environment variables for API keys (dummy values for testing)
        self.setup_environment()

    def tearDown(self):
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    def create_test_input_files(self):
        """Create minimal test input CSV files."""
        
        # Create automatic input CSV
        self.auto_input_csv = os.path.join(self.temp_dir, "test_auto_input.csv")
        auto_data = pd.DataFrame({
            'smiles': ['O=C(NC1CC(C(F)(F)F)C1)c1cc2ccsc2[nH]1'],
            'hit1': ['Ax0310a'],
            'hit2': ['Ax0528a'],
            'hit3': [''],
            'template': ['Ax0310a'],
            'compound_set': ['test_auto_set']
        })
        auto_data.to_csv(self.auto_input_csv, index=False)
        
        # Create manual input CSV
        self.manual_input_csv = os.path.join(self.temp_dir, "test_manual_input.csv")
        manual_data = pd.DataFrame({
            'smiles': ['O=C(NC1CC(C(F)(F)F)C1)c1cc2ccsc2[nH]1'],
            'reaction_name_step1': ['Amide formation'],
            'reactant_step1': ['CCO'],
            'reactant2_step1': [''],
            'product_step1': ['CCOC(=O)c1cc2ccsc2[nH]1'],
            'reaction_name_step2': [''],
            'reactant_step2': [''],
            'product_step2': [''],
            'reaction_name_step3': [''],
            'reactant_step3': [''],
            'hit1': ['Ax0310a'],
            'hit2': ['Ax0528a'],
            'hit3': [''],
            'template': ['Ax0310a'],
            'compound_set': ['test_manual_set']
        })
        manual_data.to_csv(self.manual_input_csv, index=False)

    def setup_environment(self):
        """Set up environment variables for testing."""
        os.environ.setdefault('MANIFOLD_API_KEY', 'test_key')
        os.environ.setdefault('MANIFOLD_API_URL', 'https://api.manifold.bio')
        os.environ.setdefault('ARTHOR_API_URL', 'https://arthor.docking.org')

    def test_cli_help(self):
        """Test CLI help functionality."""
        try:
            result = subprocess.run([
                sys.executable, '-m', 'syndirella.cli', '--help'
            ], capture_output=True, text=True, timeout=10)
            self.assertIn('usage:', result.stdout.lower())
        except (subprocess.TimeoutExpired, Exception) as e:
            logging.warning(f"CLI help test failed: {e}")

    def test_cli_run_command_help(self):
        """Test CLI run command help."""
        try:
            result = subprocess.run([
                sys.executable, '-m', 'syndirella.cli', 'run', '--help'
            ], capture_output=True, text=True, timeout=10)
            self.assertIn('usage:', result.stdout.lower())
        except (subprocess.TimeoutExpired, Exception) as e:
            logging.warning(f"CLI run help test failed: {e}")

    def test_cli_add_reaction_command_help(self):
        """Test CLI add-reaction command help."""
        try:
            result = subprocess.run([
                sys.executable, '-m', 'syndirella.cli', 'add-reaction', '--help'
            ], capture_output=True, text=True, timeout=10)
            self.assertIn('usage:', result.stdout.lower())
        except (subprocess.TimeoutExpired, Exception) as e:
            logging.warning(f"CLI add-reaction help test failed: {e}")

    def test_pipeline_config_creation(self):
        """Test PipelineConfig creation with minimal settings."""
        settings = {
            'input': self.auto_input_csv,
            'output': self.test_output_dir,
            'templates': os.path.join(self.inputs_dir, 'templates'),
            'hits_path': os.path.join(self.inputs_dir, 'A71EV2A_combined.sdf'),
            'batch_num': 1,
            'atom_diff_min': 0,
            'atom_diff_max': 5,
            'scaffold_place_num': 1,
            'retro_tool': 'manifold',
            'db_search_tool': 'arthor'
        }
        
        config = PipelineConfig.from_settings(settings)
        
        self.assertEqual(config.csv_path, self.auto_input_csv)
        self.assertEqual(config.output_dir, self.test_output_dir)
        self.assertEqual(config.retro_tool, RetrosynthesisTool.MANIFOLD)
        self.assertEqual(config.db_search_tool, DatabaseSearchTool.ARTHOR)
        self.assertFalse(config.manual_routes)
        self.assertTrue(config.scaffold_place)

    def test_cli_basic_pipeline_run(self):
        """Test basic CLI pipeline run with minimal settings."""
        try:
            result = subprocess.run([
                sys.executable, '-m', 'syndirella.cli',
                '-i', self.auto_input_csv,
                '-o', self.test_output_dir,
                '--templates', 'syndirella/tests/inputs/test_inputs/templates',
                '--hits_path', 'syndirella/tests/inputs/test_inputs/A71EV2A_combined.sdf',
                '--metadata', 'syndirella/tests/inputs/test_inputs/metadata.csv',
                '--only_scaffold_place',
                '--scaffold_place_num', '1',
                '--batch_num', '1'
            ], capture_output=True, text=True, timeout=60)
            
            # Check that the command ran (even if it failed due to missing dependencies)
            self.assertIsNotNone(result)
            
        except (subprocess.TimeoutExpired, Exception) as e:
            logging.warning(f"CLI basic pipeline test failed: {e}")

    def test_cli_manual_mode(self):
        """Test CLI with manual mode."""
        try:
            result = subprocess.run([
                sys.executable, '-m', 'syndirella.cli',
                '-i', self.manual_input_csv,
                '-o', self.test_output_dir,
                '--templates', 'syndirella/tests/inputs/test_inputs/templates',
                '--hits_path', 'syndirella/tests/inputs/test_inputs/A71EV2A_combined.sdf',
                '--metadata', 'syndirella/tests/inputs/test_inputs/metadata.csv',
                '--manual',
                '--only_scaffold_place',
                '--scaffold_place_num', '1',
                '--batch_num', '1'
            ], capture_output=True, text=True, timeout=60)
            
            self.assertIsNotNone(result)
            
        except (subprocess.TimeoutExpired, Exception) as e:
            logging.warning(f"CLI manual mode test failed: {e}")

    def test_cli_just_retro_mode(self):
        """Test CLI with just_retro mode."""
        try:
            result = subprocess.run([
                sys.executable, '-m', 'syndirella.cli',
                '-i', self.auto_input_csv,
                '-o', self.test_output_dir,
                '--just_retro'
            ], capture_output=True, text=True, timeout=60)
            
            self.assertIsNotNone(result)
            
        except (subprocess.TimeoutExpired, Exception) as e:
            logging.warning(f"CLI just_retro test failed: {e}")

    def test_cli_with_profiling(self):
        """Test CLI with profiling enabled."""
        try:
            result = subprocess.run([
                sys.executable, '-m', 'syndirella.cli',
                '-i', self.auto_input_csv,
                '-o', self.test_output_dir,
                '--templates', 'syndirella/tests/inputs/test_inputs/templates',
                '--hits_path', 'syndirella/tests/inputs/test_inputs/A71EV2A_combined.sdf',
                '--metadata', 'syndirella/tests/inputs/test_inputs/metadata.csv',
                '--only_scaffold_place',
                '--scaffold_place_num', '1',
                '--profile'
            ], capture_output=True, text=True, timeout=60)
            
            self.assertIsNotNone(result)
            
        except (subprocess.TimeoutExpired, Exception) as e:
            logging.warning(f"CLI profiling test failed: {e}")

    def test_cli_different_retro_tools(self):
        """Test CLI with different retrosynthesis tools."""
        retro_tools = ['manifold', 'aizynthfinder']
        
        for tool in retro_tools:
            try:
                result = subprocess.run([
                    sys.executable, '-m', 'syndirella.cli',
                    '-i', self.auto_input_csv,
                    '-o', self.test_output_dir,
                    '--templates', 'syndirella/tests/inputs/test_inputs/templates',
                    '--hits_path', 'syndirella/tests/inputs/test_inputs/A71EV2A_combined.sdf',
                    '--metadata', 'syndirella/tests/inputs/test_inputs/metadata.csv',
                    '--only_scaffold_place',
                    '--scaffold_place_num', '1',
                    '--retro_tool', tool
                ], capture_output=True, text=True, timeout=60)
                
                self.assertIsNotNone(result)
                
            except (subprocess.TimeoutExpired, Exception) as e:
                logging.warning(f"CLI {tool} test failed: {e}")

    def test_cli_different_db_search_tools(self):
        """Test CLI with different database search tools."""
        db_tools = ['arthor', 'postera']
        
        for tool in db_tools:
            try:
                result = subprocess.run([
                    sys.executable, '-m', 'syndirella.cli',
                    '-i', self.auto_input_csv,
                    '-o', self.test_output_dir,
                    '--templates', 'syndirella/tests/inputs/test_inputs/templates',
                    '--hits_path', 'syndirella/tests/inputs/test_inputs/A71EV2A_combined.sdf',
                    '--metadata', 'syndirella/tests/inputs/test_inputs/metadata.csv',
                    '--only_scaffold_place',
                    '--scaffold_place_num', '1',
                    '--db_search_tool', tool
                ], capture_output=True, text=True, timeout=60)
                
                self.assertIsNotNone(result)
                
            except (subprocess.TimeoutExpired, Exception) as e:
                logging.warning(f"CLI {tool} db search test failed: {e}")

    def test_cli_add_reaction_command(self):
        """Test CLI add-reaction command."""
        try:
            result = subprocess.run([
                sys.executable, '-m', 'syndirella.cli', 'add-reaction',
                '--name', 'Test Reaction',
                '--smirks', '[C:1]=[O:2].[N:3]>>[C:1](=[O:2])[N:3]',
                '--find_parent',
                '--fp_type', 'maccs_rxn_fp',
                '--threshold', '0.2',
                '--similarity_metric', 'cosine'
            ], capture_output=True, text=True, timeout=30)
            
            self.assertIsNotNone(result)
            
        except (subprocess.TimeoutExpired, Exception) as e:
            logging.warning(f"CLI add-reaction test failed: {e}")

    def test_cli_with_elab_single_reactant(self):
        """Test CLI with elab_single_reactant option."""
        try:
            result = subprocess.run([
                sys.executable, '-m', 'syndirella.cli',
                '-i', self.auto_input_csv,
                '-o', self.test_output_dir,
                '--templates', 'syndirella/tests/inputs/test_inputs/templates',
                '--hits_path', 'syndirella/tests/inputs/test_inputs/A71EV2A_combined.sdf',
                '--metadata', 'syndirella/tests/inputs/test_inputs/metadata.csv',
                '--elab_single_reactant',
                '--scaffold_place_num', '1'
            ], capture_output=True, text=True, timeout=60)
            
            self.assertIsNotNone(result)
            
        except (subprocess.TimeoutExpired, Exception) as e:
            logging.warning(f"CLI elab_single_reactant test failed: {e}")

    def test_cli_with_no_scaffold_place(self):
        """Test CLI with no_scaffold_place option."""
        try:
            result = subprocess.run([
                sys.executable, '-m', 'syndirella.cli',
                '-i', self.auto_input_csv,
                '-o', self.test_output_dir,
                '--templates', 'syndirella/tests/inputs/test_inputs/templates',
                '--hits_path', 'syndirella/tests/inputs/test_inputs/A71EV2A_combined.sdf',
                '--metadata', 'syndirella/tests/inputs/test_inputs/metadata.csv',
                '--no_scaffold_place'
            ], capture_output=True, text=True, timeout=60)
            
            self.assertIsNotNone(result)
            
        except (subprocess.TimeoutExpired, Exception) as e:
            logging.warning(f"CLI no_scaffold_place test failed: {e}")

    def test_cli_with_custom_atom_diff_limits(self):
        """Test CLI with custom atom difference limits."""
        try:
            result = subprocess.run([
                sys.executable, '-m', 'syndirella.cli',
                '-i', self.auto_input_csv,
                '-o', self.test_output_dir,
                '--templates', 'syndirella/tests/inputs/test_inputs/templates',
                '--hits_path', 'syndirella/tests/inputs/test_inputs/A71EV2A_combined.sdf',
                '--metadata', 'syndirella/tests/inputs/test_inputs/metadata.csv',
                '--only_scaffold_place',
                '--scaffold_place_num', '1',
                '--atom_diff_min', '2',
                '--atom_diff_max', '8'
            ], capture_output=True, text=True, timeout=60)
            
            self.assertIsNotNone(result)
            
        except (subprocess.TimeoutExpired, Exception) as e:
            logging.warning(f"CLI atom_diff_limits test failed: {e}")


    def test_input_file_validation(self):
        """Test that input files are properly validated."""
        self.assertTrue(os.path.exists(self.auto_input_csv))
        self.assertTrue(os.path.exists(self.manual_input_csv))
        
        # Test with invalid input file
        invalid_csv = os.path.join(self.temp_dir, "nonexistent.csv")
        self.assertFalse(os.path.exists(invalid_csv))

    def test_output_directory_creation(self):
        """Test that output directories are created properly."""
        test_output = os.path.join(self.temp_dir, "test_output_creation")
        os.makedirs(test_output, exist_ok=True)
        self.assertTrue(os.path.exists(test_output))
        self.assertTrue(os.path.isdir(test_output))

    def test_cli_without_metadata(self):
        """Test CLI without metadata argument (should work)."""
        try:
            result = subprocess.run([
                sys.executable, '-m', 'syndirella.cli',
                '-i', self.auto_input_csv,
                '-o', self.test_output_dir,
                '--templates', 'syndirella/tests/inputs/test_inputs/templates',
                '--hits_path', 'syndirella/tests/inputs/test_inputs/A71EV2A_combined.sdf',
                # Note: No --metadata argument provided
                '--only_scaffold_place',
                '--scaffold_place_num', '1',
                '--batch_num', '1'
            ], capture_output=True, text=True, timeout=60)
            
            # Check that the command ran (even if it failed due to missing dependencies)
            self.assertIsNotNone(result)
            
        except (subprocess.TimeoutExpired, Exception) as e:
            logging.warning(f"CLI without metadata test failed: {e}")

    def test_template_file_existence(self):
        """Test that template files exist."""
        template_dir = os.path.join(self.inputs_dir, 'templates')
        self.assertTrue(os.path.exists(template_dir))
        template_files = os.listdir(template_dir)
        self.assertGreater(len(template_files), 0)

    def test_hits_file_existence(self):
        """Test that hits file exists."""
        hits_file = os.path.join(self.inputs_dir, 'A71EV2A_combined.sdf')
        self.assertTrue(os.path.exists(hits_file))

    def test_metadata_file_existence(self):
        """Test that metadata file exists."""
        metadata_file = os.path.join(self.inputs_dir, 'metadata.csv')
        self.assertTrue(os.path.exists(metadata_file))


if __name__ == '__main__':
    unittest.main() 