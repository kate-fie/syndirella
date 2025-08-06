# python
import logging
import os
import shutil
import unittest
import tempfile
import pandas as pd
from unittest.mock import patch, MagicMock, mock_open
from rdkit import Chem

from syndirella.pipeline import (
    run_pipeline, 
    PipelineConfig, 
    assert_scaffold_placement,
    elaborate_from_cobbler_workshops,
    start_elaboration,
    elaborate_compound_with_manual_routes,
    elaborate_compound_full_auto,
    process_row
)
from syndirella.constants import RetrosynthesisTool, DatabaseSearchTool


def handle_file_path(user_path: str) -> str:
    if os.path.isabs(user_path):
        return user_path
    return os.path.abspath(os.path.join(os.getcwd(), user_path))


class TestPipelineConfig(unittest.TestCase):
    """Test the PipelineConfig dataclass."""
    
    def setUp(self):
        self.settings = {
            'input': 'test_input.csv',
            'output': 'tests/outputs/test_pipeline/test_output/',
            'templates': 'tests/inputs/test_inputs/templates',
            'hits_path': 'tests/inputs/test_inputs/A71EV2A_combined.sdf',
            'metadata': 'tests/inputs/test_inputs/metadata.csv',
            'batch_num': 1,
            'atom_diff_min': 0,
            'atom_diff_max': 1,
            'scaffold_place_num': 1,
            'long_code_column': 'Long code',
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
    """Test individual pipeline functions with extensive mocking."""
    
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        self.inputs_dir = 'syndirella/tests/inputs/test_inputs'
        self.test_scaffold = "O=C(NC1CC(C(F)(F)F)C1)c1cc2ccsc2[nH]1"
        self.test_template_path = os.path.join(self.inputs_dir, "templates", "Ax0310a_apo-desolv.pdb")
        self.test_hits_path = os.path.join(self.inputs_dir, "A71EV2A_combined.sdf")
        self.test_output_dir = os.path.join(self.temp_dir, "test_output")
        self.metadata = os.path.join(self.inputs_dir, "metadata.csv")

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
            template_path=self.test_template_path,
            hits_path=self.test_hits_path,
            hits_names=['Ax0310a', 'Ax0528a'],
            output_dir=self.test_output_dir,
            scaffold_place_num=3
        )
        
        self.assertIsNotNone(result)
        mock_slipper_fitter.assert_called_once()
        mock_fitter.check_scaffold.assert_called_once()

    def test_assert_scaffold_placement_invalid_smiles(self):
        """Test scaffold placement assertion with invalid SMILES."""
        with self.assertRaises(ValueError):
            assert_scaffold_placement(
                scaffold="invalid_smiles",
                template_path=self.test_template_path,
                hits_path=self.test_hits_path,
                hits_names=['Ax0310a', 'Ax0528a'],
                output_dir=self.test_output_dir,
                scaffold_place_num=3
            )

    @patch('syndirella.pipeline.CobblersWorkshop')
    @patch('syndirella.pipeline.Slipper')
    def test_elaborate_from_cobbler_workshops(self, mock_slipper, mock_workshop):
        """Test elaboration from cobblers workshops."""
        mock_workshop_instance = MagicMock()
        mock_workshop.return_value = mock_workshop_instance
        mock_slipper_instance = MagicMock()
        mock_slipper.return_value = mock_slipper_instance
        
        # Mock the workshop methods
        mock_workshop_instance.get_products.return_value = pd.DataFrame({
            'smiles': ['CCO', 'CCCO'],
            'compound_set': ['test1', 'test2']
        })
        
        result = elaborate_from_cobbler_workshops(
            cobbler_workshops=[mock_workshop_instance],
            template_path=self.test_template_path,
            hits_path=self.test_hits_path,
            hits=['Ax0310a', 'Ax0528a'],
            batch_num=1,
            csv_path='test.csv',
            output_dir=self.test_output_dir,
            scaffold_placements={},
            additional_info=None
        )
        
        self.assertIsNotNone(result)
        mock_workshop.assert_called_once()
        mock_slipper.assert_called_once()

    @patch('syndirella.pipeline.assert_scaffold_placement')
    def test_start_elaboration(self, mock_assert_scaffold):
        """Test start elaboration function."""
        mock_assert_scaffold.return_value = {"test": "placement"}
        
        result = start_elaboration(
            product=self.test_scaffold,
            template_path=self.test_template_path,
            hits_path=self.test_hits_path,
            hits=['Ax0310a', 'Ax0528a'],
            output_dir=self.test_output_dir,
            scaffold_place_num=3,
            scaffold_place=True
        )
        
        self.assertIsNotNone(result)
        mock_assert_scaffold.assert_called_once()

    @patch('syndirella.pipeline.start_elaboration')
    def test_elaborate_compound_with_manual_routes(self, mock_start_elaboration):
        """Test elaboration with manual routes."""
        mock_start_elaboration.return_value = (1.0, {"test": "result"})
        
        result = elaborate_compound_with_manual_routes(
            product=self.test_scaffold,
            reactants=[('CCO', 'test_reactant')],
            reaction_names=['test_reaction'],
            num_steps=1,
            hits=['Ax0310a', 'Ax0528a'],
            template_path=self.test_template_path,
            hits_path=self.test_hits_path,
            batch_num=1,
            output_dir=self.test_output_dir,
            csv_path='test.csv',
            atom_diff_min=0,
            atom_diff_max=5,
            scaffold_place_num=3,
            only_scaffold_place=False,
            scaffold_place=True,
            elab_single_reactant=False,
            retro_tool=RetrosynthesisTool.MANIFOLD,
            db_search_tool=DatabaseSearchTool.ARTHOR,
            additional_info=None
        )
        
        self.assertIsNotNone(result)
        mock_start_elaboration.assert_called_once()

    @patch('syndirella.pipeline.Cobbler')
    @patch('syndirella.pipeline.start_elaboration')
    def test_elaborate_compound_full_auto(self, mock_start_elaboration, mock_cobbler):
        """Test full auto elaboration."""
        mock_cobbler_instance = MagicMock()
        mock_cobbler.return_value = mock_cobbler_instance
        mock_start_elaboration.return_value = (1.0, {"test": "result"})
        
        result = elaborate_compound_full_auto(
            product=self.test_scaffold,
            hits=['Ax0310a', 'Ax0528a'],
            template_path=self.test_template_path,
            hits_path=self.test_hits_path,
            batch_num=1,
            output_dir=self.test_output_dir,
            csv_path='test.csv',
            atom_diff_min=0,
            atom_diff_max=5,
            scaffold_place_num=3,
            only_scaffold_place=False,
            scaffold_place=True,
            elab_single_reactant=False,
            retro_tool=RetrosynthesisTool.MANIFOLD,
            db_search_tool=DatabaseSearchTool.ARTHOR,
            additional_info=None
        )
        
        self.assertIsNotNone(result)
        mock_cobbler.assert_called_once()
        mock_start_elaboration.assert_called_once()

    @patch('syndirella.pipeline.check_inputs.get_template_path')
    @patch('syndirella.pipeline.elaborate_compound_with_manual_routes')
    def test_process_row_manual_routes(self, mock_elaborate_manual, mock_get_template):
        """Test processing row with manual routes."""
        mock_get_template.return_value = self.test_template_path
        mock_elaborate_manual.return_value = (1.0, {"test": "result"})
        
        row = pd.Series({
            'smiles': self.test_scaffold,
            'template': 'Ax0310a'
        })
        
        config = PipelineConfig.from_settings({
            'input': 'test.csv',
            'output': self.test_output_dir,
            'templates': self.inputs_dir + '/templates',
            'hits_path': self.test_hits_path,
            'metadata': self.metadata,
            'batch_num': 1,
            'atom_diff_min': 0,
            'atom_diff_max': 5,
            'scaffold_place_num': 3,
            'long_code_column': 'Long code',
            'retro_tool': 'manifold',
            'db_search_tool': 'arthor',
            'manual': True
        })
        
        result = process_row(row, config)
        
        self.assertIsNotNone(result)
        mock_get_template.assert_called_once()
        mock_elaborate_manual.assert_called_once()

    @patch('syndirella.pipeline.elaborate_compound_full_auto')
    @patch('syndirella.pipeline.check_inputs.get_template_path')
    def test_process_row_auto_routes(self, mock_get_template_path, mock_elaborate_auto):
        """Test processing row with auto routes."""
        mock_get_template_path.return_value = self.test_template_path
        mock_elaborate_auto.return_value = (1.0, {"test": "result"})
        
        row = pd.Series({
            'smiles': self.test_scaffold,
            'template': 'Ax0310a'
        })
        
        config = PipelineConfig.from_settings({
            'input': 'test.csv',
            'output': self.test_output_dir,
            'templates': self.inputs_dir + '/templates',
            'hits_path': self.test_hits_path,
            'metadata': self.metadata,
            'batch_num': 1,
            'atom_diff_min': 0,
            'atom_diff_max': 5,
            'scaffold_place_num': 3,
            'long_code_column': 'Long code',
            'retro_tool': 'manifold',
            'db_search_tool': 'arthor',
            'manual': False
        })
        
        result = process_row(row, config)
        
        self.assertIsNotNone(result)
        mock_get_template_path.assert_called_once()
        mock_elaborate_auto.assert_called_once()


class TestPipelineIntegration(unittest.TestCase):
    """Test pipeline integration with extensive mocking."""
    
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        self.inputs_dir = 'syndirella/tests/inputs/test_inputs'
        self.test_output_dir = os.path.join(self.temp_dir, "test_output")
        
        # Create test input CSV
        self.test_input_csv = os.path.join(self.temp_dir, "test_input.csv")
        test_data = pd.DataFrame({
            'smiles': ['O=C(NC1CC(C(F)(F)F)C1)c1cc2ccsc2[nH]1'],
            'hit1': ['Ax0310a'],
            'hit2': ['Ax0528a'],
            'hit3': [''],
            'template': ['Ax0310a'],
            'compound_set': ['test_set']
        })
        test_data.to_csv(self.test_input_csv, index=False)

    def tearDown(self):
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    @patch('syndirella.pipeline.Cobbler')
    @patch('syndirella.pipeline.SlipperFitter')
    @patch('syndirella.pipeline.pd.read_csv')
    def test_pipeline_with_only_scaffold_place(self, mock_read_csv, mock_slipper_fitter, mock_cobbler):
        """Test pipeline with only scaffold placement."""
        # Mock the CSV reading
        mock_read_csv.return_value = pd.DataFrame({
            'smiles': ['O=C(NC1CC(C(F)(F)F)C1)c1cc2ccsc2[nH]1'],
            'hit1': ['Ax0310a'],
            'hit2': ['Ax0528a'],
            'hit3': [''],
            'template': ['Ax0310a'],
            'compound_set': ['test_set']
        })
        
        # Mock SlipperFitter
        mock_fitter = MagicMock()
        mock_slipper_fitter.return_value = mock_fitter
        mock_fitter.check_scaffold.return_value = "test_placement"
        
        settings = {
            'input': self.test_input_csv,
            'output': self.test_output_dir,
            'templates': self.inputs_dir + '/templates',
            'hits_path': self.inputs_dir + '/A71EV2A_combined.sdf',
            'metadata': self.inputs_dir + '/metadata.csv',
            'only_scaffold_place': True,
            'scaffold_place_num': 3,
            'retro_tool': 'manifold',
            'db_search_tool': 'arthor'
        }
        
        result = run_pipeline(settings)
        
        self.assertIsNotNone(result)
        mock_slipper_fitter.assert_called()

    @patch('syndirella.pipeline.Cobbler')
    @patch('syndirella.pipeline.SlipperFitter')
    @patch('syndirella.pipeline.pd.read_csv')
    def test_pipeline_with_custom_parameters(self, mock_read_csv, mock_slipper_fitter, mock_cobbler):
        """Test pipeline with custom parameters."""
        # Mock the CSV reading
        mock_read_csv.return_value = pd.DataFrame({
            'smiles': ['O=C(NC1CC(C(F)(F)F)C1)c1cc2ccsc2[nH]1'],
            'hit1': ['Ax0310a'],
            'hit2': ['Ax0528a'],
            'hit3': [''],
            'template': ['Ax0310a'],
            'compound_set': ['test_set']
        })
        
        # Mock SlipperFitter
        mock_fitter = MagicMock()
        mock_slipper_fitter.return_value = mock_fitter
        mock_fitter.check_scaffold.return_value = "test_placement"
        
        # Mock Cobbler
        mock_cobbler_instance = MagicMock()
        mock_cobbler.return_value = mock_cobbler_instance
        mock_cobbler_instance.get_products.return_value = pd.DataFrame({
            'smiles': ['CCO'],
            'compound_set': ['test1']
        })
        
        settings = {
            'input': self.test_input_csv,
            'output': self.test_output_dir,
            'templates': self.inputs_dir + '/templates',
            'hits_path': self.inputs_dir + '/A71EV2A_combined.sdf',
            'metadata': self.inputs_dir + '/metadata.csv',
            'atom_diff_min': 0,
            'atom_diff_max': 5,
            'scaffold_place_num': 3,
            'retro_tool': 'manifold',
            'db_search_tool': 'arthor'
        }
        
        result = run_pipeline(settings)
        
        self.assertIsNotNone(result)
        mock_slipper_fitter.assert_called()
        mock_cobbler.assert_called()


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
            'metadata': 'test/metadata.csv'
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