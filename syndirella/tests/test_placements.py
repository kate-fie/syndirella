# python
import logging
import os
import shutil
import unittest
import tempfile
import pandas as pd
from unittest.mock import patch, MagicMock, mock_open
from rdkit import Chem

from syndirella.slipper.SlipperFitter import SlipperFitter
from syndirella.utils.error import *


def handle_file_path(user_path: str) -> str:
    if os.path.isabs(user_path):
        return user_path
    return os.path.abspath(os.path.join(os.getcwd(), user_path))


class TestSlipperFitterInitialization(unittest.TestCase):
    """Test SlipperFitter initialization and basic setup."""
    
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        self.template_path = 'syndirella/tests/inputs/test_inputs/templates/Ax0310a_apo-desolv.pdb'
        self.hits_path = 'syndirella/tests/inputs/test_inputs/A71EV2A_combined.sdf'
        self.hits_names = ['A71EV2A-x0526_A_147_1_A71EV2A-x0526+A+147+1', 'A71EV2A-x0540_A_147_1_A71EV2A-x0540+A+147+1']
        self.output_dir = self.temp_dir
        self.test_id = 'test_id'
        self.route_uuid = 'test_uuid'

    def tearDown(self):
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    def test_slipper_fitter_initialization(self):
        """Test basic SlipperFitter initialization."""
        slipper_fitter = SlipperFitter(
            template_path=self.template_path,
            hits_path=self.hits_path,
            hits_names=self.hits_names,
            output_dir=self.output_dir,
            id=self.test_id,
            route_uuid=self.route_uuid
        )
        
        self.assertEqual(slipper_fitter.template_path, self.template_path)
        self.assertEqual(slipper_fitter.hits_path, self.hits_path)
        self.assertEqual(slipper_fitter.hits_names, self.hits_names)
        self.assertEqual(slipper_fitter.output_dir, self.output_dir)
        self.assertEqual(slipper_fitter.id, self.test_id)
        self.assertEqual(slipper_fitter.route_uuid, self.route_uuid)
        self.assertEqual(slipper_fitter.atom_diff_min, 0)
        self.assertEqual(slipper_fitter.atom_diff_max, 10)
        self.assertEqual(slipper_fitter.n_cores, 8)
        self.assertEqual(slipper_fitter.timeout, 240)

    def test_slipper_fitter_initialization_with_scaffold_placements(self):
        """Test SlipperFitter initialization with scaffold placements."""
        scaffold_placements = {Chem.MolFromSmiles('CCO'): 'test_placement'}
        
        slipper_fitter = SlipperFitter(
            template_path=self.template_path,
            hits_path=self.hits_path,
            hits_names=self.hits_names,
            output_dir=self.output_dir,
            scaffold_placements=scaffold_placements
        )
        
        self.assertEqual(slipper_fitter.scaffold_placements, scaffold_placements)

    def test_slipper_fitter_initialization_defaults(self):
        """Test SlipperFitter initialization with default values."""
        slipper_fitter = SlipperFitter(
            template_path=self.template_path,
            hits_path=self.hits_path,
            hits_names=self.hits_names,
            output_dir=self.output_dir
        )
        
        self.assertIsNone(slipper_fitter.id)
        self.assertIsNone(slipper_fitter.route_uuid)
        self.assertIsNone(slipper_fitter.scaffold_placements)
        self.assertIsNone(slipper_fitter.final_products)
        self.assertIsNone(slipper_fitter.placements)
        self.assertIsNone(slipper_fitter.merged_placements)


class TestSlipperFitterHelperMethods(unittest.TestCase):
    """Test SlipperFitter helper methods."""
    
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        self.template_path = 'syndirella/tests/inputs/test_inputs/templates/Ax0310a_apo-desolv.pdb'
        self.hits_path = 'syndirella/tests/inputs/test_inputs/A71EV2A_combined.sdf'
        self.hits_names = ['A71EV2A-x0526_A_147_1_A71EV2A-x0526+A+147+1', 'A71EV2A-x0540_A_147_1_A71EV2A-x0540+A+147+1']
        self.output_dir = self.temp_dir
        
        self.slipper_fitter = SlipperFitter(
            template_path=self.template_path,
            hits_path=self.hits_path,
            hits_names=self.hits_names,
            output_dir=self.output_dir
        )

    def tearDown(self):
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    def test_print_diff(self):
        """Test _print_diff method."""
        orig_df = pd.DataFrame({'col1': [1, 2], 'col2': [3, 4]})
        input_df = pd.DataFrame({'col1': [1, 2], 'col2': [3, 5]})
        
        # This method should not raise any exceptions
        self.slipper_fitter._print_diff(orig_df, input_df, "test")

    @patch('syndirella.utils.fairy.generate_inchi_ID')
    def test_add_hits(self, mock_generate_id):
        """Test add_hits method."""
        mock_generate_id.return_value = 'test_id'
        
        input_df = pd.DataFrame({
            'compound_set': ['test1', 'test2'],
            'Long code': ['code1', 'code2'],
            'SMILES': ['CCO', 'CCCO']
        })
        
        result = self.slipper_fitter.add_hits(input_df)
        
        self.assertIsInstance(result, pd.DataFrame)
        self.assertIn('hits', result.columns)
        self.assertTrue(all(result['hits'].values))

    def test_prep_scaffold_input_df(self):
        """Test _prep_scaffold_input_df method."""
        scaffold = Chem.MolFromSmiles('O=C(NC1CC(C(F)(F)F)C1)c1cc2ccsc2[nH]1')
        scaffold_name = 'test_scaffold'
        
        result = self.slipper_fitter._prep_scaffold_input_df(scaffold, scaffold_name)
        
        self.assertIsInstance(result, pd.DataFrame)
        self.assertIn('hits', result.columns)
        self.assertTrue(all(result['hits'].values))


class TestSlipperFitterScaffoldChecking(unittest.TestCase):
    """Test SlipperFitter scaffold checking functionality."""
    
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        self.template_path = 'syndirella/tests/inputs/test_inputs/templates/Ax0310a_apo-desolv.pdb'
        self.hits_path = 'syndirella/tests/inputs/test_inputs/A71EV2A_combined.sdf'
        self.hits_names = ['A71EV2A-x0526_A_147_1_A71EV2A-x0526+A+147+1', 'A71EV2A-x0540_A_147_1_A71EV2A-x0540+A+147+1']
        self.output_dir = self.temp_dir
        
        self.slipper_fitter = SlipperFitter(
            template_path=self.template_path,
            hits_path=self.hits_path,
            hits_names=self.hits_names,
            output_dir=self.output_dir
        )

    def tearDown(self):
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    @patch.object(SlipperFitter, 'setup_Fragmenstein')
    @patch.object(SlipperFitter, '_place_scaffold')
    @patch.object(SlipperFitter, '_prep_scaffold_input_df')
    @patch('syndirella.utils.fairy.generate_inchi_ID')
    @patch('os.path.exists')
    def test_check_scaffold_success(self, mock_exists, mock_generate_id, mock_prep_df, mock_place, mock_setup):
        """Test successful scaffold checking."""
        # Mock setup
        mock_generate_id.return_value = 'test_id'
        mock_prep_df.return_value = pd.DataFrame({'hits': [True]})
        mock_lab = MagicMock()
        mock_setup.return_value = mock_lab
        mock_place.return_value = Chem.MolFromSmiles('O=C(NC1CC(C(F)(F)F)C1)c1cc2ccsc2[nH]1')
        mock_exists.return_value = True
        
        scaffold = Chem.MolFromSmiles('O=C(NC1CC(C(F)(F)F)C1)c1cc2ccsc2[nH]1')
        scaffold_name = 'test_scaffold'
        scaffold_place_num = 3
        
        result = self.slipper_fitter.check_scaffold(scaffold, scaffold_name, scaffold_place_num)
        
        self.assertIsNotNone(result)
        mock_setup.assert_called_once()
        mock_place.assert_called_once()

    @patch.object(SlipperFitter, 'setup_Fragmenstein')
    @patch.object(SlipperFitter, '_place_scaffold')
    @patch.object(SlipperFitter, '_prep_scaffold_input_df')
    @patch('syndirella.utils.fairy.generate_inchi_ID')
    @patch('os.path.exists')
    def test_check_scaffold_failure(self, mock_exists, mock_generate_id, mock_prep_df, mock_place, mock_setup):
        """Test failed scaffold checking."""
        # Mock setup
        mock_generate_id.return_value = 'test_id'
        mock_prep_df.return_value = pd.DataFrame({'hits': [True]})
        mock_lab = MagicMock()
        mock_setup.return_value = mock_lab
        mock_place.return_value = None  # Placement failed
        mock_exists.return_value = False
        
        scaffold = Chem.MolFromSmiles('O=C(NC1CC(C(F)(F)F)C1)c1cc2ccsc2[nH]1')
        scaffold_name = 'test_scaffold'
        scaffold_place_num = 1  # Only one attempt to avoid multiple calls
        
        result = self.slipper_fitter.check_scaffold(scaffold, scaffold_name, scaffold_place_num)
        
        self.assertIsNone(result)
        mock_setup.assert_called_once()
        mock_place.assert_called_once()

    @patch.object(SlipperFitter, 'setup_Fragmenstein')
    @patch.object(SlipperFitter, '_place_scaffold')
    @patch.object(SlipperFitter, '_prep_scaffold_input_df')
    @patch('syndirella.utils.fairy.generate_inchi_ID')
    def test_check_scaffold_invalid_input_df(self, mock_generate_id, mock_prep_df, mock_place, mock_setup):
        """Test scaffold checking with invalid input dataframe."""
        # Mock setup with invalid dataframe (no hits)
        mock_generate_id.return_value = 'test_id'
        mock_prep_df.return_value = pd.DataFrame({'hits': [False]})
        
        scaffold = Chem.MolFromSmiles('O=C(NC1CC(C(F)(F)F)C1)c1cc2ccsc2[nH]1')
        scaffold_name = 'test_scaffold'
        scaffold_place_num = 3
        
        with self.assertRaises(AssertionError):
            self.slipper_fitter.check_scaffold(scaffold, scaffold_name, scaffold_place_num)


class TestSlipperFitterPlacementMethods(unittest.TestCase):
    """Test SlipperFitter placement-related methods."""
    
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        self.template_path = 'syndirella/tests/inputs/test_inputs/templates/Ax0310a_apo-desolv.pdb'
        self.hits_path = 'syndirella/tests/inputs/test_inputs/A71EV2A_combined.sdf'
        self.hits_names = ['A71EV2A-x0526_A_147_1_A71EV2A-x0526+A+147+1', 'A71EV2A-x0540_A_147_1_A71EV2A-x0540+A+147+1']
        self.output_dir = self.temp_dir
        
        self.slipper_fitter = SlipperFitter(
            template_path=self.template_path,
            hits_path=self.hits_path,
            hits_names=self.hits_names,
            output_dir=self.output_dir
        )

    def tearDown(self):
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    @patch.object(SlipperFitter, 'setup_Fragmenstein')
    @patch.object(SlipperFitter, '_place_scaffold')
    def test_place_scaffold(self, mock_place, mock_setup):
        """Test _place_scaffold method."""
        mock_lab = MagicMock()
        mock_setup.return_value = mock_lab
        mock_place.return_value = Chem.MolFromSmiles('O=C(NC1CC(C(F)(F)F)C1)c1cc2ccsc2[nH]1')
        
        input_df = pd.DataFrame({
            'compound_set': ['test1'],
            'Long code': ['code1'],
            'SMILES': ['CCO'],
            'hits': [True]
        })
        
        result = self.slipper_fitter._place_scaffold(mock_lab, input_df)
        
        self.assertIsInstance(result, Chem.Mol)

    def test_check_intra_geom_flatness_results(self):
        """Test _check_intra_geom_flatness_results method."""
        geometries = {
            'results': {
                'bond_lengths_within_bounds': True,
                'bond_angles_within_bounds': True,
                'no_internal_clash': True
            }
        }
        flat_results = {
            'results': {
                'flatness_passes': True
            }
        }
        
        result = self.slipper_fitter._check_intra_geom_flatness_results(geometries, flat_results)
        
        self.assertIsInstance(result, bool)

    @patch('builtins.open', new_callable=mock_open, read_data="test content")
    def test_setup_Fragmenstein(self, mock_file):
        """Test setup_Fragmenstein method."""
        output_path = os.path.join(self.temp_dir, 'test_output')
        
        result = self.slipper_fitter.setup_Fragmenstein(output_path)
        
        # The result should be a Laboratory object
        from fragmenstein import Laboratory
        self.assertIsInstance(result, Laboratory)


class TestSlipperFitterFormattingMethods(unittest.TestCase):
    """Test SlipperFitter formatting and output methods."""
    
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        self.template_path = 'syndirella/tests/inputs/test_inputs/templates/Ax0310a_apo-desolv.pdb'
        self.hits_path = 'syndirella/tests/inputs/test_inputs/A71EV2A_combined.sdf'
        self.hits_names = ['A71EV2A-x0526_A_147_1_A71EV2A-x0526+A+147+1', 'A71EV2A-x0540_A_147_1_A71EV2A-x0540+A+147+1']
        self.output_dir = self.temp_dir
        
        self.slipper_fitter = SlipperFitter(
            template_path=self.template_path,
            hits_path=self.hits_path,
            hits_names=self.hits_names,
            output_dir=self.output_dir
        )

    def tearDown(self):
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    def test_fix_intxns(self):
        """Test fix_intxns method."""
        # Create a mock placements dataframe
        self.slipper_fitter.placements = pd.DataFrame({
            'compound_set': ['test1'],
            'Long code': ['code1'],
            'SMILES': ['CCO'],
            'hits': [True]
        })
        
        # This method should not raise any exceptions
        self.slipper_fitter.fix_intxns()

    def test_merge_placements(self):
        """Test merge_placements method."""
        # Create mock placements dataframes
        self.slipper_fitter.placements = pd.DataFrame({
            'compound_set': ['test1'],
            'Long code': ['code1'],
            'SMILES': ['CCO'],
            'hits': [True],
            'name': ['test1']
        })
        
        self.slipper_fitter.final_products = pd.DataFrame({
            'compound_set': ['test1'],
            'Long code': ['code1'],
            'SMILES': ['CCO'],
            'name': ['test1']
        })
        
        result = self.slipper_fitter.merge_placements()
        
        self.assertIsInstance(result, pd.DataFrame)


class TestSlipperFitterIntegration(unittest.TestCase):
    """Test SlipperFitter integration functionality."""
    
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        self.template_path = 'syndirella/tests/inputs/test_inputs/templates/Ax0310a_apo-desolv.pdb'
        self.hits_path = 'syndirella/tests/inputs/test_inputs/A71EV2A_combined.sdf'
        self.hits_names = ['A71EV2A-x0526_A_147_1_A71EV2A-x0526+A+147+1', 'A71EV2A-x0540_A_147_1_A71EV2A-x0540+A+147+1']
        self.output_dir = self.temp_dir
        
        self.slipper_fitter = SlipperFitter(
            template_path=self.template_path,
            hits_path=self.hits_path,
            hits_names=self.hits_names,
            output_dir=self.output_dir
        )

    def tearDown(self):
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    @patch.object(SlipperFitter, 'prep_products')
    @patch.object(SlipperFitter, 'place_products')
    @patch.object(SlipperFitter, 'check_intra_geometry')
    @patch.object(SlipperFitter, 'format_placements')
    @patch.object(SlipperFitter, 'merge_placements')
    def test_fit_products(self, mock_merge, mock_format, mock_check, mock_place, mock_prep):
        """Test the main fit_products method."""
        # Mock all the component methods
        mock_prep.return_value = pd.DataFrame({'test': ['data']})
        mock_place.return_value = pd.DataFrame({'test': ['data']})
        mock_check.return_value = pd.DataFrame({'test': ['data']})
        mock_merge.return_value = pd.DataFrame({'test': ['data']})
        
        result = self.slipper_fitter.fit_products()
        
        self.assertIsInstance(result, pd.DataFrame)
        mock_prep.assert_called_once()
        mock_place.assert_called_once()
        mock_check.assert_called_once()
        mock_format.assert_called_once()
        mock_merge.assert_called_once()

    @patch('os.path.exists')
    @patch('builtins.open', new_callable=mock_open, read_data='{"mRMSD": 1.5, "∆∆G": -2.0, "comRMSD": 1.0, "∆G": -1.5}')
    def test_get_scaffold_check_values(self, mock_file, mock_exists):
        """Test get_scaffold_check_values method."""
        mock_exists.return_value = True
        
        scaffold_path = os.path.join(self.temp_dir, 'test_scaffold.mol')
        
        result = self.slipper_fitter.get_scaffold_check_values(scaffold_path)
        
        self.assertIsInstance(result, tuple)
        self.assertEqual(len(result), 4)


if __name__ == '__main__':
    unittest.main()
