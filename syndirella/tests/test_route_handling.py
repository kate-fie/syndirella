import logging
import os
import unittest

from aizynth.AiZynthManager import AiZynthManager
from syndirella.Cobbler import Cobbler


# Configure logging at the application level
def setup_logging(log_level=logging.INFO, log_file=None):
    """Set up logging configuration."""
    # Create logger
    logger = logging.getLogger('syndirella')
    logger.setLevel(log_level)

    # Create formatter
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

    # Add console handler
    console_handler = logging.StreamHandler()
    console_handler.setLevel(log_level)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    # Add file handler if log_file is specified
    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(log_level)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    return logger


class TestAiZynthManagerSingle(unittest.TestCase):
    def setUp(self):
        self.logger = setup_logging(log_level=logging.INFO)

    def test_get_for_target(self, target_smiles: str = 'C1CCN(C(=O)NC2CCC2)CC1'):
        manager = AiZynthManager()
        manager.logger = self.logger
        output_path = f'{target_smiles}.json'
        manager.find_routes(target_smiles=target_smiles)
        manager.export_routes_to_dict(routes=manager.routes, output_path=output_path)
        self.assertTrue(os.path.exists(output_path), "json file was not created")


class TestCobbler(unittest.TestCase):
    def setUp(self):
        self.logger = setup_logging(log_level=logging.INFO)
        self.cobbler = Cobbler(
            scaffold_compound='C1CCN(C(=O)NC2CCC2)CC1',  # Fill in with an actual compound
            output_dir='syndirella/tests/outputs/test_route_handling',
            atom_diff_min=0,
            atom_diff_max=10,
            elab_single_reactant=False,
            retro_tool='aizynthfinder'
        )

    def test_get_routes(self):
        routes = self.cobbler.get_routes()
        assert len(routes) > 0


# TODO: Finish

class TestJustretroquery(unittest.TestCase):
    def setUp(self):
        self.logger = setup_logging(log_level=logging.INFO)

    def test_w_aizynth(self):
        pass

    def test_w_manifold(self):
        pass


if __name__ == '__main__':
    unittest.main()
