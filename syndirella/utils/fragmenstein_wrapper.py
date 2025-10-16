"""
Wrapper for fragmenstein that suppresses warnings.
"""
import warnings

# Suppress fragmenstein warnings that appear during normal operation
warnings.filterwarnings("ignore", message="PyRosetta is not installed. A mock object is loaded. Any Igor calls will fail.")
warnings.filterwarnings("ignore", category=RuntimeWarning, module="fragmenstein.igor.pyrosetta_import")
warnings.filterwarnings("ignore", message="pkg_resources is deprecated as an API")

# Now import fragmenstein with warnings suppressed
from fragmenstein import Laboratory, Wictor, Igor
from fragmenstein.laboratory.validator import place_input_validator

# Re-export everything
__all__ = ['Laboratory', 'Wictor', 'Igor', 'place_input_validator']
