########################################################################################################################

__doc__ = \
    """
    See GitHub documentation
    """
__author__ = "Kate Fieseler. [Github](https://github.com/kate-fie)"
__email__ = "kate.k.fieseler@gmail.com"
__date__ = "2025"
__license__ = "MIT"
__citation__ = ""

# Suppress fragmenstein warnings that appear during normal operation
import warnings
warnings.filterwarnings("ignore", message="PyRosetta is not installed. A mock object is loaded. Any Igor calls will fail.")
warnings.filterwarnings("ignore", category=RuntimeWarning, module="fragmenstein.igor.pyrosetta_import")
warnings.filterwarnings("ignore", message="pkg_resources is deprecated as an API")

from .utils.version import __version__

##################################################################################################################

if __name__ == "__main__":
    from .cli import main
    main()