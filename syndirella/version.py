__version__ = "2.0.0-alpha"

from typing import Dict

def get_versions() -> Dict[str, str]:
    """
    Return a dict of versions of os, python, fragmenstein etc.
    """
    import pkg_resources, sys, platform

    get_version = lambda name: pkg_resources.get_distribution(name).version

    return dict(python=sys.version,
                os_type=platform.system(),
                arc=platform.machine(),
                fragmenstein=get_version("fragmenstein"),
                pyrosetta=get_version("pyrosetta"),
                rdkit=get_version("rdkit"),
                rdkit_to_params=get_version('rdkit-to-params')
                )