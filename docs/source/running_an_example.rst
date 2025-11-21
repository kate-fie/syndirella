==================
Running an example
==================

This guide will walk you through an example to run Syndirella's pipeline through a jupyter notebook.

.. attention::

   If you have a Windows OS, use Docker Desktop by following these steps before continuing:

   1. Create a Docker account.
   2. Install Docker Desktop.
   3. In Docker Desktop: go to Add Extensions.
   4. Search for Jupyter Notebook Scientific Python Stack Extension and install it.
   5. Open newly installed extension.
   6. Open a terminal in the JupyterLab interface.

1. Install Syndirella.

.. code-block:: bash

   conda create -n syndirella python=3.10
   conda activate syndirella
   pip install "cython<3.2"  # Required to avoid Cython 3.2.0 compiler bug when building cgrtools
   pip install syndirella

**Note:** If `cgrtools` installation fails, try running:

.. code-block:: bash

    conda install -c conda-forge c-compiler cxx-compiler
    pip install --no-build-isolation cgrtools

**Troubleshooting:** If you encounter a ``TypeError: 'AttributeFilledMock' object is not iterable`` error, this is related to a PyRosetta dependency by Fragmenstein. To resolve this, you can download PyRosetta for academic and non-commercial use (see `PyRosetta License <https://github.com/RosettaCommons/rosetta/blob/main/LICENSE.PyRosetta.md>`_). 

.. code-block:: bash

    pip install pyrosetta-installer
    python -c 'import pyrosetta_installer; pyrosetta_installer.install_pyrosetta()'

2. Clone Syndirella repo.

.. code-block:: bash

   git clone https://github.com/kate-fie/syndirella.git
   cd syndirella


3. Open notebook in the syndirella environment at `syndirella/example/run_syndirella_example/run_examples.ipynb` which will take you through installation,
looking at designs, running a handful of placements, and analysing outputs.




