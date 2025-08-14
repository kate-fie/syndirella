==========
User Guide
==========

Installation and setup
======================

1. Install Syndirella:
----------------------

On Mac OS and Linux you can install from PyPI using Conda.

.. code-block:: bash

    conda create -n syndirella python=3.10
    conda activate syndirella
    pip install --upgrade syndirella

.. attention::

    Installation and usage have not been tested on Windows OS.

2. Install AiZynthFinder:
-------------------------------------

Syndirella can use AiZynthFinder for retrosynthesis functionality (see Retrosynthesis for more information).

**Automatic Setup (Recommended):**
Syndirella will automatically handle AiZynthFinder setup when you first use it. Simply run:

.. code-block:: bash

    syndirella setup-aizynth

.. attention::

    ⚠️ This will download large model files (~750MB total) if run the first time. These are required to run AiZynthFinder. This will automatically download the required data to `[syndirella_package_path]/aizynth` directory and create the configuration file automatically. 

**Manual Setup (Alternative):**
If you prefer manual setup:

.. code-block:: bash

    cd [syndirella_package_path]/aizynth
    download_public_data .
    # Update config.yml (if you prefer)
    export AIZYNTH_CONFIG_FILE="path/to/syndirella/aizynth/config.yml"

Command Line Interface
======================

Syndirella provides a command-line interface with multiple subcommands. Get help with `-h` or `--help`.

**Available Commands:**

- **setup-aizynth**: Setup AiZynthFinder data and configuration
- **run**: Run the main Syndirella pipeline  
- **add-reaction**: Add a new reaction to the library

**Main Help Output:**

.. code-block:: bash

    usage: syndirella [-h] {setup-aizynth,run,add-reaction} ...

    Run the Syndirella pipeline with specified configurations.

    Available commands:
      setup-aizynth    Setup AiZynthFinder data and configuration
      run              Run the main Syndirella pipeline  
      add-reaction     Add a new reaction to the library

    Syndirella is installed at [path_to_installation]

**Run Command Help:**

.. code-block:: bash

    usage: syndirella run [-h] -i INPUT -o OUTPUT [-t TEMPLATES] [--hits_path HITS_PATH] [--metadata METADATA] [--products PRODUCTS] [--batch_num BATCH_NUM] [--manual] [--only_scaffold_place] [--scaffold_place_num SCAFFOLD_PLACE_NUM] [--retro_tool {manifold,aizynthfinder}] [--db_search_tool {postera,arthor}] [--profile] [--atom_diff_min ATOM_DIFF_MIN] [--atom_diff_max ATOM_DIFF_MAX] [--long_code_column LONG_CODE_COLUMN] [--just_retro] [--no_scaffold_place] [--elab_single_reactant]

    Run the full Syndirella pipeline with specified input files and parameters.

    options:
      -h, --help            show this help message and exit
      -i INPUT, --input INPUT
                            Input .csv file path for the pipeline.
      -o OUTPUT, --output OUTPUT
                            Output directory for the pipeline results.
      -t TEMPLATES, --templates TEMPLATES
                            Absolute path to a directory containing the template(s).
      --hits_path HITS_PATH
                            Absolute path to hits_path for placements (.sdf or .mol).
      --metadata METADATA   Absolute path to metadata for placements.
      --products PRODUCTS   Absolute path to products for placements.
      --batch_num BATCH_NUM
                            Batch number for processing. (default: 10000)
      --manual              Use manual routes for processing. (default: False)
      --only_scaffold_place
                            Only place scaffolds. Do not continue to elaborate. (default: False)
      --scaffold_place_num SCAFFOLD_PLACE_NUM
                            Number of times to attempt scaffold placement. (default: 5)
      --retro_tool {manifold,aizynthfinder}
                            Retrosynthesis tool to use. (default: manifold)
      --db_search_tool {postera,arthor}
                            Database search tool to use. (default: postera)
      --profile             Run the pipeline with profiling. (default: False)
      --atom_diff_min ATOM_DIFF_MIN
                            Minimum atom difference between elaborations and scaffold to keep. (default: 0)
      --atom_diff_max ATOM_DIFF_MAX
                            Maximum atom difference between elaborations and scaffold to keep. (default: 10)
      --long_code_column LONG_CODE_COLUMN
                            Column name for long code in metadata csv to match to SDF name. (default: Long code)
      --just_retro          Only run retrosynthesis querying of scaffolds. (default: False)
      --no_scaffold_place   Do not place scaffolds initially before elaborating. (default: False)
      --elab_single_reactant
                            Only elaborate one reactant per elaboration series. (default: False)

**Add Reaction Command Help:**

.. code-block:: bash

    usage: syndirella add-reaction [-h] --name NAME --smirks SMIRKS [--find_parent] [--fp_type {maccs_rxn_fp,morgan_rxn_fp}] [--threshold THRESHOLD] [--similarity_metric {tanimoto,dice,cosine}]

    Add a new reaction SMIRKS to the reaction library with optional parent finding.

    options:
      -h, --help            show this help message and exit
      --name NAME           Name of the new reaction.
      --smirks SMIRKS       SMIRKS string for the reaction.
      --find_parent         If True, treat as a child reaction and find parent based on similarity. (default: False)
      --fp_type {maccs_rxn_fp,morgan_rxn_fp}
                            Fingerprint type for similarity calculation. (default: maccs_rxn_fp)
      --threshold THRESHOLD
                            Similarity threshold for finding parent reaction. (default: 0.2)
      --similarity_metric {tanimoto,dice,cosine}
                            Similarity metric for finding parent reaction. (default: tanimoto)

**Setup AiZynthFinder Command Help:**

.. code-block:: bash

    usage: syndirella setup-aizynth [-h]

    Automatically download AiZynthFinder data and create configuration file.

    options:
      -h, --help  show this help message and exit

Basic Usage
===========

Elaborate a set of scaffolds using these steps:

1. Setup fragments and protein templates
----------------------------------------

Download the fragment hits from Fragalysis. In the download folder the important files are:

::

    metadata.csv # contains short and long codes
    target_name_combined.sdf # fragment poses with long code names
    /aligned_files/fragment_name/fragment_name_apo-desolv.pdb # apo pdb used for placement

2. Create input csv
-------------------

Syndirella can be run either in *automatic* or *manual* mode.

**Automatic**:
    Scaffolds can be elaborated by routes automatically proposed by Manifold.
    An example template is at ``syndirella/syndirella_input_template.csv``.

Required headers:

    ``smiles``:
        smiles string of :term:`scaffold`.
    ``hit1``:
        string of short code (or any substring of short code found in the metadata.csv) of 1 fragment inspiring hit.
    ``template``:
        path to apo protein template to use for :term:`placement`.
    ``compound_set``:
        string or int identifier.

Not required headers:

    ``hitX``:
        string of short code of additional fragment inspiring hit.

.. note::

    Any number of fragment inspirations can be used. You just need to specify in a seperate header. Ex.
    ``hit1, hit2, hit3, hit4, hit5``.


**Manual**:
    You can set the exact route to elaborate the scaffold with the reaction names, exact reactants, and number of steps in the route.
    An example template is at ``syndirella/syndirella_input_template_manual.csv``.

Required headers:

    ``smiles``:
        smiles string of scaffold.
    ``hit1``:
        string of short code (or any substring of short code found in the metadata.csv) of 1 fragment inspiring hit.
    ``template``:
        path to apo protein template to use for :term:`placement`.
    ``compound_set``:
        string or int identifier.
    ``reaction_name_step1``:
        string of reaction name.
    ``reactant_step1``:
        smiles string of reactant.

Not required headers:

    ``reactant2_step1``:
        smiles string of second reactant in reaction step 1.
    ``product_stepX``:
        smiles string of product of step X. Only required for internal or first step to specify reactant for next step. Not required
        if step is final step of route (as the scaffold is the final product).
    ``reaction_name_stepX``:
        string of reaction name of step X.
    ``reactant_stepX``:
        smiles string of reactant that is *not* a product of previous step.
    ``hitX``:
        string of short code of additional fragment inspiring hit. Any number of hits can be used.

3. Run!
-------

Run pipeline in *automatic* mode:

.. code-block:: bash

    syndirella run --input [path_to_automatic.csv] --output [path_to_output_dir] --templates [path_to_templates_dir]
    --hits_path [path_to_fragments.sdf] --metadata [path_to_metadata.csv]


Run pipeline in *manual* mode:
    Add ``--manual`` flag.

4. Outputs
----------

**Output directory structure:**

🔑🔑🔑: Inchi key of scaffold. Example: ``ZJENMQHSGLZNHL-UHFFFAOYSA-N``

.. code-block::

    output_dir
    ├── 🔑🔑🔑-scaffold-check # scaffold check directory per scaffold
    │   └── scaffold-check
    │       ├── scaffold-check.holo_minimised.pdb
    │       ├── scaffold-check.minimised.json
    │       └── scaffold-check.minimised.mol
    ├── 🔑🔑🔑 # directory per scaffold
    │   ├── extra
    │   │   ├── 🔑🔑🔑_[route_uuid]_[rxn_name]_r[reactant_num]_[step_num]of[total_steps].pkl.gz # reactants for step
    │   │   └── continued for all steps...
    │   ├── output
    │   │   ├── 🔑🔑🔑_[route_uuid]_[num]-[stereoisomer]
    │   │   │   ├── 🔑🔑🔑_[route_uuid]_[num]-[stereoisomer].mol
    │   │   │   ├── 🔑🔑🔑_[route_uuid]_[num]-[stereoisomer].json # energy values
    │   │   └── continued for all products...
    │   ├── 🔑🔑🔑_[route_uuid]_[rxn_name]_products_[last_step]of[total_steps].pkl.gz & .csv # final products
    │   ├── 🔑🔑🔑_[route_uuid]_[rxn_name]_products_[last_step]of[total_steps]_placements.pkl.gz & .csv # merged placements with products info
    │   ├── 🔑🔑🔑_[route_uuid]_fragmenstein_placements.pkl.gz & .csv # fragmenstein output
    │   └── 🔑🔑🔑_[route_uuid]_to_hippo.pkl.gz # full routes and placements
    ├── continued for all scaffolds...
    └── [input_csv]_output_YYYYMMDD_HHMM.csv # summary stats of all scaffolds

**Important output files:**

**[input_csv]_output_YYYYMMDD_HHMM.csv:**
    Summary stats of all scaffolds. Most columns are self-explanatory. The following columns might need clarification:

    ``total_num_products_enumstereo``:
        Total number of products enumerated with stereochemistry in the final step. This is counting the number of unique
        products with stereochemistry, so if a product with same stereochemistry is generated multiple times via different routes
        it will only be counted once.

    ``total_num_unique_products``:
        Total number of unique products without stereochemistry in the final step. If a product is generated multiple times
        by different routes it will only be counted once.

**🔑🔑🔑_[route_uuid]_[rxn_name]_products_[last_step]of[total_steps]_placements.pkl.gz & .csv:**
    Merged placements with products info.

**🔑🔑🔑_[route_uuid]_to_hippo.pkl.gz:**
    Full routes and placements.

.. note::

    Placements of products are labeled succesful if:
        1. ΔΔG < 0.
        2. comRMSD < 2.0 Å.
        3. Pose of product passes `PoseBusters <https://github.com/maabuu/posebusters>`_ intrageometry checks:
            - Bond lengths: The bond lengths in the input molecule are within 0.75 of the lower and 1.25 of the upper bounds determined by distance geometry.
            - Bond angles: The angles in the input molecule are within 0.75 of the lower and 1.25 of the upper bounds determined by distance geometry.
            - Planar aromatic rings: All atoms in aromatic rings with 5 or 6 members are within 0.25 Å of the closest shared plane.
            - Planar double bonds: The two carbons of aliphatic carbon–carbon double bonds and their four neighbours are within 0.25 Å of the closest shared plane.
            - Internal steric clash: The interatomic distance between pairs of non-covalently bound atoms is above 60% of the lower bound distance apart determined by distance geometry.


Usage Option: Only Place Scaffolds (or Specifically Don't Place)
===============================================================

You can run Syndirella to only place scaffolds. It will not perform the full elaboration procedure.

.. code-block:: bash

    syndirella run --input [path_to_automatic.csv] --output [path_to_output_dir] --templates [path_to_templates_dir]
    --hits_path [path_to_fragments.sdf] --metadata [path_to_metadata.csv] --scaffold_place

You can also specify to not place the scaffold (most likely you confirmed placement using another method).

.. code-block:: bash

    syndirella run --input [path_to_automatic.csv] --output [path_to_output_dir] --templates [path_to_templates_dir]
    --hits_path [path_to_fragments.sdf] --metadata [path_to_metadata.csv] --no_scaffold_place


Usage Option: Only Get Retrosynthesis Routes of Scaffolds
========================================================

You can run Syndirella to find the Top 5 retrosynthesis routes of the scaffolds. It will identify the routes that contains
all reactions you have encoded in the RXN_SMIRKS_CONSTANTS.json file (a CAR route) and routes that don't contain those
reactions (non-CAR route).

.. code-block:: bash

    syndirella run --input [path_to_automatic.csv] --output [path_to_output_dir] --just_retro

**Output file: [input_csv_name].pkl.gz**

.. note::

    You can read this file using pandas and reading it in as a pickle.

Structure of the important columns are:

    ``routeX``:
        List of dictionaries of each step (X is an int) in the route with reaction names, reactants, and product.
    ``routeX_names``:
        List of reaction names in the route.
    ``routeX_CAR``:
        Boolean value if all reactions in route are in RXN_SMIRKS_CONSTANTS.json.
    ``routeX_non_CAR``:
        List of reaction names that are not in RXN_SMIRKS_CONSTANTS.json. Or None if all reactions are in RXN_SMIRKS_CONSTANTS.json.

If there are `NaN` values for all above columns, it means that there are no routes found for the scaffold.

Usage Option: Only Elaborate One Reactant per Series
========================================================

.. attention::

    This functionality is only provided for single step reactions.

You can have Syndirella output elaboration series for one reactant at a time. For example, if the route is a single step
amidation, there will be two elaboration series output: (1) only elaborating reactant 1 and (2) only elaborating reactant 2.

.. note::

    Each series per reactant will be handled as seperate, so they will have their own unique route uuids. If an
    alternative route is found for the original route, the alternative route will produce two seperate series as well
    for each reactant elaboration.

.. code-block:: bash

    syndirella run --input [path_to_input.csv] --output [path_to_output_dir] --templates [path_to_templates_dir]
    --hits_path [path_to_fragments.sdf] --metadata [path_to_metadata.csv] --elab_single_reactant



