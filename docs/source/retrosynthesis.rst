==============
Retrosynthesis
==============

Syndirella offers two options for retrosynthesis functionality, AiZynthFinder and Postera's Manifold.

AiZynthFinder
-------------

Before using AiZynthFinder for retrosynthesis, you need to set it up. There are two options:

**Automatic Setup (Recommended):**

Run the setup command to automatically download required data and create the configuration file:

.. code-block:: bash

    syndirella setup-aizynth

.. attention::

    ⚠️ This will download large model files (~750MB total) if run the first time (will take ~5 min). These are required to run AiZynthFinder. This will automatically download the required data to ``[syndirella_package_path]/aizynth`` directory and create the configuration file automatically.

**Manual Setup (Alternative):**

If you prefer to set up AiZynthFinder manually:

.. code-block:: bash

    cd [syndirella_package_path]/aizynth
    download_public_data .
    # Update config.yml (if you prefer)
    export AIZYNTH_CONFIG_FILE="path/to/syndirella/aizynth/config.yml"

.. note::

    To make the environment variable permanent, add the export line to your shell profile (~/.bashrc, ~/.zshrc, etc.).

Basic usage:

.. code-block:: python

    from syndirella.aizynth.AiZynthManager import AiZynthManager
    manager = AiZynthManager()
    routes = manager.perform_route_search(smiles, matching_strategy="best_overall")

Strategies:
    - "best_overall": Use most similar SMIRKS from whole library
    - "best_child": Use most specific SMIRKS
    - "parent_only": Use only original Syndirella SMIRKS
    - "all_above_threshold": Use any SMIRKS above similarity threshold

.. attention::

    An output file will be created in the output directory per scaffold listing all information about the routes found by AiZynthFinder. 

Manifold
--------

Syndirella can perform retrosynthesis search using `Postera's Manifold <https://app.postera.ai/>`_ API to perform retrosynthesis searches of input scaffolds
and superstructure searches of reactants. You can find more information on Manifold API access
`here <https://api.postera.ai/api/v1/docs/>`_.

.. attention::

    A Manifold Pro account is required for API access.

**Setup:**

1. Get your Manifold API key from your `Postera account <https://app.postera.ai/>`_.

2. Set the required environment variables:

.. code-block:: bash

   export MANIFOLD_API_KEY=[YOUR_API_KEY]
   export MANIFOLD_API_URL=[API_URL]

   # Example:
   # export MANIFOLD_API_KEY=your_api_key_here
   # export MANIFOLD_API_URL=https://api.postera.ai

.. note::

    To make the environment variables permanent, add the export lines to your shell profile:
    
    - For bash: ``~/.bashrc`` or ``~/.bash_profile``
    - For zsh: ``~/.zshrc``
    
    After adding, reload your shell configuration:
    
    .. code-block:: bash
    
       source ~/.bashrc  # or source ~/.zshrc for zsh