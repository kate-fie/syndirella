==================
Database Search
==================

Syndirella offers three options for database search functionality: Arthor, Hippo, and Postera (Manifold). All three tools support superstructure searches to find commercially available compounds that can be used as reactants in synthesis.

Arthor
------

Arthor is a free, publicly accessible database search tool that provides access to multiple chemical databases including ZINC, Enamine, Mcule, and stock databases.

**Usage:**

When running Syndirella, specify Arthor as the database search tool:

.. code-block:: bash

    syndirella run --input input.csv --output output_dir --db_search_tool arthor

**Default Vendors:**

The default vendor list includes:
- ``enamine_real`` - Enamine REAL Database
- ``stock`` - In-Stock database
- ``mcule`` - Mcule databases (includes multiple Mcule variants)
- ``zinc`` - ZINC databases (includes multiple ZINC variants)

**Available Vendor Options:**

Arthor supports a wide range of vendor databases: check the `Arthor documentation <https://arthor.docking.org/>`_ for more details.

**Checking Available Databases:**

You can programmatically check which databases are available:

.. code-block:: python

    from syndirella.database.Arthor import Arthor
    
    arthor = Arthor()
    available_dbs = arthor.get_available_databases()
    print(available_dbs)

Postera (Manifold)
------------------

Postera's Manifold API provides access to a comprehensive database of commercially available compounds. It supports superstructure searches across multiple vendor catalogs.

**Setup:**

1. Get your Manifold API key from your `Postera account <https://app.postera.ai/>`_.

2. Set the required environment variables:

.. code-block:: bash

   export MANIFOLD_API_KEY=[YOUR_API_KEY]
   export MANIFOLD_API_URL=[API_URL]

   # Example:
   # export MANIFOLD_API_KEY=your_api_key_here
   # export MANIFOLD_API_URL=https://api.postera.ai

.. attention::

    A Manifold Pro account is required for API access.

**Usage:**

When running Syndirella, specify Postera as the database search tool:

.. code-block:: bash

    syndirella run --input input.csv --output output_dir --db_search_tool manifold

**Default Vendors:**

The default vendor list includes:
- ``enamine_bb`` - Enamine Building Blocks
- ``mcule`` - Mcule database
- ``mcule_ultimate`` - Mcule Ultimate database
- ``enamine_real`` - Enamine REAL database
- ``enamine_made`` - Enamine Made database

Hippo
-----

Hippo provides local database search functionality using a HIPPO database file. This option is useful when you have your own database or want to search locally without API dependencies.

**Setup:**

Have a HIPPO database file available and initialize the Hippo search with the database path.

**Usage:**

Specify Hippo as the database search tool:

.. code-block:: bash

    syndirella run --input input.csv --output output_dir --db_search_tool hippo --reference_db /path/to/database

**Default Behavior:**

If no database search tool is specified, Syndirella will use the default tool (Arthor).
