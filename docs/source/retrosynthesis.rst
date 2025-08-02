==============
Retrosynthesis
==============

Syndirella offers two options for retrosynthesis functionality, AiZynthFinder and Postera's Manifold.

AiZynthFinder
-------------

Basic usage:

.. code-block:: python

    from aizynth.AiZynthManager import AiZynthManager
    manager = AiZynthManager()
    routes = manager.perform_route_search(smiles, matching_strategy="best_overall")

Strategies:
    - "best_overall": Use most similar SMIRKS from whole library
    - "best_child": Use most specific SMIRKS
    - "parent_only": Use only original Syndirella SMIRKS
    - "all_above_threshold": Use any SMIRKS above similarity threshold

Manifold
--------

Syndirella can perform retrosynthesis search using `Postera's Manifold <https://app.postera.ai/>`_ API to perform retrosynthesis searches of input scaffolds
and superstructure searches of reactants. You can find more information on Manifold API access
`here <https://api.postera.ai/api/v1/docs/>`_.

.. attention::

    A Manifold Pro account is required for API access.

Set the Manifold API key as an environment variable.

.. code-block:: bash

   export MANIFOLD_API_KEY=[API_KEY]