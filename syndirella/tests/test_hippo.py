#!/usr/bin/env python3
"""
Tests for HIPPO database search (graceful import, error when not installed, enum).
No hippo-db or real HIPPO database required; tests pass with or without hippo-db installed.
When hippo-db is installed, tests use a temp directory path so HIPPO can create/open a DB.
"""
import os
import tempfile
import unittest

from syndirella.database.Hippo import Hippo
from syndirella.constants import DatabaseSearchTool
from syndirella.utils.error import HippoNotInstalledError


def _hippo_test_db_path():
    """Return a path in an existing temp dir for HIPPO DB (so hippo-db can create it when installed)."""
    tmp = tempfile.mkdtemp(prefix="syndirella_hippo_test_")
    return os.path.join(tmp, "reference_db.sqlite")


class TestHippo(unittest.TestCase):
    """Test HIPPO module import, instantiation, and error behavior when hippo-db is not installed."""

    def test_graceful_import(self):
        """Importing Hippo does not raise even when hippo-db is not installed."""
        from syndirella.database.Hippo import Hippo as HippoCls
        self.assertIs(HippoCls, Hippo)

    def test_hippo_instantiation(self):
        """Hippo can be instantiated; db_path is set."""
        db_path = _hippo_test_db_path()
        h = Hippo(db_path)
        self.assertEqual(h.db_path, db_path)

    def test_perform_superstructure_search_raises_when_animal_none(self):
        """When hippo-db is not installed, perform_superstructure_search raises with install message."""
        db_path = _hippo_test_db_path()
        h = Hippo(db_path)
        if h.animal is not None:
            self.skipTest("hippo-db is installed; cannot test error path")
        with self.assertRaises(HippoNotInstalledError) as ctx:
            h.perform_superstructure_search("CCO")
        self.assertIn("syndirella[hippo]", str(ctx.exception))
        self.assertIn("pip install", str(ctx.exception))

    def test_perform_superstructure_search_non_string_raises_type_error(self):
        """perform_superstructure_search with non-string raises TypeError."""
        db_path = _hippo_test_db_path()
        h = Hippo(db_path)
        with self.assertRaises(TypeError) as ctx:
            h.perform_superstructure_search(123)
        self.assertIn("Smiles must be a string", str(ctx.exception))

    def test_database_search_tool_hippo_enum(self):
        """DatabaseSearchTool.HIPPO exists and from_string('hippo') works."""
        self.assertEqual(DatabaseSearchTool.HIPPO.value, "hippo")
        self.assertIs(DatabaseSearchTool.from_string("hippo"), DatabaseSearchTool.HIPPO)


if __name__ == "__main__":
    unittest.main()
