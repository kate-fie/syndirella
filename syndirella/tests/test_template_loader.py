"""
Tests for the data loader utility.
"""
import pytest
import tempfile
import os
from pathlib import Path
from syndirella.utils.template_loader import TemplateLoader, load_uspto_lookup


def test_template_loader_initialization():
    """Test TemplateLoader initialization."""
    loader = TemplateLoader()
    assert loader.cache_dir.exists()
    assert loader.cache_dir.is_dir()


def test_template_loader_with_custom_cache():
    """Test TemplateLoader with custom cache directory."""
    with tempfile.TemporaryDirectory() as temp_dir:
        loader = TemplateLoader(cache_dir=temp_dir)
        assert loader.cache_dir == Path(temp_dir)


def test_load_uspto_lookup():
    """Test loading USPTO lookup data."""
    try:
        data = load_uspto_lookup()
        assert isinstance(data, dict)
        # Check that it has the expected structure
        if data:  # If data is not empty
            first_key = list(data.keys())[0]
            assert isinstance(data[first_key], dict)
            assert 'uspto_template' in data[first_key]
    except FileNotFoundError:
        # This is expected if the gzipped file doesn't exist in test environment
        pytest.skip("USPTO lookup file not available for testing")


def test_cache_functionality():
    """Test that caching works correctly."""
    with tempfile.TemporaryDirectory() as temp_dir:
        loader = TemplateLoader(cache_dir=temp_dir)
        
        # Create a test gzipped JSON file
        import gzip
        import json
        
        test_data = {"test": "data", "numbers": [1, 2, 3]}
        test_gz_path = Path(temp_dir) / "test.json.gz"
        
        with gzip.open(test_gz_path, 'wt', encoding='utf-8') as f:
            json.dump(test_data, f)
        
        # Load the data
        loaded_data = loader.load_json_gz(str(test_gz_path), "test_cache")
        
        # Verify data is correct
        assert loaded_data == test_data
        
        # Verify cache file was created
        cache_file = loader.cache_dir / "test_cache.json"
        assert cache_file.exists()
        
        # Load again - should use cache
        loaded_data2 = loader.load_json_gz(str(test_gz_path), "test_cache")
        assert loaded_data2 == test_data


def test_cache_size():
    """Test cache size calculation."""
    with tempfile.TemporaryDirectory() as temp_dir:
        loader = TemplateLoader(cache_dir=temp_dir)
        
        # Initially should be 0
        assert loader.get_cache_size() == 0
        
        # Create a test cache file
        cache_file = loader.cache_dir / "test.json"
        cache_file.write_text('{"test": "data"}')
        
        # Should have some size
        assert loader.get_cache_size() > 0


def test_clear_cache():
    """Test cache clearing functionality."""
    with tempfile.TemporaryDirectory() as temp_dir:
        loader = TemplateLoader(cache_dir=temp_dir)
        
        # Create some test cache files
        cache_file1 = loader.cache_dir / "test1.json"
        cache_file2 = loader.cache_dir / "test2.json"
        
        cache_file1.write_text('{"test1": "data"}')
        cache_file2.write_text('{"test2": "data"}')
        
        assert len(list(loader.cache_dir.glob("*.json"))) == 2
        
        # Clear cache
        loader.clear_cache()
        
        assert len(list(loader.cache_dir.glob("*.json"))) == 0
