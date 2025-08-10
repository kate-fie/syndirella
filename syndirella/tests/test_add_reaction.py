#!/usr/bin/env python3
"""
Tests for the add-reaction functionality
"""
import json
import os
import tempfile
import pytest
from unittest.mock import patch

from syndirella.route.SmirksLibraryManager import SmirksLibraryManager


@pytest.fixture
def temp_constants_dir():
    """Create a temporary directory with test constants files"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create test SMIRKS data
        test_smirks_data = {
            "Amidation": {
                "smirks": "[#6:1](=[#8:2])-[#8;H1].[$([N+0&H1,N+0&H2]);!$(NC=*);!$(NS);!$(N=*);!$(N-O);!$(N-o):3]>>[#6:1](=[#8:2])-[#7X3:3]",
                "source": "manual",
                "type": "parent"
            },
            "Suzuki_coupling": {
                "smirks": "[#6:1]-[Cl,Br,I].[#6:2]-[B](-[O])(-[O])>>[#6:1]-[#6:2]",
                "source": "manual",
                "type": "parent"
            }
        }
        
        # Create test USPTO data
        test_uspto_data = {
            "730": {
                "uspto_template": "[O;D1;H0:2]=[CH;D2;+0:1]-[NH;D2;+0:3]-[c:4]>>O-[CH;D2;+0:1]=[O;D1;H0:2].[NH2;D1;+0:3]-[c:4]",
                "similarity_method": "jaccard_maccs_concat_rxn_fp",
                "mappings": [
                {
                    "reaction_name": "Amidation",
                    "similarity": 0.8125,
                    "type": "parent",
                    "source": "manual"
                }]
            }
        }
        
        # Write test files
        smirks_path = os.path.join(temp_dir, 'RXN_SMIRKS_CONSTANTS.json')
        uspto_path = os.path.join(temp_dir, 'uspto_template_lookup.json')
        
        with open(smirks_path, 'w') as f:
            json.dump(test_smirks_data, f, indent=2)
        
        with open(uspto_path, 'w') as f:
            json.dump(test_uspto_data, f, indent=2)
        
        yield temp_dir


@pytest.fixture
def smirks_manager(temp_constants_dir):
    """Create a SmirksLibraryManager instance with test data"""
    smirks_path = os.path.join(temp_constants_dir, 'RXN_SMIRKS_CONSTANTS.json')
    uspto_path = os.path.join(temp_constants_dir, 'uspto_template_lookup.json')
    return SmirksLibraryManager(smirks_path, uspto_path)


def test_add_reaction_basic(smirks_manager):
    """Test basic reaction addition without parent finding"""
    name = "Test_New_Reaction"
    smirks = "[#6:1]-[OH].[#6:2]-[Cl]>>[#6:1]-[O]-[#6:2]"
    
    result = smirks_manager.add_reaction_to_library(
        name=name,
        smirks=smirks,
        find_parent=False,
        fp_type='maccs_rxn_fp',
        threshold=0.0,
        similarity_metric='cosine'
    )
    
    print(result)

    assert result['reaction_name'] == name
    assert result['parent_reaction'] is None
    assert isinstance(result['similar_templates'], list)
    assert len(result['similar_templates']) > 0
    
    # Verify reaction was added to smirks data
    assert name in smirks_manager.smirks_data
    assert smirks_manager.smirks_data[name]['smirks'] == smirks
    assert smirks_manager.smirks_data[name]['type'] == 'parent'



def test_add_reaction_with_parent_finding(smirks_manager):
    """Test reaction addition with parent finding enabled"""
    name = "Test_Child_Reaction"
    smirks = "[#6:1](=[#8:2])-[#8;H1].[$([N+0&H1,N+0&H2]);!$(NC=*);!$(NS);!$(N=*);!$(N-O);!$(N-o):3]>>[#6:1](=[#8:2])-[#7X3:3]"
    
    result = smirks_manager.add_reaction_to_library(
        name=name,
        smirks=smirks,
        find_parent=True,
        fp_type='maccs_rxn_fp',
        threshold=0.1,  # Low threshold to ensure match
        similarity_metric='cosine'
    )
    
    assert result['reaction_name'] == name
    assert result['parent_reaction'] == "Amidation"  # Should match the existing Amidation reaction
    assert isinstance(result['similar_templates'], list)
    
    # Verify reaction was added as child
    assert name in smirks_manager.smirks_data
    assert smirks_manager.smirks_data[name]['type'] == 'child'
    assert smirks_manager.smirks_data[name]['parent'] == "Amidation"


def test_add_reaction_morgan_fingerprint(smirks_manager):
    """Test reaction addition with Morgan fingerprints"""
    name = "Test_Morgan_Reaction"
    smirks = "[#6:1]-[Cl,Br,I].[#6:2]-[B](-[O])(-[O])>>[#6:1]-[#6:2]"
    
    result = smirks_manager.add_reaction_to_library(
        name=name,
        smirks=smirks,
        find_parent=True,
        fp_type='morgan_rxn_fp',
        threshold=0.1,
        similarity_metric='jaccard'
    )
    
    assert result['reaction_name'] == name
    assert isinstance(result['similar_templates'], list)


def test_add_reaction_no_parent_found(smirks_manager):
    """Test reaction addition when no similar parent is found"""
    name = "Test_Unique_Reaction"
    smirks = "[#6:1]-[#6:2]-[#6:3]>>[#6:1]=[#6:2]-[#6:3]"  # Very different SMIRKS
    
    result = smirks_manager.add_reaction_to_library(
        name=name,
        smirks=smirks,
        find_parent=True,
        fp_type='maccs_rxn_fp',
        threshold=0.9,  # High threshold to ensure no match
        similarity_metric='cosine'
    )
    
    assert result['reaction_name'] == name
    assert result['parent_reaction'] is None
    assert smirks_manager.smirks_data[name]['type'] == 'parent'


def test_add_reaction_invalid_smirks(smirks_manager):
    """Test reaction addition with invalid SMIRKS"""
    name = "Test_Invalid_Reaction"
    invalid_smirks = "invalid_smirks_string"
    
    with pytest.raises(Exception):
        smirks_manager.add_reaction_to_library(
            name=name,
            smirks=invalid_smirks,
            find_parent=False,
            fp_type='maccs_rxn_fp',
            threshold=0.2,
            similarity_metric='cosine'
        )


def test_add_reaction_invalid_fp_type(smirks_manager):
    """Test reaction addition with invalid fingerprint type"""
    name = "Test_Invalid_FP_Reaction"
    smirks = "[#6:1]-[OH].[#6:2]-[Cl]>>[#6:1]-[O]-[#6:2]"
    
    with pytest.raises(KeyError):
        smirks_manager.add_reaction_to_library(
            name=name,
            smirks=smirks,
            find_parent=False,
            fp_type='invalid_fp_type',
            threshold=0.2,
            similarity_metric='cosine'
        )


def test_add_reaction_invalid_similarity_metric(smirks_manager):
    """Test reaction addition with invalid similarity metric"""
    name = "Test_Invalid_Metric_Reaction"
    smirks = "[#6:1]-[OH].[#6:2]-[Cl]>>[#6:1]-[O]-[#6:2]"
    
    with pytest.raises(ValueError):
        smirks_manager.add_reaction_to_library(
            name=name,
            smirks=smirks,
            find_parent=False,
            fp_type='maccs_rxn_fp',
            threshold=0.2,
            similarity_metric='invalid_metric'
        )


def test_add_reaction_duplicate_name(smirks_manager):
    """Test adding a reaction with a name that already exists"""
    name = "Amidation"  # This already exists in test data
    smirks = "[#6:1]-[OH].[#6:2]-[Cl]>>[#6:1]-[O]-[#6:2]"
    
    # Should overwrite existing reaction
    result = smirks_manager.add_reaction_to_library(
        name=name,
        smirks=smirks,
        find_parent=False,
        fp_type='maccs_rxn_fp',
        threshold=0.2,
        similarity_metric='cosine'
    )
    
    assert result['reaction_name'] == name
    assert smirks_manager.smirks_data[name]['smirks'] == smirks


def test_add_reaction_persistence(temp_constants_dir):
    """Test that added reactions are persisted to files"""
    smirks_path = os.path.join(temp_constants_dir, 'RXN_SMIRKS_CONSTANTS.json')
    uspto_path = os.path.join(temp_constants_dir, 'uspto_template_lookup.json')
    
    manager = SmirksLibraryManager(smirks_path, uspto_path)
    
    name = "Test_Persistent_Reaction"
    smirks = "[#6:1]-[OH].[#6:2]-[Cl]>>[#6:1]-[O]-[#6:2]"
    
    manager.add_reaction_to_library(
        name=name,
        smirks=smirks,
        find_parent=False,
        fp_type='maccs_rxn_fp',
        threshold=0.2,
        similarity_metric='cosine'
    )
    
    # Verify file was updated
    with open(smirks_path, 'r') as f:
        updated_data = json.load(f)
    
    assert name in updated_data
    assert updated_data[name]['smirks'] == smirks 