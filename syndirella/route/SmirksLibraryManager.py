"""
SmirksLibraryManager Module - Manages structured SMIRKS library with metadata and USPTO template lookups.
"""

import json
import logging
import os
import numpy as np
import pandas as pd
from collections import defaultdict
from typing import Dict, List, Optional, Any

from syndirella.utils.classifier import get_fp, classify_reaction, calc_cosine_similarity, calc_jaccard_similarity
from syndirella.utils.template_loader import load_uspto_lookup


class SmirksLibraryManager:
    """Manages the structured SMIRKS library with automatic loading of metadata."""

    def __init__(self,
                 smirks_library_path: str,
                 uspto_lookup_path: str = None,
                 logger: logging.Logger = None):
        self.logger = logger or logging.getLogger(__name__)
        self.smirks_library_path = smirks_library_path

        # Auto-determine USPTO lookup path if not provided
        if uspto_lookup_path is None:
            base_dir = os.path.dirname(smirks_library_path)
            self.uspto_lookup_path = os.path.join(base_dir, "uspto_template_lookup.json")
        else:
            self.uspto_lookup_path = uspto_lookup_path

        # Load everything
        self._load_libraries()

    def _load_libraries(self):
        """Load the SMIRKS library and USPTO lookup."""
        # Load main SMIRKS library (now contains metadata structure)
        with open(self.smirks_library_path, 'r') as f:
            self.smirks_data = json.load(f)

        # Extract just the SMIRKS strings for compatibility
        self.all_reactions = {}
        self.reaction_metadata = {}

        for name, data in self.smirks_data.items():
            if isinstance(data, dict) and 'smirks' in data:
                # New structured format
                self.all_reactions[name] = data['smirks']
                self.reaction_metadata[name] = {
                    'source': data.get('source', 'unknown'),
                    'type': data.get('type', 'unknown'),
                    'parent': data.get('parent', None)
                }
            else:
                # Handle old format where data is just the SMIRKS string
                self.all_reactions[name] = data
                self.reaction_metadata[name] = {
                    'source': 'manual',
                    'type': 'parent',
                    'parent': None
                }

        # Load USPTO template lookup using template loader
        try:
            self.uspto_lookup = load_uspto_lookup()
            self.logger.info("Loaded USPTO lookup from compressed file")
        except Exception as e:
            self.uspto_lookup = {}
            self.logger.warning(f"USPTO lookup file not found or could not be loaded: {e}")

        self.logger.info(f"Loaded {len(self.all_reactions)} total reactions")
        self.logger.info(f"Loaded USPTO mappings for {len(self.uspto_lookup)} template codes")

    def reload_libraries(self):
        """Reload libraries from disk (useful for development/testing)."""
        self.logger.info("Reloading SMIRKS libraries from disk...")
        self._load_libraries()

    def get_all_reactions(self) -> Dict[str, str]:
        """Get all reactions as name -> SMIRKS mapping."""
        return self.all_reactions.copy()

    def get_reaction_smirks(self, reaction_name: str) -> Optional[str]:
        """Get SMIRKS for a specific reaction."""
        return self.all_reactions.get(reaction_name)

    def get_reaction_metadata(self, reaction_name: str) -> Optional[Dict]:
        """Get metadata for a specific reaction."""
        return self.reaction_metadata.get(reaction_name)

    def get_uspto_mappings_by_code(self, template_code: str) -> Dict:
        """Get SMIRKS mappings for a USPTO template code."""
        return self.uspto_lookup.get(str(template_code), {})

    def get_uspto_mappings_by_code_and_template(self, template_code: str, template: str) -> Dict:
        """Get SMIRKS mappings for a USPTO template code and template."""
        template_data = self.uspto_lookup.get(str(template_code), {})
        if not template_data or 'uspto_template' not in template_data:
            return {}
        if template == template_data['uspto_template']:
            return template_data
        else:
            return {}

    def get_parent_reactions(self) -> Dict[str, str]:
        """Get only parent reactions."""
        return {name: smirks for name, smirks in self.all_reactions.items()
                if self.reaction_metadata[name]['type'] == 'parent'}

    def get_child_reactions(self) -> Dict[str, str]:
        """Get only child reactions."""
        return {name: smirks for name, smirks in self.all_reactions.items()
                if self.reaction_metadata[name]['type'] == 'child'}

    def get_manual_reactions(self) -> Dict[str, str]:
        """Get only manually curated reactions."""
        return {name: smirks for name, smirks in self.all_reactions.items()
                if self.reaction_metadata[name]['source'] == 'manual'}

    def get_rxn_insight_reactions(self) -> Dict[str, str]:
        """Get only reactions from RXN insight paper."""
        return {name: smirks for name, smirks in self.all_reactions.items()
                if self.reaction_metadata[name]['source'] == 'rxn-insight'}

    def is_child(self, reaction_name: str) -> bool:
        """Check if a reaction is a child reaction."""
        return reaction_name in self.get_child_reactions()

    def get_children_of_parent(self, parent_name: str) -> List[str]:
        """Get all child reactions for a given parent."""
        return [name for name, metadata in self.reaction_metadata.items()
                if metadata['parent'] == parent_name]

    def get_parent_of_child(self, child_name: str) -> Optional[str]:
        """Get the parent reaction for a given child."""
        metadata = self.get_reaction_metadata(child_name)
        return metadata['parent'] if metadata else None

    def get_parent_child_relationships(self) -> Dict[str, List[str]]:
        """Get parent -> children mapping."""
        relationships = defaultdict(list)

        for reaction_name, metadata in self.reaction_metadata.items():
            if metadata['parent']:
                relationships[metadata['parent']].append(reaction_name)

        return dict(relationships)

    def get_reaction_family(self, reaction_name: str) -> Dict[str, List[str]]:
        """Get the complete family (parent + all siblings) for a given reaction."""
        metadata = self.get_reaction_metadata(reaction_name)
        if not metadata:
            return {'parent': None, 'siblings': []}

        if metadata['type'] == 'parent':
            # This is a parent, get all its children
            children = self.get_children_of_parent(reaction_name)
            return {'parent': reaction_name, 'siblings': children}
        else:
            # This is a child, get parent and all siblings
            parent = metadata['parent']
            if parent:
                siblings = [name for name in self.get_children_of_parent(parent)
                            if name != reaction_name]
                return {'parent': parent, 'siblings': siblings}
            else:
                return {'parent': None, 'siblings': []}

    def find_reactions_by_pattern(self, pattern: str, search_in: str = 'name') -> List[str]:
        """Find reactions containing a pattern."""
        pattern_lower = pattern.lower()
        matches = []

        for name, smirks in self.all_reactions.items():
            if search_in in ['name', 'both'] and pattern_lower in name.lower():
                matches.append(name)
            elif search_in in ['smirks', 'both'] and pattern_lower in smirks.lower():
                matches.append(name)

        return matches

    def get_uspto_template_info(self, template_code: str) -> Dict:
        """Get detailed information about a USPTO template."""
        template_data = self.get_uspto_mappings_by_code(template_code)
        if not template_data:
            return {}

        info = {
            'template_code': template_code,
            'uspto_template': template_data.get('uspto_template', ''),
            'similarity_method': template_data.get('similarity_method', ''),
            'num_mappings': len(template_data.get('mappings', [])),
            'mapped_reactions': [m['reaction_name'] for m in template_data.get('mappings', [])],
            'best_similarity': max([m['similarity'] for m in template_data.get('mappings', [])
                                    if m['similarity'] is not None], default=0)
        }

        return info

    def get_library_stats(self) -> Dict:
        """Get comprehensive statistics about the library."""
        base_reactions = sum(1 for m in self.reaction_metadata.values() if m['source'] == 'manual')
        extended_reactions = sum(1 for m in self.reaction_metadata.values() if m['source'] == 'rxn-insight')
        parent_count = sum(1 for m in self.reaction_metadata.values() if m['type'] == 'parent')
        child_count = sum(1 for m in self.reaction_metadata.values() if m['type'] == 'child')

        # USPTO template stats
        templates_with_mappings = len([t for t in self.uspto_lookup.values() if t.get('mappings')])
        total_mappings = sum(len(t.get('mappings', [])) for t in self.uspto_lookup.values())

        return {
            "total_reactions": len(self.all_reactions),
            "manual_reactions": base_reactions,
            "rxn-insight_reactions": extended_reactions,
            "parent_reactions": parent_count,
            "child_reactions": child_count,
            "uspto_templates_total": len(self.uspto_lookup),
            "uspto_templates_with_mappings": templates_with_mappings,
            "total_uspto_mappings": total_mappings,
            "parent_child_relationships": len(self.get_parent_child_relationships())
        }

    def validate_library_integrity(self) -> Dict:
        """Validate the integrity of the library structure."""
        issues = []
        warnings = []

        # Check that all child reactions have valid parents
        for name, metadata in self.reaction_metadata.items():
            if metadata['type'] == 'child':
                parent = metadata['parent']
                if not parent:
                    issues.append(f"Child reaction '{name}' has no parent specified")
                elif parent not in self.all_reactions:
                    issues.append(f"Child reaction '{name}' has invalid parent '{parent}'")
                elif self.reaction_metadata[parent]['type'] != 'parent':
                    issues.append(f"Child reaction '{name}' parent '{parent}' is not marked as parent type")

        # Check for orphaned parents (parents with no children)
        relationships = self.get_parent_child_relationships()
        for name, metadata in self.reaction_metadata.items():
            if metadata['type'] == 'parent' and name not in relationships:
                warnings.append(f"Parent reaction '{name}' has no child reactions")

        # Check USPTO template consistency
        for template_code, template_data in self.uspto_lookup.items():
            mappings = template_data.get('mappings', [])
            for mapping in mappings:
                reaction_name = mapping['reaction_name']
                if reaction_name not in self.all_reactions:
                    issues.append(f"USPTO template {template_code} maps to non-existent reaction '{reaction_name}'")

        return {
            'valid': len(issues) == 0,
            'issues': issues,
            'warnings': warnings,
            'summary': f"{len(issues)} issues, {len(warnings)} warnings found"
        }

    def export_summary(self, output_path: str = None) -> Dict:
        """Export a comprehensive summary of the library."""
        stats = self.get_library_stats()
        validation = self.validate_library_integrity()
        relationships = self.get_parent_child_relationships()

        summary = {
            'library_stats': stats,
            'validation': validation,
            'parent_child_relationships': relationships,
            'sample_reactions': {
                'parents': list(self.get_parent_reactions().keys())[:5],
                'children': list(self.get_child_reactions().keys())[:5]
            },
            'sample_uspto_templates': list(self.uspto_lookup.keys())[:5]
        }

        if output_path:
            with open(output_path, 'w') as f:
                json.dump(summary, f, indent=2)

        return summary

    def add_reaction_to_library(self, name: str, smirks: str, find_parent: bool = False, 
                               fp_type: str = 'maccs_rxn_fp', threshold: float = 0.2,
                               similarity_metric: str = 'cosine') -> Dict[str, Any]:
        """Add a new reaction to the library with similarity analysis."""
        # Save original data for potential rollback
        original_smirks_data = self.smirks_data.copy()
        original_uspto_data = self.uspto_lookup.copy()
        
        try:
            # Calculate fingerprint for new reaction
            fp_type_clean = fp_type.replace('_rxn_fp', '').upper()
            new_reaction_fp = get_fp(smirks, fp=fp_type_clean, concatenate=True)
            
            # Convert existing reactions to DataFrame for similarity calculation
            reaction_data = []
            for rxn_name, rxn_info in self.smirks_data.items():
                if isinstance(rxn_info, dict) and 'smirks' in rxn_info:
                    rxn_smirks = rxn_info['smirks']
                else:
                    rxn_smirks = rxn_info
                    
                rxn_fp = get_fp(rxn_smirks, fp=fp_type_clean, concatenate=True)
                reaction_data.append({
                    'name': rxn_name,
                    'smirks': rxn_smirks,
                    'reaction_fp': rxn_fp
                })
            
            smirks_df = pd.DataFrame(reaction_data)
            
            # Find similar reactions if find_parent is True
            parent_reaction = None
            if find_parent:
                similar_reactions = classify_reaction(
                    new_reaction_fp, 
                    smirks_df, 
                    similarity_metric, 
                    top_n=5, 
                    threshold=threshold
                )
                
                if similar_reactions:
                    best_match = similar_reactions[0]
                    self.logger.info(f"Found similar reaction: {best_match['name']} (similarity: {best_match['similarity']:.3f})")
                    parent_reaction = best_match['name']
                else:
                    self.logger.info("No similar parent reaction found above threshold.")
            
            # Add new reaction to smirks_data (in memory only)
            new_reaction_info = {
                "smirks": smirks,
                "source": "manual",
                "type": "child" if find_parent and parent_reaction else "parent",
                "parent": parent_reaction if find_parent and parent_reaction else None
            }
            
            self.smirks_data[name] = new_reaction_info
            
            # Find similar USPTO templates
            similar_templates = []
            self.logger.info(f"Finding similar USPTO templates for {name} with {similarity_metric} similarity metric")
            for template_code, template_info in self.uspto_lookup.items():
                if 'uspto_template' in template_info:
                    template_rxn_fp = get_fp(template_info['uspto_template'], fp=fp_type_clean, concatenate=True)
                    if similarity_metric == 'cosine':
                        sim = calc_cosine_similarity(new_reaction_fp, template_rxn_fp)
                    elif similarity_metric == 'jaccard':
                        sim = calc_jaccard_similarity(new_reaction_fp, template_rxn_fp)
                    else:
                        raise ValueError(f"Unknown similarity metric: {similarity_metric}")

                    if sim >= threshold:
                        similar_templates.append({
                            'template_code': template_code,
                            'similarity': sim,
                            'mapped_reaction': template_info['uspto_template']
                        })
            
            # Sort by similarity
            similar_templates.sort(key=lambda x: x['similarity'], reverse=True)
            
            # Add new reaction to similar USPTO templates (in memory only)
            for template_info in similar_templates:
                template_code = template_info['template_code']
                similarity = template_info['similarity']
                
                if template_code not in self.uspto_lookup:
                    self.uspto_lookup[template_code] = {
                        'uspto_template': '',
                        'similarity_method': similarity_metric,
                        'mappings': []
                    }
                
                new_mapping = {
                    'reaction_name': name,
                    'similarity': similarity,
                    'fp_type': fp_type,
                    'threshold': threshold
                }
                
                self.uspto_lookup[template_code]['mappings'].append(new_mapping)
                self.logger.info(f"Added reaction '{name}' to USPTO template {template_code} (similarity: {similarity:.3f})")
            
            # Validate the updated data structure
            self._validate_updated_data(name, new_reaction_info, similar_templates)
            
            # Only after validation passes, save to files
            self._save_updated_files()
            
            self.logger.info(f"Successfully added reaction '{name}' with {len(similar_templates)} USPTO template matches")
            
            return {
                'reaction_name': name,
                'parent_reaction': parent_reaction,
                'similar_templates': similar_templates
            }
            
        except Exception as e:
            # Rollback changes on any error
            self.logger.error(f"Error adding reaction '{name}': {e}")
            self.logger.info("Rolling back changes...")
            self.smirks_data = original_smirks_data
            self.uspto_lookup = original_uspto_data
            raise
    
    def _validate_updated_data(self, name: str, new_reaction_info: Dict, similar_templates: List[Dict]):
        """Validate the updated data structure before saving."""
        # Validate SMIRKS data
        if name not in self.smirks_data:
            raise ValueError(f"Reaction '{name}' was not properly added to smirks_data")
        
        if self.smirks_data[name] != new_reaction_info:
            raise ValueError(f"Reaction '{name}' data does not match expected structure")
        
        # Validate parent-child relationships
        if new_reaction_info['type'] == 'child' and new_reaction_info['parent']:
            parent_name = new_reaction_info['parent']
            if parent_name not in self.smirks_data:
                raise ValueError(f"Parent reaction '{parent_name}' not found in smirks_data")
            if self.smirks_data[parent_name]['type'] != 'parent':
                raise ValueError(f"Parent reaction '{parent_name}' is not marked as parent type")
        
        # Validate USPTO data structure
        for template_info in similar_templates:
            template_code = template_info['template_code']
            if template_code not in self.uspto_lookup:
                raise ValueError(f"USPTO template '{template_code}' not found in lookup data")
            
            # Check that the new mapping was added
            mappings = self.uspto_lookup[template_code]['mappings']
            found_mapping = False
            for mapping in mappings:
                if mapping['reaction_name'] == name:
                    found_mapping = True
                    break
            
            if not found_mapping:
                raise ValueError(f"Reaction '{name}' mapping not found in USPTO template '{template_code}'")
        
        self.logger.info("Data validation passed")
    
    def _save_updated_files(self):
        """Save the updated data to files."""
        # Save updated SMIRKS data
        with open(self.smirks_library_path, 'w') as f:
            json.dump(self.smirks_data, f, indent=2)
        
        self.logger.info(f"Updated SMIRKS data saved to {self.smirks_library_path}")
        
        # Save updated USPTO data to compressed file
        import gzip
        compressed_path = self.uspto_lookup_path.replace('.json', '.json.gz')
        with gzip.open(compressed_path, 'wt', encoding='utf-8') as f:
            json.dump(self.uspto_lookup, f, indent=2)
        self.logger.info(f"Updated USPTO data saved to compressed file: {compressed_path}")
        
        # Also update the cache
        from syndirella.utils.template_loader import _template_loader
        cache_path = _template_loader.cache_dir / "uspto_template_lookup.json"
        with open(cache_path, 'w') as f:
            json.dump(self.uspto_lookup, f, indent=2)
        self.logger.info(f"Updated USPTO cache at: {cache_path}")
