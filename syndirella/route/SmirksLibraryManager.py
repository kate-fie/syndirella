"""
SmirksLibraryManager Module - Manages structured SMIRKS library with metadata and USPTO template lookups.
"""

import json
import logging
import os
from collections import defaultdict
from typing import Dict, List, Optional


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

        # Load USPTO template lookup
        if os.path.exists(self.uspto_lookup_path):
            with open(self.uspto_lookup_path, 'r') as f:
                self.uspto_lookup = json.load(f)
        else:
            self.uspto_lookup = {}
            self.logger.warning(f"USPTO lookup file not found: {self.uspto_lookup_path}")

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
            "base_reactions": base_reactions,
            "extended_reactions": extended_reactions,
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

    def check_new_smirks(self):
        """Checks if a new SMIRKS was added to the constants that is not included in the uspto lookup."""
        # TODO: Add this functionality
        return NotImplementedError
