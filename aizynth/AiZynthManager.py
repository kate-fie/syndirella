"""
AiZynthManager Module - Enhanced AiZynthFinder integration with structured SMIRKS library support.
"""

import json
import logging
import os
from collections import defaultdict
from typing import List, Dict, Optional, Tuple

import numpy as np
from aizynthfinder.aizynthfinder import AiZynthFinder
from aizynthfinder.chem.reaction import FixedRetroReaction
from rdkit import Chem
from rdkit.Chem import rdChemReactions
from sklearn.cluster import DBSCAN

from syndirella.error import AiZynthFinderError, NoSynthesisRoute
from syndirella.route.SmirksLibraryManager import SmirksLibraryManager


class AiZynthManager:
    def __init__(self,
                 config_file: str = os.getenv('AIZYNTH_CONFIG_FILE', 'config.yml'),
                 stock: str = 'zinc',
                 expansion_policy: str = 'uspto',
                 filter_policy: str = 'uspto',
                 min_route_score: float = 0.5,
                 clustering_eps: float = 0.2,
                 clustering_min_samples: int = 2,
                 similarity_threshold: float = 0.2,
                 smirks_library_path: str = None,
                 logger: logging.Logger = None):
        # Setup logging
        self.logger = logger or logging.getLogger(__name__)

        # Initialize AiZynthFinder
        self.finder = self._set_finder(config_file, stock, expansion_policy, filter_policy)

        # Set parameters
        self.min_route_score = min_route_score
        self.clustering_eps = clustering_eps
        self.clustering_min_samples = clustering_min_samples
        self.similarity_threshold = similarity_threshold

        # Initialize SMIRKS library manager
        self._initialize_smirks_manager(smirks_library_path)

    def _initialize_smirks_manager(self, smirks_library_path: str = None):
        """Initialize the SMIRKS library manager."""
        if smirks_library_path is None:
            base_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                                    'syndirella', 'constants')
            smirks_library_path = os.path.join(base_dir, 'RXN_SMIRKS_CONSTANTS.json')

        # Initialize manager
        self.smirks_manager = SmirksLibraryManager(
            smirks_library_path=smirks_library_path,
            logger=self.logger
        )

        # Get all available reactions
        self.all_reactions = self.smirks_manager.get_all_reactions()

        # Log library stats
        stats = self.smirks_manager.get_library_stats()
        self.logger.info(f"SMIRKS library loaded: {stats['base_reactions']} base + "
                         f"{stats['extended_reactions']} extended = {stats['total_reactions']} total reactions")
        self.logger.info(
            f"Parent/Child breakdown: {stats['parent_reactions']} parents, {stats['child_reactions']} children")
        self.logger.info(f"USPTO template mappings: {stats['uspto_templates_total']} codes mapped")

    def _set_finder(self, config_file: str, stock: str, expansion_policy: str, filter_policy: str):
        """Set up AiZynthFinder instance."""
        try:
            self.logger.info(f"Loading AiZynthFinder configuration file {config_file}")
            finder = AiZynthFinder(configfile=config_file)
            finder.stock.select(stock)
            finder.expansion_policy.select(expansion_policy)
            finder.filter_policy.select(filter_policy)
            return finder
        except Exception as e:
            self.logger.error(f"Error creating AiZynthFinder instance: {e}")
            raise e

    def _extract_uspto_template_from_reaction(self, reaction) -> Tuple[Optional[int], Optional[str]]:
        """Extract USPTO template code from AiZynthFinder reaction metadata."""
        if not hasattr(reaction, 'metadata') or not reaction.metadata:
            return None, None

        template_code = None
        # Common keys where template code might be stored
        for key in ['template_code', 'template_id']:
            if key in reaction.metadata:
                try:
                    template_code = int(reaction.metadata[key])
                except (ValueError, TypeError):
                    continue

        template = None
        if 'template' in reaction.metadata:
            template = str(reaction.metadata['template'])

        return template_code, template

    def _match_uspto_template_to_list_of_smirks(self,
                                                uspto_template_code: int,
                                                strategy: str = "best_overall",
                                                uspto_template: str = None) -> Optional[List[Tuple[str, str]]]:
        """Match USPTO template to list of SMIRKS using different strategies."""
        if not uspto_template:
            self.logger.debug(f"Finding mappings with just the USPTO template code {uspto_template} ")
            template_data = self.smirks_manager.get_uspto_mappings_by_code(uspto_template_code)
        else:
            template_data = self.smirks_manager.get_uspto_mappings_by_code_and_template(uspto_template_code,
                                                                                        uspto_template)

        if not template_data or 'mappings' not in template_data:
            self.logger.debug(f"No USPTO mapping found for template {uspto_template_code}")
            return None

        mappings = template_data['mappings']

        if strategy == "parent_only":
            # Find just highest similarity parent mappings
            mappings = [m for m in mappings if m["type"] == "parent"]
        elif strategy == "best_child":
            # Find just highest similarity child mappings
            mappings = [m for m in mappings if m["type"] == "child"]
        elif strategy == "best_overall":
            mappings = mappings

        if mappings:
            sorted_mappings = sorted(
                mappings,
                key=lambda x: x["similarity"] if x["similarity"] is not None else 0,
                reverse=True  # descending similarity
            )
            result = []
            for m in sorted_mappings:
                reaction_name = m["reaction_name"]
                smirks = self.smirks_manager.get_reaction_smirks(reaction_name)
                if smirks:
                    result.append((reaction_name, smirks))
            if result:
                return result

        return None

    def _validate_reaction_with_smirks(self, reaction: FixedRetroReaction, smirks: str) -> bool:
        """Validate if a reaction matches a SMIRKS pattern by attempting to apply it."""
        try:
            reaction = reaction.copy()
            pattern = rdChemReactions.ReactionFromSmarts(smirks)

            # Extract reactants and products
            reactant_mols: List[Chem.Mol] = [mol.rd_mol for tuple in reaction.reactants for mol in tuple]
            product_inchi_key: str = reaction.mol.inchi_key

            if any(mol is None for mol in reactant_mols):
                return False

            # Apply reaction
            try:
                result_tuples = []
                result_tuples.extend(pattern.RunReactants(tuple(reactant_mols)))
                # Try swapped order if there are exactly 2 reactants
                if len(reactant_mols) == 2:
                    swapped = [reactant_mols[1], reactant_mols[0]]
                    result_tuples.extend(pattern.RunReactants(tuple(swapped)))
                if not result_tuples:
                    return False
                # Check if any result matches the expected products
                for result_tuple in result_tuples:
                    result_inchi_keys = [Chem.MolToInchiKey(mol) for mol in result_tuple]
                    if product_inchi_key in result_inchi_keys:
                        return True

            except Exception:
                return False

            return False

        except Exception as e:
            self.logger.debug(f"Error validating reaction {reaction_smiles} with SMIRKS {smirks}: {e}")
            return False

    def _match_reaction_to_smirks(self, reaction, strategy: str = "best_overall", validate: bool = True) -> Optional[
        str]:
        """Match a reaction to SMIRKS using USPTO template codes."""
        # Try USPTO template matching first
        uspto_code, uspto_template = self._extract_uspto_template_from_reaction(reaction)

        if uspto_code is not None:
            matches: List = self._match_uspto_template_to_list_of_smirks(uspto_code, strategy,
                                                                         uspto_template=uspto_template)
            if matches:
                for match in matches:
                    reaction_name, smirks = match
                    # Optionally validate the match
                    if validate:
                        if self._validate_reaction_with_smirks(reaction, smirks):
                            self.logger.info(f"Matched and validated USPTO template {uspto_code} to {reaction_name}")
                            return reaction_name
                        else:
                            self.logger.debug(
                                f"USPTO template {uspto_code} matched {reaction_name} but failed validation")
                    else:
                        self.logger.info(f"Matched USPTO template {uspto_code} to {reaction_name} (no validation)")
                        return reaction_name
        return None

    def find_routes(self, target_smiles: str, top_n: int = 10) -> bool:
        """Find synthesis routes for a target molecule."""
        try:
            self.logger.debug(f'Initializing AiZynthFinder instance')
            self.finder.target_smiles = target_smiles
            self.finder.prepare_tree()
            self.finder.tree_search()
            self.finder.build_routes()
            if not self.finder.extract_statistics()['is_solved']:
                self.logger.error(f"Could not find any routes for {target_smiles}")
                return False

            # Store top N routes
            self.routes: list[dict] = list(self.finder.routes)[:top_n]

            self.logger.info(f"Found {len(self.routes)} routes for {target_smiles}")
            return True
        except Exception as e:
            raise AiZynthFinderError(
                message=f'AiZynthFinder failure when finding a route: {e}',
                smiles=target_smiles
            )

    def filter_routes(self) -> List:
        """Filter routes by score and branching criteria."""
        if not hasattr(self, 'routes'):
            self.logger.error("No routes available. Call find_routes first.")
            return []

        self.logger.info('Filtering routes...')
        self.filtered_routes = [
            route for route in self.routes
            if route['route_metadata']['is_solved'] and
               not route['reaction_tree'].is_branched()  # don't include if the route is branched
        ]

        self.logger.info(f"Filtered {len(self.routes)} routes to {len(self.filtered_routes)} routes")
        return self.filtered_routes

    def _compute_route_fingerprints(self, routes: List) -> np.ndarray:
        """Compute OHE fingerprints for routes based on the reactions they contain."""
        # Create a set of all unique reactions across all routes
        all_reactions = set()
        for route in routes:
            for reaction in route['reaction_tree'].reactions():
                # Use reaction SMILES as identifier
                rxn_smiles = reaction.reaction_smiles()
                all_reactions.add(rxn_smiles)

        # Create a mapping from reaction to index
        reaction_to_idx = {rxn: i for i, rxn in enumerate(all_reactions)}

        # Create fingerprints for each route
        fingerprints = np.zeros((len(routes), len(all_reactions)))
        for i, route in enumerate(routes):
            for reaction in route['reaction_tree'].reactions():
                rxn_smiles = reaction.reaction_smiles()
                idx = reaction_to_idx[rxn_smiles]
                fingerprints[i, idx] = 1

        return fingerprints

    def cluster_routes(self) -> Dict[int, List]:
        """Cluster routes based on their reactions."""
        if not hasattr(self, 'filtered_routes'):
            self.logger.error("No filtered routes available. Call filter_routes first.")
            return {}

        if not self.filtered_routes:
            self.logger.warning("No routes to cluster")
            return {}

        # Compute fingerprints for routes
        fingerprints = self._compute_route_fingerprints(self.filtered_routes)

        # Cluster routes using DBSCAN
        clustering = DBSCAN(
            eps=self.clustering_eps,
            min_samples=self.clustering_min_samples,
            metric='jaccard'
        ).fit(fingerprints)

        # Group routes by cluster
        clusters = defaultdict(list)
        for i, label in enumerate(clustering.labels_):
            clusters[int(label)].append(self.filtered_routes[i])

        self.clusters = dict(clusters)

        self.logger.info(f"Clustered routes into {len(self.clusters)} clusters")
        return self.clusters

    def assign_reactions_to_smirks(self, strategy: str = "best_overall", validate: bool = True) -> Dict[int, List]:
        """Assign reactions in routes to SMIRKS patterns and filter out routes with unassignable reactions."""
        if not hasattr(self, 'clusters'):
            self.logger.error("No clustered routes available. Call cluster_routes first.")
            return {}

        if not self.all_reactions:
            self.logger.warning("SMIRKS library is empty. All routes will be discarded.")
            return {}

        valid_clusters = defaultdict(list)

        for cluster_id, routes in self.clusters.items():
            for route in routes:
                route_valid = True
                route['reaction_mappings'] = {}

                for reaction in route['reaction_tree'].reactions():
                    reaction_name: str = self._match_reaction_to_smirks(reaction, strategy, validate)

                    if reaction_name is None:
                        route_valid = False
                        break

                    route['reaction_mappings'][reaction.reaction_smiles()] = reaction_name

                if route_valid:
                    valid_clusters[cluster_id].append(route)

        self.valid_clusters = dict(valid_clusters)

        total_routes = sum(len(routes) for routes in self.clusters.values())
        valid_routes = sum(len(routes) for routes in self.valid_clusters.values())

        self.logger.info(
            f"Assigned reactions using {strategy} strategy: {valid_routes}/{total_routes} routes are valid")

        return self.valid_clusters

    def get_best_routes(self, max_routes_per_cluster: int = 1) -> List:
        """Get the best routes from each cluster."""
        if not hasattr(self, 'valid_clusters'):
            self.logger.error("No valid routes availablse. Call assign_reactions_to_smirks first.")
            return []

        if not self.valid_clusters:
            self.logger.info("No valid clusters found. Returning empty list.")
            return []

        best_routes = []

        for cluster_id, routes in self.valid_clusters.items():
            sorted_routes = sorted(routes, key=lambda r: r['score']['state score'], reverse=True)
            best_routes.extend(sorted_routes[:max_routes_per_cluster])

        best_routes = self._format_routes_as_dict(best_routes)

        self.logger.info(f"Selected {len(best_routes)} best routes from {len(self.valid_clusters)} clusters")
        return best_routes

    def _format_routes_as_dict(self, routes: List[dict]) -> List[dict]:
        """Format routes to return as a dictionary"""
        formatted_routes = []
        for route in routes:
            reactions = []
            for rxn in route['reaction_tree'].reactions():
                reactions.append({'name': route['reaction_mappings'][rxn.reaction_smiles()],
                                  'reactantSmiles': tuple([umol.smiles for umol in rxn.reactants[0]]),
                                  'productSmiles': rxn.mol.smiles})
            formatted_routes.append({
                'reactions': reactions
            })
        return formatted_routes

    def perform_route_search(self,
                             target_smiles: str,
                             matching_strategy: str = "best_child",
                             validate_matches: bool = True,
                             top_n: int = 20,
                             max_routes_per_cluster: int = 1) -> List:
        """Perform full route search with the new SMIRKS library structure."""
        self.logger.info(f"Starting route search for {target_smiles}")
        self.logger.info(f"Using {len(self.all_reactions)} total reactions ({matching_strategy} strategy)")

        if not self.find_routes(target_smiles, top_n):
            raise NoSynthesisRoute(message='No synthesis route found by AiZynthFinder.',
                                   smiles=target_smiles)

        self.filter_routes()
        self.cluster_routes()  # Why is this done??
        self.assign_reactions_to_smirks(matching_strategy, validate_matches)
        best_routes: List = self.get_best_routes(max_routes_per_cluster)

        self.logger.info(f"Completed pipeline for {target_smiles} using {matching_strategy} strategy. "
                         f"Found {len(best_routes)} best routes.")
        return best_routes

    def get_library_info(self) -> Dict:
        """Get detailed information about the current library."""
        stats = self.smirks_manager.get_library_stats()
        relationships = self.smirks_manager.get_parent_child_relationships()

        return {
            **stats,
            "example_reactions": list(self.all_reactions.keys())[:10],
            "parent_child_examples": dict(list(relationships.items())[:5])
        }

    def export_routes_to_dict(self, routes: List[dict], output_path: str = None) -> list:
        """Export routes to dict and can save to .json if specified."""
        if output_path:
            assert '.json' in output_path, self.logger.error(
                f"Output file {output_path} does not contain .json extension")
        try:
            stored_routes = []
            for i, route in enumerate(routes):
                rxn_tree = route['reaction_tree']
                rxn_names = []
                rxn_metadata = []
                for reaction in rxn_tree.reactions():
                    rxn_names.append(reaction.smiles)
                    rxn_metadata.append(reaction.metadata)
                route_data = {
                    'smiles': self.finder.target_smiles,
                    'route_id': i,
                    'score': route['score']['state score'],
                    'num_reactions': len(list(rxn_tree.reactions())),
                    'reactions': rxn_names,
                    'metadata': rxn_metadata,
                    'smirks_mappings': '|'.join([
                        f"{r.reaction_smiles()}:{route.reaction_mappings.get(r.reaction_smiles(), 'N/A')}"
                        for r in rxn_tree.reactions()
                    ]) if hasattr(route, 'reaction_mappings') else 'N/A'
                }
                stored_routes.append(route_data)
            if output_path:
                with open(output_path, 'w') as json_file:
                    json.dump(stored_routes, json_file, indent=4)
            return stored_routes
        except Exception as e:
            self.logger.error(f"Failed to export routes to dict: {e}")
            return []
