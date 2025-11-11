"""
AiZynthManager Module - Enhanced AiZynthFinder integration with structured SMIRKS library support.
"""

import json
import logging
import os
import subprocess
import sys
from collections import defaultdict
from typing import List, Dict, Optional, Tuple
import glob2

import numpy as np
from aizynthfinder.aizynthfinder import AiZynthFinder
from aizynthfinder.chem.reaction import FixedRetroReaction
from rdkit import Chem
from rdkit.Chem import rdChemReactions
from sklearn.cluster import DBSCAN

from syndirella.utils.error import AiZynthFinderError, NoSynthesisRoute, USPTOTemplateValidationError
from syndirella.route.SmirksLibraryManager import SmirksLibraryManager


class AiZynthManager:
    def __init__(self,
                 config_file: str = None,
                 stock: str = 'zinc',
                 expansion_policy: str = 'uspto',
                 filter_policy: str = 'uspto',
                 min_route_score: float = 0.5,
                 clustering_eps: float = 0.2,
                 clustering_min_samples: int = 2,
                 similarity_threshold: float = 0.2,
                 smirks_library_path: str = None,
                 logger: logging.Logger = None,
                 auto_setup: bool = True):
        # Setup logging
        self.logger = logger or logging.getLogger(__name__)

        # Determine config file path
        if config_file is None:
            config_file = os.getenv('AIZYNTH_CONFIG_FILE', 'config.yml')
        
        # Auto-setup AiZynthFinder if requested
        if auto_setup:
            config_file = self._setup_aizynthfinder(config_file)

        self.config_file = config_file

        # Initialize AiZynthFinder
        self.finder = self._set_finder(config_file, stock, expansion_policy, filter_policy)

        # Set parameters
        self.min_route_score = min_route_score
        self.clustering_eps = clustering_eps
        self.clustering_min_samples = clustering_min_samples
        self.similarity_threshold = similarity_threshold

        # Initialize SMIRKS library manager
        self._initialize_smirks_manager(smirks_library_path)

    def _setup_aizynthfinder(self, config_file: str) -> str:
        """
        Automatically setup AiZynthFinder by downloading data and creating config if needed.
        """
        # Find the syndirella installation directory
        try:
            import syndirella
            syndirella_dir = os.path.dirname(syndirella.__file__)
        except ImportError:
            # Fallback to relative path if import fails
            syndirella_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        
        aizynth_dir = os.path.join(syndirella_dir, 'aizynth')
        os.makedirs(aizynth_dir, exist_ok=True)
        
        # Check if data already exists
        aizynth_files = glob2.glob(os.path.join(aizynth_dir, '*.onnx'))
        if not aizynth_files:
            self.logger.info("Downloading AiZynthFinder public data...")
            try:
                # Run download_public_data
                result = subprocess.run(
                    ['download_public_data', aizynth_dir],
                    capture_output=True,
                    text=True,
                    cwd=aizynth_dir
                )
                if result.returncode != 0:
                    self.logger.error(f"Failed to download AiZynthFinder data: {result.stderr}")
                    return config_file
                self.logger.info("AiZynthFinder data downloaded successfully")
            except FileNotFoundError:
                self.logger.error("download_public_data command not found. Please install AiZynthFinder properly.")
                return config_file
        
        # Create config file if it doesn't exist
        config_path = os.path.join(aizynth_dir, 'config.yml')
        if not os.path.exists(config_path):
            self.logger.info("Creating AiZynthFinder config file...")
            self._create_config_file(config_path, aizynth_dir)
        else:
            self.logger.info(f"AiZynthFinder config file found: {config_path}")
        
        # Set environment variable for current session
        os.environ['AIZYNTH_CONFIG_FILE'] = config_path
        
        return config_path

    def _create_config_file(self, config_path: str, aizynth_dir: str):
        """Create a basic AiZynthFinder config file."""
        config_content = f"""expansion:
                        uspto:
                            - {os.path.join(aizynth_dir, 'uspto_model.onnx')}
                            - {os.path.join(aizynth_dir, 'uspto_templates.csv.gz')}
                        ringbreaker:
                            - {os.path.join(aizynth_dir, 'uspto_ringbreaker_model.onnx')}
                            - {os.path.join(aizynth_dir, 'uspto_ringbreaker_templates.csv.gz')}
                        filter:
                        uspto: {os.path.join(aizynth_dir, 'uspto_filter_model.onnx')}
                        stock:
                        zinc: {os.path.join(aizynth_dir, 'zinc_stock.hdf5')}
                        """
        
        with open(config_path, 'w') as f:
            f.write(config_content)

    def _initialize_smirks_manager(self, smirks_library_path: str = None):
        """Initialize the SMIRKS library manager."""
        if smirks_library_path is None:
            # Fix the path construction to correctly locate the constants directory
            base_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                                    'constants')
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
        self.logger.info(f"SMIRKS library loaded: {stats['manual_reactions']} manual + "
                         f"{stats['rxn-insight_reactions']} rxn-insight = {stats['total_reactions']} total reactions")
        self.logger.info(
            f"Parent/Child breakdown: {stats['parent_reactions']} parents, {stats['child_reactions']} children")
        self.logger.info(f"USPTO template mappings: {stats['uspto_templates_total']} codes mapped")

        # Validate that USPTO template lookup is properly loaded
        if stats['uspto_templates_total'] == 0:
            self.logger.warning("USPTO template lookup appears to be empty. This may cause issues with reaction validation.")

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
        total_routes = sum(len(routes) for routes in self.clusters.values())
        routes_with_unassignable_reactions = 0

        for cluster_id, routes in self.clusters.items():
            for route in routes:
                route_valid = True
                route['reaction_mappings'] = {}

                for reaction in route['reaction_tree'].reactions():
                    reaction_name: str = self._match_reaction_to_smirks(reaction, strategy, validate)

                    if reaction_name is None:
                        route_valid = False
                        routes_with_unassignable_reactions += 1
                        break

                    route['reaction_mappings'][reaction.reaction_smiles()] = reaction_name

                if route_valid:
                    valid_clusters[cluster_id].append(route)

        self.valid_clusters = dict(valid_clusters)
        valid_routes = sum(len(routes) for routes in self.valid_clusters.values())

        # Check if we have routes but not all could be assigned to USPTO templates
        if total_routes > 0 and valid_routes == 0:
            raise USPTOTemplateValidationError(
                message=f"Found {total_routes} routes but none could be assigned to USPTO templates.",
                smiles=getattr(self, 'finder', None) and getattr(self.finder, 'target_smiles', None),
                num_routes_found=total_routes,
                num_routes_with_valid_templates=valid_routes
            )

        self.logger.info(
            f"Assigned reactions using {strategy} strategy: {valid_routes}/{total_routes} routes are valid")

        return self.valid_clusters

    def get_best_routes(self, max_routes_per_cluster: int = 1) -> List:
        """Get the best routes from each cluster."""
        if not hasattr(self, 'valid_clusters'):
            self.logger.error("No valid routes available. Call assign_reactions_to_smirks first.")
            return []

        if not self.valid_clusters:
            self.logger.info("No valid clusters found. Returning empty list.")
            return []

        best_routes = []

        for cluster_id, routes in self.valid_clusters.items():
            sorted_routes = sorted(routes, key=lambda r: r['score']['state score'], reverse=True)
            best_routes.extend(sorted_routes[:max_routes_per_cluster])

        # Sort all routes by score overall (not just within clusters) to get true top N
        best_routes = sorted(best_routes, key=lambda r: r['score']['state score'], reverse=True)
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
                             max_routes_per_cluster: int = 1,
                             save_analysis: bool = False,
                             analysis_output_path: str = None) -> List:
        """Perform full route search with the new SMIRKS library structure."""
        self.logger.info(f"Starting route search for {target_smiles}")
        self.logger.info(f"Using {len(self.all_reactions)} total reactions ({matching_strategy} strategy)")

        if not self.find_routes(target_smiles, top_n):
            raise NoSynthesisRoute(message='No synthesis route found by AiZynthFinder.',
                                   smiles=target_smiles)

        self.filter_routes()
        self.cluster_routes()
        self.assign_reactions_to_smirks(matching_strategy, validate_matches)
        
        # Store original routes before formatting for detailed analysis
        original_routes = []
        for cluster_id, routes in self.valid_clusters.items():
            original_routes.extend(routes)
        
        best_routes: List = self.get_best_routes(max_routes_per_cluster)

        # Check if no routes were returned after all processing
        if not best_routes:
            raise NoSynthesisRoute(
                message='No valid synthesis routes found after USPTO template validation.',
                smiles=target_smiles
            )

        # Save route analysis data if requested
        if save_analysis:
            if analysis_output_path is None:
                analysis_output_path = f"{target_smiles}_route_analysis.json"
            self.export_route_analysis_to_json(target_smiles, original_routes, analysis_output_path)

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

    def export_route_analysis_to_json(self, target_smiles: str, original_routes: List[dict], output_path: str = None) -> List[dict]:
        """
        Export detailed route analysis data to JSON format.
        """
        route_analyses = []
        
        for route in original_routes:
            original_route = route if original_routes else None
            
            route_analysis = {
                'input_smiles': target_smiles,
                'original_aizynth_route': None
            }
            
            # Add original route object if available
            if original_route:
                # Convert the original route to a serializable format
                try:
                    if isinstance(original_route, dict):
                        # Create a serializable version of the original route
                        serializable_route = {}
                        for key, value in original_route.items():
                            if key == 'reaction_tree':
                                # Convert ReactionTree to a serializable format
                                try:
                                    reaction_smiles = [rxn.reaction_smiles() for rxn in value.reactions()]
                                    metadata = []
                                    for rxn in value.reactions():
                                        rxn_metadata = {
                                            'reaction_smiles': rxn.reaction_smiles(),
                                            'template_code': rxn.metadata.get('template_code'),
                                            'template': rxn.metadata.get('template'),
                                            'metadata': dict(rxn.metadata)
                                        }
                                        metadata.append(rxn_metadata)
                                    serializable_route[key] = {
                                        'reaction_smiles': reaction_smiles,
                                        'metadata': metadata
                                    }
                                except Exception as e:
                                    self.logger.warning(f"Could not serialize reaction_tree: {e}")
                                    serializable_route[key] = str(value)
                            elif key == 'reaction_mappings':
                                serializable_route[key] = value
                            elif key == 'score':
                                serializable_route[key] = value
                            else:
                                # Try to convert other objects to string
                                try:
                                    serializable_route[key] = str(value)
                                except:
                                    serializable_route[key] = None
                        
                        route_analysis['original_aizynth_route'] = serializable_route
                    else:
                        # Try to convert to dict if it's an object
                        route_analysis['original_aizynth_route'] = original_route.to_dict()
                except Exception as e:
                    self.logger.warning(f"Could not serialize original route: {e}")
            
            route_analyses.append(route_analysis)
        
        # Save to file if output_path is provided
        if output_path:
            try:
                with open(output_path, 'w') as f:
                    json.dump(route_analyses, f, indent=2)
                self.logger.info(f"Route analysis exported to {output_path}")
            except Exception as e:
                self.logger.error(f"Failed to save route analysis to {output_path}: {e}")
        
        return route_analyses
