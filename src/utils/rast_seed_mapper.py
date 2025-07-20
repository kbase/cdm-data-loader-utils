"""
RAST to SEED Role Mapper

This module provides utilities for mapping RAST (Rapid Annotation using Subsystem Technology)
annotations to SEED role ontology identifiers.

The mapper handles complex multi-function annotations and provides multiple fallback strategies
to achieve comprehensive mapping coverage.

Example usage:
    from rast_seed_mapper import RASTSeedMapper, map_rast_to_seed
    
    # Use default bundled ontology
    mapper = RASTSeedMapper()
    
    # Or specify custom ontology file
    mapper = RASTSeedMapper("path/to/seed.json")
    
    # Single annotation mapping
    seed_id = map_rast_to_seed("Alpha-ketoglutarate permease")
    
    # Batch processing
    results = mapper.map_annotations(["annotation1", "annotation2"])
"""

import json
import logging
import gzip
import shutil
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Union

__version__ = "1.0.0"
__author__ = "KBase CDM Team"

logger = logging.getLogger(__name__)

# Default path to bundled SEED ontology
DEFAULT_ONTOLOGY_PATH = Path(__file__).parent.parent / "data" / "seed_ontology.json"


class RASTSeedMapper:
    """
    Maps RAST annotations to SEED role ontology identifiers.
    
    This mapper handles:
    - Direct exact matches
    - Multi-function annotations with various separators (/, @, ;)
    - Different SEED ontology formats (URL-based and clean IDs)
    
    Attributes:
        seed_mapping (Dict[str, str]): Internal mapping from annotations to seed.role IDs
        multi_func_separators (List[str]): Separators used in multi-function annotations
    """
    
    def __init__(self, seed_ontology_path: Optional[Union[str, Path]] = None):
        """
        Initialize the mapper with a SEED ontology file.
        
        Args:
            seed_ontology_path: Path to SEED ontology JSON file. If None, uses bundled ontology.
            
        Raises:
            FileNotFoundError: If the ontology file doesn't exist
            ValueError: If the file format is not supported
            json.JSONDecodeError: If the JSON file is malformed
        """
        self.seed_mapping: Dict[str, str] = {}
        self.multi_func_separators = [' / ', ' @ ', '; ']
        
        # Use default path if none provided
        if seed_ontology_path is None:
            path = DEFAULT_ONTOLOGY_PATH
            # Check for compressed version if uncompressed doesn't exist
            if not path.exists() and path.with_suffix('.json.gz').exists():
                self._decompress_ontology(path.with_suffix('.json.gz'), path)
        else:
            path = Path(seed_ontology_path)
        
        if not path.exists():
            raise FileNotFoundError(f"Ontology file not found: {path}")
            
        if path.suffix not in ['.json', '.gz']:
            raise ValueError(f"Unsupported file format: {path.suffix}. Expected .json or .json.gz")
            
        self._load_seed_ontology(path)
    
    def _decompress_ontology(self, gz_path: Path, json_path: Path) -> None:
        """Decompress gzipped ontology file"""
        logger.info(f"Decompressing ontology from {gz_path}")
        with gzip.open(gz_path, 'rb') as f_in:
            with open(json_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        
    def _load_seed_ontology(self, path: Path) -> None:
        """Load SEED ontology from JSON file."""
        try:
            with open(path, 'r', encoding='utf-8') as f:
                data = json.load(f)
        except json.JSONDecodeError as e:
            logger.error(f"Failed to parse JSON file: {e}")
            raise
        
        # Extract nodes from JSON-LD format
        graphs = data.get("graphs", [])
        if not graphs:
            logger.warning("No graphs found in ontology file")
            return
            
        nodes = graphs[0].get("nodes", [])
        
        for node in nodes:
            label = node.get("lbl")
            node_id = node.get("id")
            
            if not label or not node_id:
                continue
                
            # Parse different ID formats
            seed_role_id = self._parse_seed_role_id(node_id)
            if seed_role_id:
                self.seed_mapping[label] = seed_role_id
                
        logger.info(f"Loaded {len(self.seed_mapping)} SEED role mappings")
        
    def _parse_seed_role_id(self, raw_id: str) -> Optional[str]:
        """
        Parse SEED role ID from various formats.
        
        Handles:
        - URL format: https://pubseed.theseed.org/RoleEditor.cgi?page=ShowRole&Role=0000000001563
        - Clean format: seed.role:0000000001563
        - Other formats
        
        Args:
            raw_id: Raw ID string from ontology
            
        Returns:
            Cleaned seed.role ID or None if parsing fails
        """
        if not raw_id:
            return None
            
        # URL format with Role parameter
        if "Role=" in raw_id:
            try:
                role_number = raw_id.split("Role=")[-1]
                return f"seed.role:{role_number}"
            except IndexError:
                logger.warning(f"Failed to parse role number from URL: {raw_id}")
                return None
                
        # Already in clean format
        if raw_id.startswith("seed.role:"):
            return raw_id
            
        # Try generic parsing
        try:
            ontology_part = raw_id.split('/')[-1]
            if '_' in ontology_part:
                return ontology_part.replace("_", ":", 1)
        except Exception as e:
            logger.warning(f"Failed to parse ID: {raw_id}, error: {e}")
            
        return None
    
    def split_multi_function(self, annotation: str) -> List[str]:
        """
        Split multi-function annotations into individual components.
        
        Args:
            annotation: RAST annotation string
            
        Returns:
            List of individual function annotations
        """
        if not annotation:
            return []
            
        # Start with the full annotation
        parts = [annotation]
        
        # Split by each separator in order
        for separator in self.multi_func_separators:
            new_parts = []
            for part in parts:
                split_parts = part.split(separator)
                new_parts.extend(p.strip() for p in split_parts if p.strip())
            parts = new_parts
            
        return parts
    
    def map_annotation(self, annotation: str) -> Optional[str]:
        """
        Map a RAST annotation to its SEED role ID.
        
        Tries multiple strategies:
        1. Direct exact match
        2. Individual parts of multi-function annotations
        
        Args:
            annotation: RAST annotation string
            
        Returns:
            seed.role ID if found, None otherwise
        """
        if not annotation:
            return None
            
        # Try direct match first
        if annotation in self.seed_mapping:
            return self.seed_mapping[annotation]
            
        # Try splitting multi-function annotations
        parts = self.split_multi_function(annotation)
        
        if len(parts) > 1:
            # Try each part
            for part in parts:
                if part in self.seed_mapping:
                    return self.seed_mapping[part]
                    
        return None
    
    def map_annotations(self, annotations: List[str]) -> List[Tuple[str, Optional[str]]]:
        """
        Map multiple RAST annotations to SEED role IDs.
        
        Args:
            annotations: List of RAST annotation strings
            
        Returns:
            List of tuples (annotation, seed_role_id or None)
        """
        return [(ann, self.map_annotation(ann)) for ann in annotations]
    
    def get_stats(self, annotations: List[str]) -> Dict[str, Union[int, float, List[str]]]:
        """
        Calculate mapping statistics for a set of annotations.
        
        Args:
            annotations: List of RAST annotation strings
            
        Returns:
            Dictionary containing:
            - total: Total number of annotations
            - mapped: Number successfully mapped
            - unmapped: Number not mapped
            - coverage: Percentage of annotations mapped
            - unmapped_examples: First 10 unmapped annotations
        """
        results = self.map_annotations(annotations)
        
        total = len(annotations)
        mapped = sum(1 for _, seed_id in results if seed_id is not None)
        unmapped = total - mapped
        
        unmapped_annotations = [ann for ann, seed_id in results if seed_id is None]
        
        return {
            'total': total,
            'mapped': mapped,
            'unmapped': unmapped,
            'coverage': (mapped / total * 100) if total > 0 else 0.0,
            'unmapped_examples': unmapped_annotations[:10]
        }


# Convenience functions

def map_rast_to_seed(annotation: str, seed_ontology_path: Optional[Union[str, Path]] = None) -> Optional[str]:
    """
    Map a single RAST annotation to a SEED role ID.
    
    This is a convenience function for one-off mappings. For batch processing,
    use the RASTSeedMapper class directly to avoid reloading the ontology.
    
    Args:
        annotation: RAST annotation string
        seed_ontology_path: Path to SEED ontology JSON file. If None, uses bundled ontology.
        
    Returns:
        seed.role ID if found, None otherwise
        
    Example:
        >>> seed_id = map_rast_to_seed("Alpha-ketoglutarate permease")
        >>> print(seed_id)
        'seed.role:0000000010501'
    """
    mapper = RASTSeedMapper(seed_ontology_path)
    return mapper.map_annotation(annotation)


def map_rast_batch(annotations: List[str], 
                   seed_ontology_path: Optional[Union[str, Path]] = None) -> List[Tuple[str, Optional[str]]]:
    """
    Map a batch of RAST annotations to SEED role IDs.
    
    Args:
        annotations: List of RAST annotation strings
        seed_ontology_path: Path to SEED ontology JSON file. If None, uses bundled ontology.
        
    Returns:
        List of tuples (annotation, seed_role_id or None)
        
    Example:
        >>> annotations = ["Alpha-ketoglutarate permease", "Unknown function"]
        >>> results = map_rast_batch(annotations)
        >>> for ann, seed_id in results:
        ...     print(f"{ann} -> {seed_id}")
    """
    mapper = RASTSeedMapper(seed_ontology_path)
    return mapper.map_annotations(annotations)