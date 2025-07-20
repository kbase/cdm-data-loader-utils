"""
Unit tests for RAST to SEED role mapper
"""

import pytest
import json
import os
from pathlib import Path
import sys

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from utils.rast_seed_mapper import RASTSeedMapper, map_rast_to_seed, map_rast_batch


@pytest.fixture
def test_ontology_path():
    """Path to test SEED ontology file"""
    return Path(__file__).parent.parent / "src" / "data" / "seed_ontology.json"


@pytest.fixture
def mapper(test_ontology_path):
    """Create a mapper instance for testing"""
    if test_ontology_path.with_suffix('.json.gz').exists():
        # If only compressed version exists, decompress it first
        import gzip
        import shutil
        gz_path = test_ontology_path.with_suffix('.json.gz')
        with gzip.open(gz_path, 'rb') as f_in:
            with open(test_ontology_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
    
    return RASTSeedMapper(test_ontology_path)


@pytest.fixture
def example_annotations():
    """Load example annotations"""
    example_file = Path(__file__).parent.parent / "src" / "data" / "example_rast_annotations.json"
    with open(example_file, 'r') as f:
        data = json.load(f)
    
    # Flatten all examples
    all_annotations = []
    for category in data['annotations']:
        all_annotations.extend(category['examples'])
    
    return all_annotations, data.get('expected_mappings', {})


class TestRASTSeedMapper:
    """Test cases for RAST to SEED mapper"""
    
    def test_initialization(self, test_ontology_path):
        """Test mapper initialization"""
        mapper = RASTSeedMapper(test_ontology_path)
        assert len(mapper.seed_mapping) > 0
        assert isinstance(mapper.seed_mapping, dict)
    
    def test_invalid_file(self):
        """Test initialization with invalid file"""
        with pytest.raises(FileNotFoundError):
            RASTSeedMapper("nonexistent_file.json")
    
    def test_simple_annotations(self, mapper):
        """Test mapping of simple annotations"""
        test_cases = [
            ("Alpha-ketoglutarate permease", "seed.role:0000000010501"),
            ("Thioredoxin 2", "seed.role:0000000049856"),
            ("Unknown function", "seed.role:0000000031207"),
        ]
        
        for annotation, expected in test_cases:
            result = mapper.map_annotation(annotation)
            assert result == expected, f"Failed to map '{annotation}'"
    
    def test_multi_function_slash(self, mapper):
        """Test multi-function annotations with / separator"""
        annotation = "GMP synthase [glutamine-hydrolyzing], amidotransferase subunit (EC 6.3.5.2) / GMP synthase [glutamine-hydrolyzing], ATP pyrophosphatase subunit (EC 6.3.5.2)"
        result = mapper.map_annotation(annotation)
        assert result == "seed.role:0000000002981"
    
    def test_multi_function_at(self, mapper):
        """Test multi-function annotations with @ separator"""
        annotation = "Uracil permease @ Uracil:proton symporter UraA"
        result = mapper.map_annotation(annotation)
        assert result == "seed.role:0000000002527"
    
    def test_multi_function_semicolon(self, mapper):
        """Test multi-function annotations with ; separator"""
        annotation = "Lead, cadmium, zinc and mercury transporting ATPase (EC 3.6.3.3) (EC 3.6.3.5); Copper-translocating P-type ATPase (EC 3.6.3.4)"
        result = mapper.map_annotation(annotation)
        assert result == "seed.role:0000000010083"
    
    def test_split_multi_function(self, mapper):
        """Test splitting of multi-function annotations"""
        test_cases = [
            ("A / B / C", ["A", "B", "C"]),
            ("A @ B", ["A", "B"]),
            ("A; B", ["A", "B"]),
            ("A / B @ C", ["A", "B", "C"]),
            ("Single function", ["Single function"]),
        ]
        
        for annotation, expected_parts in test_cases:
            parts = mapper.split_multi_function(annotation)
            assert parts == expected_parts
    
    def test_edge_cases(self, mapper):
        """Test edge cases"""
        # Empty string
        assert mapper.map_annotation("") is None
        assert mapper.map_annotation(None) is None
        
        # Whitespace
        assert mapper.map_annotation("   ") is None
        
        # Unknown annotations
        assert mapper.map_annotation("This does not exist") is None
    
    def test_batch_mapping(self, mapper):
        """Test batch mapping functionality"""
        annotations = [
            "Alpha-ketoglutarate permease",
            "Unknown function",
            "Nonexistent annotation",
            "Thioredoxin 2"
        ]
        
        results = mapper.map_annotations(annotations)
        
        assert len(results) == 4
        assert results[0] == ("Alpha-ketoglutarate permease", "seed.role:0000000010501")
        assert results[1] == ("Unknown function", "seed.role:0000000031207")
        assert results[2] == ("Nonexistent annotation", None)
        assert results[3] == ("Thioredoxin 2", "seed.role:0000000049856")
    
    def test_stats(self, mapper):
        """Test statistics calculation"""
        annotations = [
            "Alpha-ketoglutarate permease",
            "Unknown function",
            "This does not exist",
            "Thioredoxin 2"
        ]
        
        stats = mapper.get_stats(annotations)
        
        assert stats['total'] == 4
        assert stats['mapped'] == 3
        assert stats['unmapped'] == 1
        assert stats['coverage'] == 75.0
        assert len(stats['unmapped_examples']) == 1
        assert "This does not exist" in stats['unmapped_examples']
    
    def test_convenience_functions(self, test_ontology_path):
        """Test standalone convenience functions"""
        # Single mapping
        result = map_rast_to_seed("Alpha-ketoglutarate permease", test_ontology_path)
        assert result == "seed.role:0000000010501"
        
        # Batch mapping
        annotations = ["Alpha-ketoglutarate permease", "Unknown function"]
        results = map_rast_batch(annotations, test_ontology_path)
        assert len(results) == 2
        assert results[0][1] == "seed.role:0000000010501"
        assert results[1][1] == "seed.role:0000000031207"
    
    def test_all_examples(self, mapper, example_annotations):
        """Test all example annotations"""
        all_annotations, expected_mappings = example_annotations
        
        # Map all annotations
        results = mapper.map_annotations(all_annotations)
        
        # Check expected mappings
        for annotation, expected_id in expected_mappings.items():
            result = mapper.map_annotation(annotation)
            assert result == expected_id, f"Failed to map '{annotation}' to '{expected_id}', got '{result}'"
        
        # Calculate overall stats
        stats = mapper.get_stats(all_annotations)
        print(f"\nExample annotations coverage: {stats['coverage']:.2f}%")
        print(f"Mapped: {stats['mapped']}/{stats['total']}")
        
        if stats['unmapped_examples']:
            print("\nUnmapped examples:")
            for ex in stats['unmapped_examples']:
                if ex.strip():  # Skip empty strings
                    print(f"  - {ex}")


class TestIDParsing:
    """Test SEED role ID parsing functionality"""
    
    def test_parse_url_format(self, mapper):
        """Test parsing of URL-based IDs"""
        url_id = "https://pubseed.theseed.org/RoleEditor.cgi?page=ShowRole&Role=0000000001234"
        parsed = mapper._parse_seed_role_id(url_id)
        assert parsed == "seed.role:0000000001234"
    
    def test_parse_clean_format(self, mapper):
        """Test parsing of clean IDs"""
        clean_id = "seed.role:0000000001234"
        parsed = mapper._parse_seed_role_id(clean_id)
        assert parsed == "seed.role:0000000001234"
    
    def test_parse_invalid_format(self, mapper):
        """Test parsing of invalid formats"""
        assert mapper._parse_seed_role_id("") is None
        assert mapper._parse_seed_role_id(None) is None
        assert mapper._parse_seed_role_id("invalid_format") is None


if __name__ == "__main__":
    pytest.main([__file__, "-v"])