# SEED Ontology Data

This directory contains the SEED role ontology used for mapping RAST annotations to seed.role identifiers.

## Current File

- **`seed_ontology.json.gz`** - Compressed SEED role ontology in JSON-LD format

## File Generation

The current `seed_ontology.json.gz` file was generated from `seed.obo` using the following process:

1. **Source**: The original `seed.obo` file was obtained from the SEED/ModelSEED project
2. **Conversion**: The OBO file was converted to JSON using ROBOT (http://robot.obolibrary.org/):
   ```bash
   robot convert --input seed.obo --output seed.json
   ```
3. **Compression**: The resulting JSON file was compressed using gzip to reduce repository size:
   ```bash
   gzip -9 seed.json
   ```

The OBO format includes idspace definitions like:
```
idspace: seed.role https://pubseed.theseed.org/RoleEditor.cgi?page=ShowRole&Role=
```

This results in JSON nodes with URL-based IDs that the mapper handles automatically.

## Future Plans

### Official Source
We are currently working to identify the official, versioned source for SEED ontology files. Once established, this file will be updated from:
- Official SEED OWL file location (TBD)
- Official SEED OBO file location (TBD)

### Automated Updates
Future versions will implement:
1. Automated fetching from the official source
2. Version tracking and changelog
3. Conversion pipeline from OWL/OBO to JSON as needed
4. Regular updates synchronized with SEED releases

### Versioning
When the official source is established, we will:
- Track the SEED ontology version in the filename (e.g., `seed_ontology_v2024.1.json.gz`)
- Maintain a version history file
- Document any custom modifications or additions

## File Format

The JSON file follows the JSON-LD format with this structure:
```json
{
  "graphs": [{
    "nodes": [
      {
        "id": "https://pubseed.theseed.org/RoleEditor.cgi?page=ShowRole&Role=0000000010501",
        "lbl": "Alpha-ketoglutarate permease",
        "type": "CLASS"
      }
    ]
  }]
}
```

The mapper automatically extracts the role number from the URL to create clean `seed.role:XXXXXXXXXX` identifiers.

## Updating the Ontology

To update the ontology file when an official source becomes available:

1. Download the latest OWL or OBO file from the official source
2. Convert to JSON using ROBOT:
   ```bash
   # From OBO
   robot convert --input seed.obo --output seed_ontology.json
   
   # From OWL
   robot convert --input seed.owl --output seed_ontology.json
   ```
3. Compress the file:
   ```bash
   gzip -9 seed_ontology.json
   ```
4. Replace the existing file in this directory
5. Update this README with the new version information
6. Run tests to ensure compatibility:
   ```bash
   pytest tests/test_rast_seed_mapper.py
   ```

## Notes

- The compressed file is automatically decompressed on first use by the mapper
- The JSON format was chosen for fast parsing and broad compatibility
- The mapper supports both URL-based and clean ID formats for flexibility