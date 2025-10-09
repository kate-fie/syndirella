<div align="center">
  <img src="logos/Full.png" alt="Syndirella Logo" width="400"/>
</div>

Syndirella (Synthesis Directed Elaborations) is a tool for generating and scoring synthetically practical elaborations of molecules designed from fragment screens. Starting from the retrosynthetic analysis of fragment merged molecules, superstructures of the original reactants are found by catalog search and filtered by selectivity issues. The elaborated final products are defined by reaction SMIRKS patterns and energy minimized in the protein with restraints to experimental data.

## Documentation

**[Full Documentation](https://syndirella.readthedocs.io/)** - Complete user guide, API reference, and examples

[![Documentation Status](https://readthedocs.org/projects/syndirella/badge/?version=latest)](https://syndirella.readthedocs.io/en/latest/?badge=latest)

## Quick Start

### Installation

```bash
conda create -n syndirella python=3.10
conda activate syndirella
pip install cgrtools --use-pep517
pip install syndirella
```

### Basic Usage

1. **Setup AiZynthFinder** (required for retrosynthesis):

   ⚠️ WARNING: This step downloads large model files (~750MB total) if run the first time. These are required to run AiZynthFinder. 

   Files will be downloaded to: `[syndirella_package_path]/aizynth/`

   ```bash
   syndirella setup-aizynth
   ```

2. **Run the pipeline**:
   ```bash
   syndirella run -i input.csv -o output_directory -t templates/ --hits_path fragments.sdf --metadata metadata.csv
   ```

### Example Input

Create a CSV file with your fragment data:

```csv
smiles,hit1,hit2,hit3,template,compound_set
O=C(NC1CC(C(F)(F)F)C1)c1cc2ccsc2[nH]1,Ax0556a,Ax0566a,,Ax0310a,my_compound_set
CC(=O)Nc1cc(CC(=O)NCC(NC(=O)CCl)c2cccnc2)c(NC(C)=O)nn1,1346a,,,Ax1346a,my_compound_set
```

## Requirements

- Python 3.10
- RDKit
- AiZynthFinder
- Fragmenstein
- Additional dependencies (see [pyproject.toml](pyproject.toml))

## Learn More

- **[Examples](examples/)** - Jupyter notebooks and sample data
- **[User Guide](https://syndirella.readthedocs.io/en/latest/user_guide_w_aizynth.html)** - Detailed usage instructions
- **[Reaction Constants](https://syndirella.readthedocs.io/en/latest/configuration.html)** - Advanced reaction settings and options

## License

This project is licensed under the MIT License.