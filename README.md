# Identifying natural ligands for the WT1 protein


## Project Overview

This project aims to identify potential natural ligands for the WT1 (hereditary Wilms' tumor 1) protein by filtering and analyzing molecules from the COCONUT database. The project includes several steps:
1. **Data Mining:** Collecting and processing molecular data.
2. **Molecular Dynamics (MD) Modeling:** Simulating molecular behavior in a biological context.
3. **Docking:** Predicting the preferred orientation of molecules when bound to the WT1 protein.

## Installation and Setup

### Prerequisites

- **Python 3.8.19**
- **Conda** (for managing dependencies)


### Setting Up the Environment

1. Clone the repository:

``` bash
git clone https://github.com/MolkaAttia-git/identifying-natural-ligands-for-the-WT1-protein.git
cd identifying-natural-ligands-for-the-WT1-protein
```

2. Create a Conda environment using the \`.yml\` file:

```bash
conda env create -f environment.yml
conda activate WT1
```

### Dependencies

The dependencies for this project are listed in the \`environment.yml\` file:

```yaml
name: WT1_Project
channels: 
 - anaconda
 - conda-forge
 - defaults
dependencies :
- python=3.8.19 
- pandas=2.0.3 
- rdkit=2024.03.5
- numpy=1.24.4 
```

### Filtering Molecules

The script \`WT1_inhibators.py\` filters molecules based on Lipinski's Rule of Five and other drug-likeness criteria. It also calculates the similarity of each molecule to known WT1 inhibitors (Shikonin and Trichostatin A) using MACCS Keys fingerprints.

```bash
python src/WT1_inhibators.py
```

### Output

- **Filtered Molecules (SDF format)**: \`filtered_molecules_with_similarity.sdf\`
- **Molecule Similarity Data (CSV format)**: \`filtered_molecule_data.csv\`

