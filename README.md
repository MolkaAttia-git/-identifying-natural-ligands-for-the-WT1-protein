# Identifying natural ligands for the WT1 protein
In this project our goal is to identify  natural ligands for the WT1 (hereditary Wilms' tumor 1) protein

## 1. Downloading data
Download Natural Products Structures in SDF formats from the COCONUT database (Collection of Open Natural Products) dedicated to natural products.

## Installing libraries
- Python 3.9.19 
- Numpy 2.0.1 
- Pandas 2.2.2
- Rdkit 2024.3.5

## Instructions
First we install Anaconda to create a virtual environment called `WT1_inhibitors` by using the folowing command :
```
conda create -n WT1_inhibitors
```
To activate the environment we use the folowing command :
```
conda activate WT1_inhibitors
```
### Importing libraries
```
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors , Lipinski
import pandas as pd
import numpy as np
```
