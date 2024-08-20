# Identifying natural ligands for the WT1 protein

## I. Project Overview
This project is focused on identify natural ligands for the WT1 (hereditary Wilms' tumor 1) protein

## II.Getting Started
### 1.Prerequisites
- Python 3.8.19
  
To run the scripts and notebooks in this repository, you'll need to have the following Python packages installed:

- Pandas  2.0.3 
- Rdkit 2024.03.5

You can install the required packages via pip:
```
pip install rdkit pandas
pip install rdkit Rdkit
```
## 2.Instructions
First we install Anaconda to create a virtual environment called `WT1_inhibitors` by using the folowing command :
```
conda create -n WT1_inhibitors
```
To activate the environment we use the folowing command :
```
conda activate WT1_inhibitors
```
### 3.Importing libraries
```
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors , Lipinski
import pandas as pd
import numpy as np
```
### III. Identifying the known WT1 inhibitors : 

 We will use Shikonin and TSA 


