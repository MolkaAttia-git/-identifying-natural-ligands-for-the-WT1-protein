# Identifying natural ligands for the WT1 protein

## I. Project Overview
This project is focused on identify natural ligands for the WT1 (hereditary Wilms' tumor 1) protein

## II.Getting Started
### 1.Prerequisites
- Python 3.9.19
  
To run the scripts and notebooks in this repository, you'll need to have the following Python packages installed:

- Pandas 2.2.2
- Rdkit 2024.3.5

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

 We will use Shikonin and TSA since they are small molecules 
 
| Inhibitor | Link |
| ------ | ----------- |
| Shikonin   | [path to data files to supply the data that will be passed into templates.](https://pubmed.ncbi.nlm.nih.gov/36500358/ ) |
| Trichostatin A (TSA) | [engine to be used for processing templates. Handlebars is the default.](https://pubmed.ncbi.nlm.nih.gov/18535006/) |

We tested 3 fingerprint methods ; 

- **MACCS Keys**
- **RDKit Fingerprints**:
- **Morgan Fingerprints**

on the first 500 molecules of the coconut db, and with a similarity threshold of 0.7 , the results were the following ; 

Number of similar molecules found (Morgan): 0
Number of similar molecules found (RDKit): 0
Number of similar molecules found (MACCS): 2

⇒ we will be using the MACCS fingerprint

## IV.DATA mining 

### 1.the criteria that determine the drug-likeness of the ligands :

- **Lipinski Rule of 5 (Ro5)**, also known as Lipinski Alert Index : 
predicts that poor absorption or permeation is more likely when:

▪	there are more than 5 H-bond donors (nHDon)
▪	there are more than 10 H-bond acceptors (N + O)
▪	molecular weight (MW) is over 500
▪	Moriguchi's logP (MLogP) is over 4.15

- **NP-likeness score : 
the likelihood of the compound to be a natural product,**
- Alogp
- Alogp2
- Fsp3

After applying the code we have :
The Number of molecules that passed the filters: `3050`
