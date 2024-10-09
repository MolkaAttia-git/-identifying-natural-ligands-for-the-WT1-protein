from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski, DataStructs, MACCSkeys
import csv

# Drug-likeness Filters

def lipinski_filter(mol):
    mol_wt = Descriptors.MolWt(mol)
    mol_logp = Crippen.MolLogP(mol)
    num_h_donors = Lipinski.NumHDonors(mol)
    num_h_acceptors = Lipinski.NumHAcceptors(mol)
    return (mol_wt <= 500 and mol_logp <= 5 and
            num_h_donors <= 5 and num_h_acceptors <= 10)

def filter_by_h_bond_donors(mol, min_donors=0, max_donors=5):
    num_h_donors = Lipinski.NumHDonors(mol)
    return min_donors <= num_h_donors <= max_donors

def filter_by_h_bond_acceptors(mol, min_acceptors=0, max_acceptors=10):
    num_h_acceptors = Lipinski.NumHAcceptors(mol)
    return min_acceptors <= num_h_acceptors <= max_acceptors

def filter_by_alogp(mol, min_alogp=-3, max_alogp=5):
    alogp = Crippen.MolLogP(mol)
    return min_alogp <= alogp <= max_alogp

def filter_by_rotatable_bonds(mol, max_rotatable_bonds=10):
    num_rotatable_bonds = Lipinski.NumRotatableBonds(mol)
    return num_rotatable_bonds <= max_rotatable_bonds

def filter_by_tpsa(mol, max_tpsa=140):
    tpsa = Descriptors.TPSA(mol)
    return tpsa <= max_tpsa

# Load the SDF files for JDQ443 and Trichostatin A
JDQ443_mol = Chem.SDMolSupplier('"D:\\KRAS inhibitors\\JDQ443.sdf"')[0]
Divarasib_mol = Chem.SDMolSupplier('"D:\\KRAS inhibitors\\divarasib.sdf"')[0]

# Calculate MACCS Keys fingerprints
JDQ443_fp = MACCSkeys.GenMACCSKeys(JDQ443_mol)
Divarasib_fp = MACCSkeys.GenMACCSKeys(Divarasib_mol)

# Load the COCONUT database molecules
sdf_file_path = 'D:\\coconut-08-2024.sdf\\coconut-08-2024.sdf'
coconut_supplier = Chem.SDMolSupplier(sdf_file_path)

# List to store similar molecules and their data
similar_molecules = []
similar_molecule_data = []

# Threshold for similarity
similarity_threshold = 0.7

for mol in coconut_supplier:
    if mol is None:
        continue
    
    # Calculate MACCS Keys fingerprints for the current molecule
    mol_fp = MACCSkeys.GenMACCSKeys(mol)
    
    # Calculate similarity with both JDQ443 and Divarasib
    JDQ443_similarity = DataStructs.TanimotoSimilarity(JDQ443_fp, mol_fp)
    Divarasib_similarity = DataStructs.TanimotoSimilarity(Divarasib_fp, mol_fp)
    
    # Check if similarity exceeds the threshold
    if JDQ443_similarity >= similarity_threshold or Divarasib_similarity >= similarity_threshold:
        similar_molecules.append(mol)

# Apply the drug-likeness filters and generate SMILES
filtered_molecules = []

for mol in similar_molecules:
    if lipinski_filter(mol):
        if (filter_by_h_bond_donors(mol) and 
            filter_by_rotatable_bonds(mol) and 
            filter_by_tpsa(mol) and
            filter_by_h_bond_acceptors(mol) and
            filter_by_alogp(mol)):
            filtered_molecules.append(mol)
            
            # Store the COCONUT_ID, similarity measures, and SMILES
            if mol.HasProp('COCONUT_ID'):
                coconut_id = mol.GetProp('COCONUT_ID')
                smiles = Chem.MolToSmiles(mol)  # Generate SMILES
                filtered_molecule_data = [coconut_id, JDQ443_similarity, Divarasib_similarity, smiles]
                similar_molecule_data.append(filtered_molecule_data)

# Write the filtered molecules to a new SDF file
output_sdf_path = 'D:/KRAS_filtered_molecules_with_similarity.sdf'
writer = Chem.SDWriter(output_sdf_path)

for mol in filtered_molecules:
    writer.write(mol)

writer.close()

# Write the data to a CSV file with SMILES
output_csv_path = 'D:/filtered_molecule_data.csv'
with open(output_csv_path, mode='w', newline='') as file:
    writer = csv.writer(file, delimiter=',')
    writer.writerow(['COCONUT_ID', 'JDQ443_Similarity', 'Divarasib_Similarity', 'SMILES'])
    writer.writerows(similar_molecule_data)
    
# Print the number of molecules that passed the filters
print(f"Number of molecules that passed the filters: {len(filtered_molecules)}")
print(f"Filtered molecules saved to {output_sdf_path}")
print(f"Molecule data saved to {output_csv_path}") 
