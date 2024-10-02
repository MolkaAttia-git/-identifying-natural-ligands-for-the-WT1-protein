from Bio import PDB

class CleanSelect(PDB.Select):
    def accept_residue(self, residue):
        # Remove water molecules (HOH) and heteroatoms (HetAtm)
        if residue.id[0] == ' ':
            return True  # Keep only standard amino acids
        return False

    def accept_atom(self, atom):
        # Keep only the most occupied conformation (altloc ' ')
        if atom.altloc == ' ' or atom.altloc == 'A':
            return True
        return False

def clean_pdb(input_file, output_file):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("1d8d", input_file)

    # Create PDBIO object to write the cleaned structure
    io = PDB.PDBIO()
    io.set_structure(structure)

    # Save cleaned structure
    io.save(output_file, CleanSelect())

input_file = "D:\\1d8d.pdb"
output_file = "D:\\cleaned_1d8d.pdb"

# Clean the PDB file
clean_pdb(input_file, output_file)

