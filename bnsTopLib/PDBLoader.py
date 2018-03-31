
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBList import PDBList

import sys

class PDBLoader():
    def __init__(self, pdb_path):
        self.useChains=False
        if "pdb:"in pdb_path:
            pdbl= PDBList(pdb='tmpPDB')
            try:
                pdb_path = pdb_path[4:].upper()
                self.id=pdb_path
                self.pdb_path=pdbl.retrieve_pdb_file(pdb_path)
                self.parser = MMCIFParser()
                self.useChains=True
                self.format='cif'
            except IOError:
                print ("#ERROR: fetching PDB "+pdb_path)
                sys.exit(2)
        else:
            self.pdb_path = pdb_path
            if '.pdb' in pdb_path:
                self.parser = PDBParser(PERMISSIVE=1)
                self.format='pdb'
            elif '.cif' in pdb_path:
                self.parser = MMCIFParser()
                self.format='cif'
            else:
                print ('#ERROR: unknown filetype')
                sys.exit(2)

    def loadStructure(self):
        try:
            st = self.parser.get_structure('st', self.pdb_path)
        except OSError:
            print ("#ERROR: parsing PDB")
            sys.exit(2)
        #====== Internal residue renumbering =========================================
        i=1
        for r in st.get_residues():
            r.index = i
            i=i+1
        #Atom renumbering for mmCIF, 
        if self.format == 'cif':
            i=1
            for at in st.get_atoms():
                at.serial_number = i
                if hasattr(at,'selected_child'):
                    at.selected_child.serial_number=i
                i=i+1
# Checking for models
        if len(st) > 1:
            print ("#WARNING: Several Models found, using only first")
# Using Model 0 any way TODO: revise for bioiunits
        st = st[0]
        self.numAts=0
        for at in st.get_atoms():
            self.numAts=self.numAts+1
	


        return st
