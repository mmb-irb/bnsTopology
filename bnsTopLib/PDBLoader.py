
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBList import PDBList

class PDBLoader():
    def __init__(self, pdb_path):
        if "pdb:"in pdb_path:
            pdbl= PDBList(pdb='tmpPDB')
            try:
                pdb_path = pdb_path[4:].upper()
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
        return st
