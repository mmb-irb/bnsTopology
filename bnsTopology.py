#! /usr/bin/python
#
# Simple parser to extact topology of NA from simulation snapshot
# Chain ids or TER not required 
#
__author__ = "gelpi"
__date__ = "$29-ago-2017 16:14:26$"

from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.PDBParser import PDBParser

import os
import sys
import argparse

COVLNK = 2.0
HBLNK  = 3.5
#verbose = True
#debug = True

aa1c = {
'DA' :'A', 'DC' :'C', 'DG' :'G', 'DT' :'T',
 'A' :'A',  'C' :'C',  'G' :'G',  'U' :'U',
'DA3':'A', 'DC3':'C', 'DG3':'G', 'DT3':'T',
 'A3':'A',  'C3':'C',  'G3':'G',  'U3':'U',
'DA5':'A', 'DC5':'C', 'DG5':'G', 'DT5':'T',
 'A5':'A',  'C5':'C',  'G5':'G',  'U5':'U',
'MRA':'A',
'ALA':'A', 'CYS':'C', 'ASP':'D', 'GLU':'E', 'PHE':'F', 'GLY':'G', 'HIS':'H',
'HID':'H', 'HIE':'H', 'ILE':'I', 'LYS':'K', 'LEU':'L', 'MET':'M', 'ASN':'N',
'PRO':'P', 'GLN':'Q', 'ARG':'R', 'SER':'S', 'THR':'T', 'VAL':'V', 'TRP':'W',
'TYR':'Y'
}

backbone_link_atoms = set({"P","O3'","O3*","N","C"})
    
class bnsTopology:
    
    def __init__(self, pdb_path, debug=False):
        self.pdb_path = pdb_path
        parser = PDBParser(PERMISSIVE=1)
        self.st = parser.get_structure('st', self.pdb_path)
        self.debug = debug
        
    def run(self):
# 1. Detecting Chains        
        bckats = []
        for at in self.st.get_atoms():
            if at.id in backbone_link_atoms:
                bckats.append(at)
        
        nbsearch = NeighborSearch(bckats)        

        sets = []
        for at1, at2 in nbsearch.search_all(COVLNK):
            if at1.get_parent() == at2.get_parent():
                continue
            
            atom1 = Atom(at1)
            atom2 = Atom(at2)
            
            if self.debug:
                print ("#DEBUG: ", atom1.atid(), atom2.atid())

            i = findInSetList(sets,atom1.resid());
            j = findInSetList(sets,atom2.resid())

            if i == -1 and j == -1:
                s = residueset()
                s.add(atom1.resid())
                s.add(atom2.resid())
                sets.append(s)
            
            elif i != -1 and j != -1 and i != j: 
                sets[i].union(sets[j])
                del sets[j]
            
            elif j == -1:
                sets[i].add(atom2.resid())
            
            elif i == -1:
                sets[j].add(atom1.resid())            
            
            if self.debug:
                for s in sorted(sets,key=lambda s: int(s.ini)):
                    print ("#DEBUG:" ,s)
        
        print ("#INFO: Found ",len(sets)," chain(s)")
        
        for s in sorted(sets,key=lambda s: int(s.ini)):
            print (s.ini,":",s.getSequence())
# Nucleotide list
        print ("#INFO: ResidueList")
        for s in sorted(sets,key=lambda s: int(s.ini)):
            print (','.join(s.getResidueList()))
            
def findInSetList(sets,item):
    i = 0
    while i < len(sets) and item not in sets[i].res:
        i = i + 1
    if i == len(sets):
        return -1
    else:
        return i
    
class residueset:
    def __init__(self):
        self.ini=0
        self.res = set()
        
    def add(self,resid):
        self.res.add(resid)
        (rn, ch, n) = resid.split(':')   
        if self.ini == 0:
            self.ini = n
        else:
            self.ini=min(self.ini,n)
    
    def union(self,other):
        self.res = self.res.union(other.res)
        self.ini = min(self.ini,other.ini)
              
    def setini(self):
        for r in self.res:
            (rn, ch, n) = r.split(':')   
            self.ini = min (self.ini,n)
    
    def getSequence(self):
        seq={}
        for r in self.res:
            (rn,ch,n)=r.split(':')
            seq[int(n)]=rn
        ss=''    
        for i in sorted(seq.keys()):
            id =seq[i].rstrip().lstrip()
            if id not in aa1c:
                ss= ss + 'X'
            else:
                ss=ss+aa1c[id]
        return ss
    
    def getResidueList(self):
        seq=[]
        for r in self.res:
            (rn,ch,n)=r.split(':')
            seq.append(str(n)+"-"+aa1c[rn.rstrip().lstrip()])
        return sorted(seq)
            
    
    def getAtoms(self, st):
        atlist=[]
        for r in self.res:
            (rn,ch,n)=r.split(':')
            res = st.get_residue(())
    
    def __str__(self):
        return str(self.ini) + ":" + self.getSequence()
    
class Atom():
    def __init__(self, at):
        self.at = at
        self.res = at.get_parent()
        self.ch = self.res.get_parent()
    
    def resid(self):
        return self.res.get_resname() + ":" + self.ch.id + ":" + str(self.res.id[1])
    
    def atid(self):
        return self.resid() + ":" + self.at.id
    
    def resnum(self):
        (rn,ch,n) = self.resid().split(':')
        return n
    
    def __str__(self):
        return self.atid()

def main():
    parser = argparse.ArgumentParser(prog='bnsTopology', description='Basic topology builder for BNS')
    parser.add_argument('--debug', '-d', action='store_true', help='debug', dest='debug')
    parser.add_argument('pdb_path')
    args = parser.parse_args()
    debug = args.debug
    pdb_path = args.pdb_path
    #print (args)
    if not pdb_path:
        parser.print_help()
        sys.exit(2)
    bnsTopology (args.pdb_path, debug=args.debug).run() 

if __name__ == "__main__":
    main()


