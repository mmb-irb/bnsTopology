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

COVLNK = 2.0
HBLNK  = 3.5
verbose = True

aa1c = {
'DA' :'A', 'DC' :'C', 'DG' :'G', 'DT' :'T',
 'A' :'A',  'C' :'C',  'G' :'G',  'U' :'U',
'DA3':'A', 'DC3':'C', 'DG3':'G', 'DT3':'T',
 'A3':'A',  'C3':'C',  'G3':'G',  'U3':'U',
'DA5':'A', 'DC5':'C', 'DG5':'G', 'DT5':'T',
 'A5':'A',  'C5':'C',  'G5':'G',  'U5':'U',
'ALA':'A', 'CYS':'C', 'ASP':'D', 'GLU':'E', 'PHE':'F', 'GLY':'G', 'HIS':'H',
'HID':'H', 'HIE':'H', 'ILE':'I', 'LYS':'K', 'LEU':'L', 'MET':'M', 'ASN':'N',
'PRO':'P', 'GLN':'Q', 'ARG':'R', 'SER':'S', 'THR':'T', 'VAL':'V', 'TRP':'W',
'TYR':'Y'
}
    
class bnsTopology:
    
    def __init__(self, pdb_path):
        self.pdb_path = pdb_path
        
    def run(self):
        bckats = []
        sets = []
        parser = PDBParser(PERMISSIVE=1)
        st = parser.get_structure('struc', self.pdb_path)
# Detecting Chains        
        for at in st.get_atoms():
            if at.id in set({"P","O3'","N","C"}):
                bckats.append(at)
        nbsearch = NeighborSearch(bckats)        
        for at1, at2 in nbsearch.search_all(COVLNK):
            if at1.get_parent() == at2.get_parent():
                continue
            atom1 = atom(at1)
            atom2 = atom(at2)
            print atom1, atom2
            i = j = 0

            while i < len(sets) and atom1.resid() not in sets[i].res:
                i = i + 1
            
            while j < len(sets) and atom2.resid() not in sets[j].res:
                j = j + 1
            
            if i == len(sets) and j == len(sets):
                sets.append(residueset())
                sets[len(sets)-1].add(atom1.resid())
                sets[len(sets)-1].add(atom2.resid())
            
            elif i < len(sets) and j < len(sets) and i != j: 
                sets[i].union(sets[j])
                del sets[j]
            
            elif i < len(sets):
                sets[i].add(atom2.resid())
            
            elif j < len(sets):
                sets[j].add(atom1.resid())            
            
            if verbose:
                for s in sorted(sets,key=lambda s: int(s.ini)):
                    print s
        
        print "Found ",len(sets)," chain(s)"
        
        for s in sorted(sets,key=lambda s: int(s.ini)):
            print s.ini,":",s.getSequence()
            
class residueset:
    def __init__(self,usechains=False):
        self.ini=0
        self.res = set()
        self.usechains = usechains
        
    def add(self,resid):
        self.res.add(resid)
        if self.usechains:
            (rn, ch, n) = resid.split(':')   
        else:
            (rn, n) = resid.split(':')   
        if self.ini == 0:
            self.ini = n
        else:
            self.ini=min(self.ini,n)
    
    def union(self,other):
        self.res = self.res.union(other.res)
        self.ini = min(self.ini,other.ini)
           
    
    def setini(self):
        for r in self.res:
            if self.usechains:
                (rn, ch, n) = r.split(':')   
            else:
                (rn, n) = r.split(':')   
            self.ini = min (self.ini,n)
    
    def getSequence(self):
        seq={}
        for r in self.res:
            if self.usechains:
               (rn,ch,n)=r.split(':')
            else:
               (rn,n)=r.split(':')
            seq[int(n)]=rn
        ss=''    
        for i in sorted(seq.keys()):
            id =seq[i].rstrip().lstrip()
            if id not in aa1c:
                ss= ss + 'X'
                print "Warning: Unknown Residue ",id
            else:
                ss=ss+aa1c[id]
        return ss
    
    def __str__(self):
        return str(self.ini) + ":" + self.getSequence()
    
class atom:
    def __init__(self, at, usechains=False):
        self.at = at
        self.res = at.get_parent()
        self.ch = self.res.get_parent()
        self.usechains=usechains
    
    def resid(self):
        if self.usechains:
            return self.res.get_resname() + ":" + self.ch.id + ":" + str(self.res.id[1])
        else:
            return self.res.get_resname() + ":" + str(self.res.id[1])
    
    def atid(self):
        return self.resid() + ":" + self.at.id
    
    def resnum(self):
        if self.usechains:
            (rn,ch,n) = self.resid().split(':')
        else:
            (rn,n) = self.resid().split(':')
        return n
    
    def __str__(self):
        return self.atid()
    


if __name__ == "__main__":
    bnsTopology (sys.argv[1]).run()

