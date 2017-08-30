#! /usr/bin/python

# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

__author__ = "gelpi"
__date__ = "$29-ago-2017 16:14:26$"

import sys
import os
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.NeighborSearch import NeighborSearch

COVLNK = 2.0
aa1c={
'DA':'A',
'DC':'C',
'DG':'G',
'DT':'T',
'A':'A',
'C':'C',
'G':'G',
'U':'U',
'DA3':'A',
'DC3':'C',
'DG3':'G',
'DT3':'T',
'A3':'A',
'C3':'C',
'G3':'G',
'U3':'U',
'DA5':'A',
'DC5':'C',
'DG5':'G',
'DT5':'T',
'A5':'A',
'C5':'C',
'G5':'G',
'U5':'U',
}
    
class naTopology (object):
    
    
    def __init__(self, pdb_path):
        self.pdb_path = pdb_path
        
    def run(self):
        bckats=[]
        sets=[]
        parser = PDBParser(PERMISSIVE=1)
        st = parser.get_structure('struc', self.pdb_path)
        for at in st.get_atoms():
            if (at.id == 'P') or (at.id == "O3'"):
                bckats.append(at)
        nbsearch = NeighborSearch(bckats)        
        for at1, at2 in nbsearch.search_all(COVLNK):
            if at1.get_parent() == at2.get_parent():
                continue
            print atid(at1),atid(at2)
            i=j=0
            while i < len(sets) and resid(at1) not in sets[i]:
                i = i +  1
            while j < len(sets) and resid(at2) not in sets[j]:
                j = j + 1
            if i == len(sets) and j == len(sets):
                sets.append(set())
                sets[len(sets)-1].add(resid(at1))
                sets[len(sets)-1].add(resid(at2))
            elif i < len(sets) and j < len(sets) and i != j: 
                sets[i]=sets[i].union(sets[j])
                del sets[j]
            elif i < len(sets):
                sets[i].add(resid(at2))
            elif j < len(sets):
                sets[j].add(resid(at1))
        print len(sets)
#        print (sets)
        for s in sets:
            seq={}
            for r in s:
                (rn,ch,n)=r.split(':')
                seq[int(n)]=rn
            ss=''    
            print seq
            for i in sorted(seq.keys()):
                 ss=ss+aa1c[seq[i].rstrip().lstrip()]
            print ss
            

def atid (at):
    return resid(at)+":"+at.id
    
def resid(at):
    res = at.get_parent()
    ch = res.get_parent()
    return res.get_resname()+":"+ch.id+":"+str(res.id[1])
    
        
        

if __name__ == "__main__":
    naTopology (sys.argv[1]).run()

