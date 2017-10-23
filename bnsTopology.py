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

backbone_link_atoms = set({".P",".O3'",".O3*",".N",".C"})
hbonds = {
    'WC': [
            ["G.N1", "C.N3"],
            ["G.N2", "C.O2"],
            ["G.O6", "C.N4"],
            ["A.N6", "T.O4"],
            ["A.N1", "T.N3"],
            ["A.N6", "U.O4"],
            ["A.N1", "U.N3"]
    ],
    'HG': []
}

class ChainList():
    def __init__(self):
        self.chains=[]
        self.n=0
   
    def find(self,item):
        i = 0
        while i < self.n and item not in self.chains[i].residues:
            i = i + 1
        if i == self.n:
            return -1
        else:
            return i        
    def append(self,item):
        self.chains.append(item)
        self.n = len(self.chains)

    def delete(self,i):
        del self.chains[i]
        self.n = len(self.chains)

    def getSortedChains(self):
        return sorted(self.chains,key=lambda s: int(s.ini))
        

class Chain():
    def __init__(self):
        self.ini=0
        self.fin=0
        self.residues = set()
        
    def add(self,r):
        self.residues.add(r)
        if self.ini == 0:
            self.ini = int(r.resNum())
        else:
            self.ini = min(self.ini,int(r.resNum()))
        if self.fin == 0:
            self.fin = int(r.resNum())
        else:
            self.fin = max(self.fin,int(r.resNum()))
    
    def union(self,other):
        self.residues = self.residues.union(other.residues)
        self.ini = min(self.ini,other.ini)
        self.fin = max(self.fin,other.fin)

    def getSequence(self):
        seq=self._getResidues()
        ss=''    
        for i in sorted(seq.keys()):
            ss=ss+seq[i]._getOneLetterResidueCode()
        return ss
   
    def getResidueIdList(self):
        seq=self._getResidues()
        seql = []
        for i in sorted(seq.keys()):
            seql.append(str(i)+"-"+seq[i]._getOneLetterResidueCode())
        return seql
    
    def _getResidues(self):
        seq={}
        for r in self.residues:
            seq[r.resNum()] = r
        #print (seq)
        return seq
    def __str__(self):
        return str(self.ini) + "-" + str(self.fin)+ ":" + self.getSequence()

class Residue():
    oneLetterResidueCode = {
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
    def __init__(self, r, useChains=False):
        self.residue = r
        self.useChains = useChains
    
    def resid(self, compact=False):
        if self.useChains:
            ch = ":"+self.residue.get_parent().id
        else:
            ch = ''
        if compact:
            return self._getOneLetterResidueCode()+ch+ str(self.residue.id[1])
        else:
            return self.residue.get_resname() + ch + ':' + str(self.residue.id[1])       

    def bnsid(self):
        return str(self.resNum())+"-"+self._getOneLetterResidueCode()
    
    def resNum(self):
        return self.residue.id[1]

    def _getOneLetterResidueCode(self):
        id = self.residue.get_resname().rstrip().lstrip()
        if not id in Residue.oneLetterResidueCode:
            return 'X'
        else:
            return Residue.oneLetterResidueCode[id]
    
    def __hash__(self):
        return hash(self.resid())
    
    def __eq__(self, other):        
        if self.useChains:
            id = self.residue.get_parent().id + str(self.residue.id[1])
            otherid = other.residue.get_parent().id + str(other.residue.id[1])
        else:
            id = self.residue.id[1]
            otherid = other.residue.id[1]
        return id == otherid

    def __lt__(self, other):        
        if self.useChains:
            id = self.residue.get_parent().id + str(self.residue.id[1])
            otherid = other.residue.get_parent().id + str(other.residue.id[1])
        else:
            id = self.residue.id[1]
            otherid = other.residue.id[1]
        return id < otherid
    def __str__(self):
        return self.resid()
            
class Atom():
    def __init__ (self,at,useChains=False):
         self.at=at
         self.useChains=useChains
    
    def atid(self, compact=False):
        return self.resid(compact)+"."+self.at.id
    
    def resid(self, compact=False):
        return Residue(self.at.get_parent(),self.useChains).resid(compact)
  
    def attype(self):
        return Residue(self.at.get_parent(),self.useChains)._getOneLetterResidueCode()+'.'+self.at.id
    
    def _hbscore(self,other):    
        d = self.at-other.at
        return 2.6875-0.625*d
    
    def __str__(self):
        return self.atid()

class BPair():
    def __init__(self,r1,r2,score):
        self.r1=r1
        self.r2=r2
        types =[self.r1._getOneLetterResidueCode(), self.r2._getOneLetterResidueCode()]
        self.type=''.join(sorted(types))
        self.score=score
    
    def bpid(self):
        return str(self.r1.resNum())+"-"+self.r1._getOneLetterResidueCode()+self.r2._getOneLetterResidueCode()
    
    def comps(self):
        return [self.r1.bnsid(),self.r2.bnsid()]
    
    def __eq__(self,other):
        return self.r1==other.r1 and self.r2 == other.r2
    
    def __lt__(self,other):
        return self.r1<other.r1
    
    def __str__(self):
        return self.bpid()

class BPStep():
    def __init__(self,bp1,bp2):
        self.bp1 =bp1
        self.bp2 =bp2
        bps = [self.bp1.r1._getOneLetterResidueCode() + self.bp2.r1._getOneLetterResidueCode(), self.bp2.r2._getOneLetterResidueCode() + self.bp1.r2._getOneLetterResidueCode()]
        self.type= ''.join(sorted(bps))
        
    def stepid(self):
        return str(self.bp1.r1.resNum())+"-"+self.bp1.r1._getOneLetterResidueCode() + self.bp2.r1._getOneLetterResidueCode() + self.bp2.r2._getOneLetterResidueCode() + self.bp1.r2._getOneLetterResidueCode()
    
    def comps(self):
        return [self.bp1.bpid(),self.bp2.bpid()]
    
    def resNum(self):
        return self.bp1.r1.resNum()
    
    def __eq__(self,other):
        return self.bp1==other.bp1 and self.bp2 == other.bp2
    
    def __lt__(self,other):
        return self.bp1<other.bp1
    
    def __str__(self):
      return self.stepid()  
    
    def __hash__(self):
        return hash(self.stepid())
        
def main():
    
    parser = argparse.ArgumentParser(prog='bnsTopology', description='Basic topology builder for BNS')
    parser.add_argument('--debug', '-d', action='store_true', help='debug', dest='debug')
    parser.add_argument('--usechains', action='store_true', help='Use Chain ids', dest='usechains')
    parser.add_argument('pdb_path')
    
    args = parser.parse_args()
    
    debug = args.debug
    pdb_path = args.pdb_path
    useChains = args.usechains
    
    if not pdb_path:
        parser.print_help()
        sys.exit(2)        

    parser = PDBParser(PERMISSIVE=1)
    
    try:
        st = parser.get_structure('st', pdb_path)
    except OSError:
        print ("#ERROR: loading PDB")
        sys.exit(2)

# Checking for models
    if len(st) > 1:
        print ("#WARNING: Several Models found, using only first")
# Using Model 0 any way
    st = st[0]   
#Checking for chain ids in input
    if not useChains and len(st) > 1:
        print ("#WARNING: Input PDB contains more than one chain ids when chains ids not used, consider renumbering or ")

# 1. Detecting Chains        
    bckats = []
    for at in st.get_atoms():
        if '.'+at.id in backbone_link_atoms:
            bckats.append(at)
        
    nbsearch = NeighborSearch(bckats)        

    covLinkPairs = set()    
    chList = ChainList()
    for at1, at2 in nbsearch.search_all(COVLNK):
        if at1.get_parent() == at2.get_parent():
            continue
# Defining residues as res1 < res2
        if Residue(at1.get_parent(),useChains) < Residue(at2.get_parent(),useChains):
            res1 = Residue(at1.get_parent(), useChains)
            res2 = Residue(at2.get_parent(), useChains)
        else:
            res2 = Residue(at1.get_parent(), useChains)
            res1 = Residue(at2.get_parent(), useChains)
        covLinkPairs.add ((res1,res2))
        i = chList.find(res1)
        j = chList.find(res2)
        if i == -1 and j == -1:
            s = Chain()
            s.add(res1)
            s.add(res2)
            chList.append(s)
        elif i != -1 and j != -1 and i != j: 
            chList.chains[i].union(chList.chains[j])
            chList.delete(j)
        elif j == -1:
            chList.chains[i].add(res2)
        elif i == -1:
            chList.chains[j].add(res1)            
     
        if debug:
            for s in chList.getSortedChains():
                print ("#DEBUG:" ,s)
    print ("#INFO: Found ",chList.n," chain(s)")
    for s in chList.getSortedChains():
        print (s)
# Nucleotide list
    print ("#INFO: Residue Ids List")
    for s in chList.getSortedChains():
        print (s.ini,'-',s.fin,':', ','.join(s.getResidueIdList()))
 
    if debug:
        print ("#DEBUG: covalently linked residue pairs")
        for r in sorted(covLinkPairs, key=lambda i: i[0].resNum()):
            print ("#DEBUG: ",r[0].resid(),r[1].resid())
#Base Pairs
    print ("#INFO: Base pairs found")
      
    hbAtoms = set()
    for hb in hbonds['WC']:
        for rat in hb:
            hbAtoms.add(rat)
    wcats = []    
    for at in st.get_atoms():
        if Atom(at).attype() in hbAtoms:
            wcats.append(at)
    wc_nbsearch = NeighborSearch(wcats)        
    wcsets = []
     
    nhbs={}
    for at1, at2 in wc_nbsearch.search_all(HBLNK):
        if at1.get_parent() == at2.get_parent():
            continue    
# Defining atoms being res1 < res2
        if Residue(at1.get_parent(),useChains) < Residue(at2.get_parent(),useChains):
            atom1 = Atom(at1, useChains)
            atom2 = Atom(at2, useChains)
        else:
            atom1 = Atom(at2, useChains)
            atom2 = Atom(at1, useChains)
        res1 = Residue(atom1.at.get_parent(),useChains)
        res2 = Residue(atom2.at.get_parent(),useChains)
        if [atom1.attype(), atom2.attype()] not in hbonds['WC'] and [atom2.attype(), atom1.attype()] not in hbonds['WC']:
            continue
        if (res1,res2) in covLinkPairs:    
            continue
        if res1 not in nhbs:
            nhbs[res1] = {}
        if res2 not in nhbs[res1]:
            nhbs[res1][res2]=0
        nhbs[res1][res2]=nhbs[res1][res2]+atom1._hbscore(atom2)
        if debug:
            print ("#DEBUG: ", res1,res2,atom1._hbscore(atom2))
          
    if debug:
        print ("#DEBUG: HB count per pair of residues")
        for r1 in nhbs.keys():
            for r2 in nhbs[r1].keys():
                print ("#DEBUG: ", r1, r2, nhbs[r1][r2])

    bps = []
    for r1 in nhbs.keys():
        maxv=0.
        pair=''
        for r2 in nhbs[r1].keys():            
            if nhbs[r1][r2]>maxv:
                pair=r2
                maxv=nhbs[r1][r2]                
        bps.append(BPair(r1,pair,maxv))
    
    bpsref={}
    
    for bp in sorted(bps):
        print (bp, '(',bp.type,'):',','.join(bp.comps()), '(',str(bp.score),')')
        bpsref[bp.r1.resNum()]=bp
            
    print ("#INFO: Base Pair steps")
    
    bpsteps=[]
    for bp1 in sorted(bps):
        for bp2 in sorted(bps):
            if bp1 < bp2:
                if bp1.r1.resNum() == bp2.r1.resNum()-1 and bp1.r2.resNum() == bp2.r2.resNum()+1:
                    bpsteps.append(BPStep(bp1,bp2))
    
    bpstpref={}
    for bpstp in sorted(bpsteps):
        print (bpstp, '(',bpstp.type,'):', ','.join(bpstp.comps()))
        bpstpref[bpstp.resNum()]=bpstp
  
    print ("#INFO: Helical segments")
    
    fragList=ChainList()
    
    for bpst1 in sorted(bpsteps):
        for bpst2 in sorted(bpsteps):    
            if bpst1 < bpst2:
                if bpst1.bp1.r1.resNum() == bpst2.bp1.r1.resNum()-1:
                    i = fragList.find(bpst1)
                    j = fragList.find(bpst2)
                    if i == -1 and j == -1:
                        s = Chain()
                        s.add(bpst2)
                        s.add(bpst1)
                        fragList.append(s)
                    elif i != -1 and j != -1 and i != j: 
                        fragList.chains[i].union(fragList.chains[j])
                        fragList.delete(j)
                    elif j == -1:
                        fragList.chains[i].add(bpst2)
                    elif i == -1:
                        fragList.chains[j].add(bpst1)            
                          
    for fr in fragList.getSortedChains():
        frag={}
        for i in range (fr.ini,fr.fin+1):
            for bb in bpstpref[i].comps():
                frag[bb]=1
        print (",".join(sorted(frag.keys())))
    
if __name__ == "__main__":
    main()


