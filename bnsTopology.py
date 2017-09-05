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
import inspect

COVLNK = 2.0
HBLNK  = 3.5
#verbose = True
#debug = True



backbone_link_atoms = set({":P",":O3'",":O3*",":N",":C"})
hbonds = {
    'WC': [
            ["G:N1", "C:N3"],
            ["G:N2", "C:O2"],
            ["G:O6", "C:N4"],
            ["A:N6", "T:O4"],
            ["A:N1", "T:N3"],
            ["A:N6", "U:O4"],
            ["A:N1", "U:N3"]
    ],
    'HG': []
}


    
class bnsTopology:
    
    def __init__(self, pdb_path, debug=False, useChains=False):
        self.pdb_path = pdb_path
        parser = PDBParser(PERMISSIVE=1)
        try:
            self.st = parser.get_structure('st', self.pdb_path)
        except OSError:
            print ("#ERROR: loading PDB")
            sys.exit(2)
        self.debug = debug
        if useChains == 'auto':
            self.useChains=''
            for at in self.st.get_atoms():
                useChains = at.get_parent().get_parent().id or useChains
            self.useChains
        else:
            self.useChains=useChains
#        print (self.st)
        
    def run(self):
# 1. Detecting Chains        
        bckats = []
        for at in self.st.get_atoms():
            if ':'+at.id in backbone_link_atoms:
                bckats.append(at)
        
        nbsearch = NeighborSearch(bckats)        

        
        self.covLinkPairs = set()
        
        self.chList = ChainList()

        for at1, at2 in nbsearch.search_all(COVLNK):
            if at1.get_parent() == at2.get_parent():
                continue
#            
            if at1.get_parent().id[1] < at2.get_parent().id[1]:
                atom1 = Atom(at1,useChains=self.useChains)
                atom2 = Atom(at2,useChains=self.useChains)
            else:
                atom1 = Atom(at2,useChains=self.useChains)
                atom2 = Atom(at1,useChains=self.useChains)

            self.covLinkPairs.add (atom1.resid(compact=False)+','+atom2.resid(compact=False))
            
            if self.debug:
                print ("#DEBUG: ", atom1.atid(), atom2.atid())

            i = self.chList.find(atom1);
            j = self.chList.find(atom2)

            if i == -1 and j == -1:
                s = Chain()
                s.add(atom1)
                s.add(atom2)
                self.chList.append(s)
            
            elif i != -1 and j != -1 and i != j: 
                chList.chains[i].union(sets[j])
                del self.chList.chains[j]
            
            elif j == -1:
                self.chList.chains[i].add(atom2.resid())
            
            elif i == -1:
                self.chList.chains[j].add(atom1.resid())            
            
            if self.debug:
                for s in self.chList.getSortedChains():
                    print ("#DEBUG:" ,s)
        
        print ("#INFO: Found ",self.chList.n," chain(s)")
        
        for s in self.chains.getSortedChains():
            print (s.ini,":",s.getSequence())
# Nucleotide list
        print ("#INFO: ResidueList")
        for s in sorted(self.chains,key=lambda s: int(s.ini)):
            print (','.join(s.getResidueIdList()))
 
        if self.debug:
            print ("#DEBUG: linked residue pairs")
            print ("#DEBUG: ",self.covLinkPairs)
#Base Pairs
#        print ("#INFO: Base pairs found")
#        
#        hbAtoms = set()
#        for hb in hbonds['WC']:
#            for rat in hb:
#                hbAtoms.add(rat)
#
#        wcats = []
#        for at in self.st.get_atoms():
#            if getoneletter(at.get_parent().get_resname())+':'+at.id in hbAtoms:
#                wcats.append(at)
#
#        wc_nbsearch = NeighborSearch(wcats)        
#
#        wcsets = []
#        
#        hbs0=[]
#        nhbs={}
#        for at1, at2 in wc_nbsearch.search_all(HBLNK):
#
#            if at1.get_parent() == at2.get_parent():
#                continue    
#
#            if at1.get_parent().id[1] < at2.get_parent().id[1]:
#                atom1 = Atom(at1)
#                atom2 = Atom(at2)
#            else:
#                atom1 = Atom(at2)
#                atom2 = Atom(at1)
#
#            if [atom1.attype(), atom2.attype()] not in hbonds['WC'] and [atom2.attype(), atom1.attype()] not in hbonds['WC']:
#                continue
#            if atom1.resid()+','+atom2.resid() in covlinkres or atom2.resid()+','+atom1.resid() in covlinkres:
#                continue
#            hbs0.append([atom1.atid(),atom2.atid(), atom1.at-atom2.at])
#
#            if atom1.resid() not in nhbs:
#                nhbs[atom1.resid()] = {}
#            if atom2.resid() not in nhbs[atom1.resid()]:
#                nhbs[atom1.resid()][atom2.resid()]=0
#
#            nhbs[atom1.resid()][atom2.resid()]=nhbs[atom1.resid()][atom2.resid()]+_hbscore(atom1.at,atom2.at)
#            
#
#        if self.debug:
#            print ("#DEBUG: candidate WC HBs sorted by distance") 
#            for h in sorted(hbs0,key=lambda d: d[2]):
#                print ("#DEBUGHBS: ",h)
#        if self.debug:
#            print ("#DEBUG: HB count per pair of residues")
#            for r1 in nhbs.keys():
#                for r2 in nhbs[r1].keys():
#                    if r1 == r2:
#                        continue
#                    print ("#DEBUG: ", r1, r2, nhbs[r1][r2])
#                
#        bps={}
#        for r1 in nhbs.keys():
#            maxv=0.
#            pair=''
#            for r2 in nhbs[r1].keys():
#                if nhbs[r1][r2]>maxv:
#                    pair=r2
#                    maxv=nhbs[r1][r2]
#            bps[r1] = r2
#            
#        for bp in bps.keys():
#            print (_bpid([bp,bps[bp]], comps=True))
#            
#        print ("#INFO: Base Pair steps")
#        
#            
#                
#        
#        
#def getoneletter(id):        
#    id = id.rstrip().lstrip()
#    if not id in aa1c:
#        return 'X'
#    else:
#        return aa1c[id]
#    
#def _hbscore(at1,at2):    
#    d = at1-at2
#    return 2.6875-0.625*d
#    
#def _bpid(bp, comps=False):
#    (rn1,ch1,n1)= bp[0].split(':')
#    (rn2,ch2,n2)= bp[1].split(':')
#    rn1 = getoneletter(rn1)
#    rn2 = getoneletter(rn2)
#    bpid = str(n1)+"-"+rn1+rn2
#    if comps:
#        return bpid+"|"+str(n1)+'-'+rn1+','+str(n2)+'-'+rn2
#    else:
#        return bpid
#    
#    
#    

class Atom():
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

    def __init__ (self,at, useChains=False):
         self.at=at
         self.res = at.get_parent()
         self.ch = self.res.get_parent()
         self.mod = self.ch.get_parent()
         self.useChains=useChains

    def resid(self, compact=True):
        if self.useChains:
            ch = self.ch.id+":"
        else:
            ch = ''
        if compact:
            return self._getOneletterResidueCode()+ch+ str(self.res.id[1])
        else:
            return self.res.get_resname() + ":" + ch + str(self.res.id[1])
    
    def atid(self, compact=True):
        return self.resid(compact)+"."+self.at.id

    def _getOneletterResidueCode(self):
        id = self.res.get_resname().rstrip().lstrip()
        if not id in Atom.oneLetterResidueCode:
            return 'X'
        else:
            return Atom.oneLetterResidueCode[id]
#    
#    def attype(self):
#        return getoneletter(self.at.get_parent().get_resname())+':'+self.at.id
#    
    def resNum(self):
        if self.useChains:
            return self.resid(False).split(':')[2]   
        else:
            return self.resid(False).split(':')[1]

    def resName(self):
        return self.resid(False).split(':')[0]   
        
    def __str__(self):
        return self.atid()

class ChainList:
    def __init__(self):
        self.chains=[]
        self.n=0
    
    def find(self,item):
        i = 0
        while i < self.n and item.resid(False) not in self.chains[i].residues:
            i = i + 1
        if i == self.n:
            return -1
        else:
            return i
        
    def append(self,item):
        self.chains.append(item)
        self.n = len(self.chains)

    def getSortedChains(self):
        return sorted(self.chains,key=lambda s: int(s.ini))
        
class Chain:
    def __init__(self):
        self.ini=0
        self.residues = set()
        
    def add(self,atom):
        self.residues.add(atom.resid())
        if self.ini == 0:
            self.ini = int(atom.resNum())
        else:
            self.ini = min(self.ini,int(atom.resNum()))
#    
#    def union(self,other):
#        self.res = self.res.union(other.res)
#        self.ini = min(self.ini,other.ini)
#              
#    def setini(self):
#        for r in self.res:
#            (rn, ch, n) = r.split(':')   
#            self.ini = min (self.ini,int(n))

    def getSequence(self):
        seq=self._getResidues()
##        for r in self.res:
##            (rn,ch,n)=r.split(':')
##            seq[int(n)]=rn
        ss=''    
        for i in sorted(seq.keys()):
            ss=ss+getoneletter(seq[i])
        return ss
    
    def getResidueIdList(self):
        seq=self._getResidues()
        seql = []
        for i in sorted(seq.keys()):
            seql.append(str(i)+"-"+getoneletter(seq[i]))
        return seql
    
    def _getResidues(self):
        seq={}
        for r in self.residues:
            seq[r.resNum()] = r.resName()
        print (seq)
        return seq
#    
##    def getAtoms(self, st):
##        atlist=[]
##        for r in self.res:
##            (rn,ch,n)=r.split(':')
##            res = st.get_residue(())
#    
    def __str__(self):
        return str(self.ini) + ":" + self.getSequence()


def main():
    parser = argparse.ArgumentParser(prog='bnsTopology', description='Basic topology builder for BNS')
    parser.add_argument('--debug', '-d', action='store_true', help='debug', dest='debug')
    parser.add_argument('--usechains', action='store_true', help='Use Chain ids', dest='usechains')
    parser.add_argument('pdb_path')
    args = parser.parse_args()
    debug = args.debug
    pdb_path = args.pdb_path
    #print (args)
    if not pdb_path:
        parser.print_help()
        sys.exit(2)        
    bnsTopology (
        args.pdb_path, 
        debug=args.debug, 
        useChains=args.usechains).run() 

if __name__ == "__main__":
    main()


