#! /usr/bin/python
#
# Simple parser to extact topology of NA from simulation snapshot
# Chain ids or TER not required 
#
__author__ = "gelpi"
__date__ = "$29-ago-2017 16:14:26$"

from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.PDBParser import PDBParser

import bnsTopLib 

import os
import sys
import argparse


COVLNK = 2.0
HBLNK  = 3.5
BPTHRESDEF = 1.0

def main():
  
    parser = argparse.ArgumentParser(
                prog='bnsTopology', 
                description='Basic topology builder for BNS'
            )

    parser.add_argument(
        '--debug', '-d', 
        action='store_true', 
        dest='debug',
        help='Produce DEBUG output'
    )

    parser.add_argument(
        '--usechains', 
        action='store_true', 
        dest='usechains',
        help='Use PDB file chain ids'
    )

    parser.add_argument(
        '--graphml', 
        action='store_true', 
        dest='graphml',
        help='Produce GraphML output file'
    )

    parser.add_argument(
        '--bpthres', 
        type = float,
        action='store', 
        help='BP Score min value ('+str(BPTHRESDEF)+')', 
        dest='bpthres',
        default = BPTHRESDEF,
        
    )

    parser.add_argument('pdb_path')
    
    args = parser.parse_args()
    
    debug = args.debug
    pdb_path = args.pdb_path
    useChains = args.usechains
    graphml = args.graphml
    bpthres = args.bpthres
    
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
        print ("#WARNING: Input PDB contains more than one chain ids ",\
               "when chains ids not used, consider renumbering or ")
# Graphml output
    if graphml:
        xml = bnsTopLib.graphml.GraphmlWriter()

# 1. Detecting Chains        
    bckats = []
    for at in st.get_atoms():
        if '.'+at.id in bnsTopLib.Structure.backbone_link_atoms:
            bckats.append(at)
        
    nbsearch = NeighborSearch(bckats)        

    covLinkPairs = set()    
    chList = bnsTopLib.Chains.ChainList()
    for at1, at2 in nbsearch.search_all(COVLNK):
        if at1.get_parent() == at2.get_parent():
            continue
# Defining residues as res1 < res2
        if bnsTopLib.Structure.Residue(at1.get_parent(),useChains) < bnsTopLib.Structure.Residue(at2.get_parent(),useChains):
            res1 = bnsTopLib.Structure.Residue(at1.get_parent(), useChains)
            res2 = bnsTopLib.Structure.Residue(at2.get_parent(), useChains)
        else:
            res2 = bnsTopLib.Structure.Residue(at1.get_parent(), useChains)
            res1 = bnsTopLib.Structure.Residue(at2.get_parent(), useChains)

        covLinkPairs.add ((res1,res2))
        
        i = chList.find(res1)
        j = chList.find(res2)
        
        if i == -1 and j == -1:
            s = bnsTopLib.Chains.Chain()
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
    i=0
    for s in chList.getSortedChains():
        print (s.ini,'-',s.fin,':', ','.join(s.getResidueIdList()))
        if graphml:
            for res in s.residues:
                xml.addResidue(i,res)
        i=i+1        
 
    if debug:
        print ("#DEBUG: covalently linked residue pairs")
        for r in sorted(covLinkPairs, key=lambda i: i[0].resNum()):
            print ("#DEBUG: ",r[0].resid(),r[1].resid())
    if graphml:
        for r in sorted(covLinkPairs, key=lambda i: i[0].resNum()):
            xml.addBond('ch',r[0],r[1])
                

#Base Pairs
    print ("#INFO: Base pairs found")
      
    hbAtoms = set()
    for hb in bnsTopLib.Structure.hbonds['WC']:
        for rat in hb:
            hbAtoms.add(rat)
    wcats = []    
    for at in st.get_atoms():
        if bnsTopLib.Structure.Atom(at).attype() in hbAtoms:
            wcats.append(at)
    wc_nbsearch = NeighborSearch(wcats)        
    wcsets = []
     
    nhbs={}
    for at1, at2 in wc_nbsearch.search_all(HBLNK):
        if at1.get_parent() == at2.get_parent():
            continue    
# Defining atoms being res1 < res2
        if bnsTopLib.Structure.Residue(at1.get_parent(),useChains) < bnsTopLib.Structure.Residue(at2.get_parent(),useChains):
            atom1 = bnsTopLib.Structure.Atom(at1, useChains)
            atom2 = bnsTopLib.Structure.Atom(at2, useChains)
        else:
            atom1 = bnsTopLib.Structure.Atom(at2, useChains)
            atom2 = bnsTopLib.Structure.Atom(at1, useChains)
        res1 = bnsTopLib.Structure.Residue(atom1.at.get_parent(),useChains)
        res2 = bnsTopLib.Structure.Residue(atom2.at.get_parent(),useChains)
        if [atom1.attype(), atom2.attype()] not in bnsTopLib.Structure.hbonds['WC'] and \
           [atom2.attype(), atom1.attype()] not in bnsTopLib.Structure.hbonds['WC']:
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
        if maxv > bpthres:
            bps.append(bnsTopLib.Structure.BPair(r1,pair,maxv))
    
    bpsref={}
    
    for bp in sorted(bps):
        print (bp, '(',bp.type,'):',','.join(bp.comps()), '(',str(bp.score),')')
        bpsref[bp.r1.resNum()]=bp
        if graphml:
            xml.addBond("bp",bp.r1, bp.r2)

# Bpair steps from neighbour bps, relays on residue renumbering
    print ("#INFO: Base Pair steps")
    
    bpsteps=[]
    for bp1 in sorted(bps):
        for bp2 in sorted(bps):
            if bp1 < bp2:
                if bp1.r1.resNum() == bp2.r1.resNum()-1 and \
                   bp1.r2.resNum() == bp2.r2.resNum()+1:
                    bpsteps.append(bnsTopLib.Structure.BPStep(bp1,bp2))
    
    bpstpref={}
    for bpstp in sorted(bpsteps):
        print (bpstp, '(',bpstp.type,'):', ','.join(bpstp.comps()))
        bpstpref[bpstp.resNum()]=bpstp
  
# Continuous helical segments from stretches of overlapping bsteps
    print ("#INFO: Helical segments")
    
    fragList=bnsTopLib.Chains.ChainList()
    
    for bpst1 in sorted(bpsteps):
        i = fragList.find(bpst1)
        if i == -1:
            s = bnsTopLib.Chains.Chain()
            s.add(bpst1)
            fragList.append(s)
        for bpst2 in sorted(bpsteps):    
            if bpst1 < bpst2:
                if bpst1.bp1.r1.resNum() == bpst2.bp1.r1.resNum()-1:
                    i = fragList.find(bpst1)
                    j = fragList.find(bpst2)
                    if i == -1 and j == -1:
                        s = bnsTopLib.Chains.Chain()
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
    
    if graphml:
        xml.save("output.graphml")
        
if __name__ == "__main__":
    main()
