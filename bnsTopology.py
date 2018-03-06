#! /usr/bin/python
#
# Simple parser to extact topology of NA from simulation snapshot
# Chain ids or TER not required
#
__author__ = "gelpi"
__date__ = "$29-ago-2017 16:14:26$"

from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBList import PDBList

import bnsTopLib

import sys
import argparse
import urllib



COVLNK = 2.0
HBLNK  = 3.5
BPTHRESDEF = 1.0


def sameResidue(at1,at2):
    return at1.get_parent() == at2.get_parent()

def getOrderedResiduePair(at1,at2,useChains=False):
# Defining residues as res1 < res2
    res1 = bnsTopLib.StructureWrapper.Residue(at1.get_parent(), useChains)
    res2 = bnsTopLib.StructureWrapper.Residue(at2.get_parent(), useChains)
    if res1 < res2:
        return [res1,res2]
    else:
        return [res2,res1]

def getOrderedAtomPair(at1,at2,useChains=False):
    atom1 = bnsTopLib.StructureWrapper.Atom(at1, useChains)
    atom2 = bnsTopLib.StructureWrapper.Atom(at2, useChains)
    if atom1 < atom2:
        return [atom1,atom2]
    else:
        return [atom2,atom1]

def main():

    argparser = argparse.ArgumentParser(
                prog='bnsTopology',
                description='Basic topology builder for BNS'
            )

    argparser.add_argument(
        '--debug', '-d',
        action='store_true',
        dest='debug',
        help='Produce DEBUG output'
    )

    argparser.add_argument(
        '--usechains',
        action='store_true',
        dest='usechains',
        help='Use PDB file chain ids'
    )

    argparser.add_argument(
        '--graphml',
        dest='graphml',
        help='Produce GraphML output file'
    )

    argparser.add_argument(
        '--json',
        dest='json',
        help='Produce Json output file'
    )

    argparser.add_argument(
        '--bpthres',
        type = float,
        action='store',
        help='BP Score min value ('+str(BPTHRESDEF)+')',
        dest='bpthres',
        default = BPTHRESDEF,
    )
    
    argparser.add_argument(
        '--contacts',
        action='store_true',
        dest='contacts',
        help='Calculate polar contacts between chains'
    )

    argparser.add_argument('pdb_path')

    args = argparser.parse_args()

    debug = args.debug
    pdb_path = args.pdb_path
    useChains = args.usechains
    graphml = args.graphml
    json = args.json
    bpthres = args.bpthres
    contacts = args.contacts

    if not pdb_path:
        argparser.print_help()
        sys.exit(2)

    if "pdb:"in pdb_path:
        pdbl= PDBList(pdb='tmpPDB')
        try:
            pdb_path = pdb_path[4:].upper()
            pdb_path=pdbl.retrieve_pdb_file(pdb_path)
            parser = MMCIFParser()
            useChains=True
            format='cif'
        except IOError:
            print ("#ERROR: fetching PDB "+pdb_path)
            sys.exit(2)
    else:
        try:
            if '.pdb' in pdb_path:
                parser = PDBParser(PERMISSIVE=1)
                format='pdb'
            elif '.cif' in pdb_path:
                parser = MMCIFParser(PERMISSIVE=1)
                format='cif'
            else:
                print ('#ERROR: unknown filetype')
                sys.exit(2)
        except OSError:
            print ("#ERROR: loading PDB")
            sys.exit(2)
    try:
        if '.pdb' in pdb_path:
            parser = PDBParser(PERMISSIVE=1)
            format='pdb'
        elif '.cif' in pdb_path:
            parser = MMCIFParser()
            format='cif'
        else:
            print ('#ERROR: unknown filetype')
            sys.exit(2)
    except OSError:
        print ("#ERROR: parsing PDB")
        sys.exit(2)
    try:
        st = parser.get_structure('st', pdb_path)
    except OSError:
        print ("#ERROR: parsing PDB")
        sys.exit(2)
# Checking for models
    if len(st) > 1:
        print ("#WARNING: Several Models found, using only first")
    
# Using Model 0 any way
    st = st[0]

#Checking for chain ids in input
    if not useChains and len(st) > 1:
        print ("#WARNING: Input PDB contains more than one chain ids ",\
               "when chains ids not used, consider renumbering ")

####### Internal renumbering 

    if useChains:
        i=1
        for r in st.get_residues():
           r.index = i
           i=i+1
    if format == 'cif':
        i=1
        for at in st.get_atoms():
            at.serial_number = i
            i=i+1

# json output
    if json:
        jsondata = bnsTopLib.JSONWriter.JSONWriter()
        JSchainData = {}
        JScovLinks=[]
        JSContacts=[]
        JSbps=[]
        JShelFrags=[]
        

# Graphml output
    if graphml:
        xml = bnsTopLib.graphml.GraphmlWriter()

# 1. Detecting Chains
    bckats = []
    for at in st.get_atoms():
        if '.'+at.id in bnsTopLib.StructureWrapper.backbone_link_atoms:
            bckats.append(at)

    nbsearch = NeighborSearch(bckats)

    covLinkPairs = []

    for at1, at2 in nbsearch.search_all(COVLNK):
        if sameResidue(at1,at2):
            continue
        covLinkPairs.append (getOrderedResiduePair(at1,at2,useChains))

    chList = bnsTopLib.ResidueSet.ResidueSetList(covLinkPairs)
    
    print ("#INFO: Found ",chList.n," chain(s)")
    if json:
        JSchainData['nOfChains']=chList.n
        JSchainData['chains'] = []
        
    for s in chList.getSortedSets():
        print (s)
        if json:
            JSchainData['chains'].append({'iniRes' : str(s.inir), 'finRes' : str(s.finr), 'sequence' : s.getSequence()})
    if json:
        jsondata.insert('chains',JSchainData)

# Nucleotide list
    print ("#INFO: Residue Ids List")
    i=0
    for s in chList.getSortedSets():
        print (s.inir,'-',s.finr,':', ','.join(s.getResidueIdList()))
        if graphml:
            for res in s.items:
                xml.addResidue(i,res)
        i=i+1

    if debug:
        print ("#DEBUG: covalently linked residue pairs")
        for r in sorted(covLinkPairs, key=lambda i: i[0].residue.index):
            print ("#DEBUG: ",r[0].resid(),r[1].resid())
    if graphml:
        for r in sorted(covLinkPairs, key=lambda i: i[0].residue.index):
            xml.addBond('ch',r[0],r[1])
    if json:
        for r in sorted(covLinkPairs, key=lambda i: i[0].residue.index):
            JScovLinks.append([r[0].resid(True),r[1].resid(True)])
        jsondata.insert('covLinks',JScovLinks)
# Contacts
    if contacts:
        print ("#INFO: Getting interchain contacts")
        conts={}
        for ch1 in chList.getSortedChains():
            conts[ch1]={}
            for ch2 in chList.getSortedChains():
                if ch2.ini <= ch1.ini:
                    continue
                conts[ch1][ch2]=[]
                ats1 = []
                s1 = ch1._getResidues()
                s2 = ch2._getResidues()
                for r1 in s1:
                    for r2 in s2:
                        cont = s1[r1].getClosestContact(s2[r2],HBLNK)
                        if cont:
                            [at1,at2,d] = cont
                            [atom1,atom2] = getOrderedAtomPair(at1,at2,useChains)
                            print (atom1.atid(True), atom2.atid(True),d)
                            conts[ch1][ch2].append([atom1,atom2,d])
                            if json:
                                JSContacts.append({'ats':[atom1.atid(True), atom2.atid(True)],'distance':float(d)})
                
        if json:
            jsondata.insert('contacts',JSContacts)

#Base Pairs
    print ("#INFO: Base pairs found")

    hbAtoms = set()
    for hb in bnsTopLib.StructureWrapper.hbonds['WC']:
        for rat in hb:
            hbAtoms.add(rat)
    wcats = []
    for at in st.get_atoms():
        if bnsTopLib.StructureWrapper.Atom(at).attype() in hbAtoms:
            wcats.append(at)
    wc_nbsearch = NeighborSearch(wcats)

    nhbs={}
    for at1, at2 in wc_nbsearch.search_all(HBLNK):
        if sameResidue(at1,at2):
            continue

        [atom1,atom2] = getOrderedAtomPair(at1,at2,useChains)

        res1 = bnsTopLib.StructureWrapper.Residue(atom1.at.get_parent(),useChains)
        res2 = bnsTopLib.StructureWrapper.Residue(atom2.at.get_parent(),useChains)

        if [atom1.attype(), atom2.attype()] not in bnsTopLib.StructureWrapper.hbonds['WC'] and \
           [atom2.attype(), atom1.attype()] not in bnsTopLib.StructureWrapper.hbonds['WC']:
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
            bps.append(bnsTopLib.StructureWrapper.BPair(r1,pair,maxv))

    bpsref={}

    for bp in sorted(bps):
        print (bp, '(',bp.type,'):',','.join(bp.compsIds()), '(',str(bp.score),')')
        bpsref[bp.r1.residue.index]=bp
        if graphml:
            xml.addBond("bp",bp.r1, bp.r2)
        if json:
            JSbps.append({'id': bp.bpid(),'type':bp.type,'score':float(bp.score), 'comps':bp.compsIds()})
    if json:
        jsondata.insert('bpList',JSbps)

# Bpair steps from neighbour bps, relays on residue renumbering
    print ("#INFO: Base Pair steps")

    bpsteps=[]
    for bp1 in sorted(bps):        
        for bp2 in sorted(bps):
            if bp1 < bp2:
                if bp1.r1.residue.index == bp2.r1.residue.index-1 and \
                   bp1.r2.residue.index == bp2.r2.residue.index+1:
                    bpsteps.append(bnsTopLib.StructureWrapper.BPStep(bp1,bp2))

    bpstpref={}
    if json:
        JSbpsteps=[]
    for bpstp in sorted(bpsteps):
        print (bpstp, '(',bpstp.type,'):', ','.join(bpstp.compsIds()))
        bpstpref[bpstp.bp1.r1.residue.index]=bpstp
        if json:
            JSbpsteps.append({'id':bpstp.stepid(),'type':bpstp.type,'comps':bpstp.compsIds()})
    if json:
        jsondata.insert('bpStepsList',JSbpsteps)
# Continuous helical segments from stretches of overlapping bsteps
    print ("#INFO: Helical segments")

    bpstepPairs=[]
    for bpst1 in sorted(bpsteps):
        for bpst2 in sorted(bpsteps):
            if bpst1 < bpst2:
                if bpst1.bp1.r1.residue.index == bpst2.bp1.r1.residue.index-1:
                    bpstepPairs.append([bpst1,bpst2])

    fragList=bnsTopLib.ResidueSet.BPSSetList(bpstepPairs)
   
    for fr in fragList.getSortedSets():
        frag={}
        for i in range (fr.ini,fr.fin+1):
            for bb in bpstpref[i].comps():
                frag[bb]=1
        seq=[]
        for bp in sorted(frag.keys(), key=bnsTopLib.StructureWrapper.BPair.__index__):
           seq.append(bp.bpid())
        print (','.join(seq))
        if json:
            JShelFrags.append(seq)
    if json:
        jsondata.insert('HelicalFrags',JShelFrags)
    if graphml:
        xml.save(graphml)
    if json:
        jsondata.save(json)
        
if __name__ == "__main__":
    main()
