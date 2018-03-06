#! /usr/bin/python
#
# Simple parser to extact topology of NA from simulation snapshot
# Chain ids or TER not required
#
__author__ = "gelpi"
__date__ = "$29-ago-2017 16:14:26$"

from Bio.PDB.NeighborSearch import NeighborSearch


import bnsTopLib

import sys

COVLNK = 2.0
HBLNK  = 3.5
INTDIST = 6.0
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

    args = bnsTopLib.cmdLine({'BPTHRESDEF': BPTHRESDEF}).parse_args()
# json output
    if args.json:
        jsondata = bnsTopLib.JSONWriter()
        jsondata.data['intdist']= INTDIST
# Graphml output
    if args.graphml:
        xml = bnsTopLib.GraphmlWriter()
# ============================================================================
    loader = bnsTopLib.PDBLoader(args.pdb_path)
    args.useChains= args.useChains or loader.useChains      
    if args.json:
        jsondata.data['useChains']=args.useChains
        jsondata.data['inputFormat']=loader.format
    st = loader.loadStructure()
#Checking for chain ids in input TODO:make automatic
    if not args.useChains and len(st) > 1:
        print ("#WARNING: Input PDB contains more than one chain ids ",\
               "when chains ids not used, consider renumbering ")

#==============================================================================
#Detecting Covalent Pairs
    bckats = []
    for at in st.get_atoms():
        if '.'+at.id in bnsTopLib.StructureWrapper.backbone_link_atoms:
            bckats.append(at)

    nbsearch = NeighborSearch(bckats)

    covLinkPairs = []

    for at1, at2 in nbsearch.search_all(COVLNK):
        if sameResidue(at1,at2):
            continue
        covLinkPairs.append (getOrderedResiduePair(at1,at2,args.useChains))
    if args.debug:
        print ("#DEBUG: covalently linked residue pairs")
        for r in sorted(covLinkPairs, key=lambda i: i[0].residue.index):
            print ("#DEBUG: ",r[0].resid(),r[1].resid())
# Building Chains
    chList = bnsTopLib.ResidueSet.ResidueSetList(covLinkPairs)
    print ("#INFO: Found ",chList.n," chain(s)")
    if args.json:
        jsondata.data['NOfChains']=chList.n
        
    for s in chList.getSortedSets():
        print ("#CH ",s)
        if args.json:
            jsondata.data['chains'].append(
                {
                'iniRes' : str(s.inir), 
                'finRes' : str(s.finr), 
                'sequence' : s.getSequence()
                }
            )

# Residue list
    print ("#INFO: Residue Ids List")
    i=0
    for s in chList.getSortedSets():
        print ('#CH',s.inir,'-',s.finr,':', ','.join(s.getResidueIdList()))
        if args.graphml:
            for res in s.items:
                xml.addResidue(i,res)
        i=i+1

    if args.graphml:
        for r in sorted(covLinkPairs, key=lambda i: i[0].residue.index):
            xml.addBond('ch',r[0],r[1])
    if args.json:
        for r in sorted(covLinkPairs, key=lambda i: i[0].residue.index):
            jsondata.data['covLinks'].append([r[0].resid(True),r[1].resid(True)])
# Contacts & interfaces
    if args.contacts:
        print ("#INFO: Getting interchain contacts")
        conts={}
        intList={}
        interfPairs={}
        for ch1 in chList.getSortedSets():
            conts[ch1]={}
            intList[ch1]={}
            interfPairs[ch1]={}
            for ch2 in chList.getSortedSets():
                if ch2.ini <= ch1.ini:
                    continue
                conts[ch1][ch2]=[]
                interfPairs[ch1][ch2]=[]
                int
                ats1 = []
                s1 = ch1._getResidues()
                s2 = ch2._getResidues()
                for r1 in s1:
                    for r2 in s2:
                        cont = s1[r1].getClosestContact(s2[r2],INTDIST)
                        if cont:
                            [at1,at2,d] = cont
                            [atom1,atom2] = getOrderedAtomPair(at1,at2,args.useChains)
                            [res1,res2] = getOrderedResiduePair(at1,at2,args.useChains)
                            if d <= HBLNK:
                                print ("#CT ", atom1.atid(True), atom2.atid(True),d)
                                conts[ch1][ch2].append([atom1,atom2,d])
                                if args.json:
                                    jsondata.data['contacts'].append({'ats':[atom1.atid(True), atom2.atid(True)],'distance':float(d)})
                                interfPairs[ch1][ch2].append([res1,res2])
                intList[ch1][ch2] = bnsTopLib.ResidueSet.ResidueSetList(interfPairs[ch1][ch2])
        
        print ("#INFO: Interface residues at "+ str(INTDIST) +"A")
        for ch1 in chList.getSortedSets():
            for ch2 in chList.getSortedSets():
                if ch2.ini <= ch1.ini:
                    continue
                for s in intList[ch1][ch2].getSortedSets():
                    print ("#INTRES ", ','.join(s.getResidueIdList()))
                    if args.json:
                        jsondata.data['interfaces'].append(','.join(s.getResidueIdList()))
                
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

        [atom1,atom2] = getOrderedAtomPair(at1,at2,args.useChains)

        res1 = bnsTopLib.StructureWrapper.Residue(atom1.at.get_parent(),args.useChains)
        res2 = bnsTopLib.StructureWrapper.Residue(atom2.at.get_parent(),args.useChains)

        if [atom1.attype(), atom2.attype()] not in bnsTopLib.StructureWrapper.hbonds['WC'] and \
           [atom2.attype(), atom1.attype()] not in bnsTopLib.StructureWrapper.hbonds['WC']:
            continue

        if [res1,res2] in covLinkPairs:
            continue

        if res1 not in nhbs:
            nhbs[res1] = {}

        if res2 not in nhbs[res1]:
            nhbs[res1][res2]=0

        nhbs[res1][res2]=nhbs[res1][res2]+atom1._hbscore(atom2)

        if args.debug:
            print ("#DEBUG: ", res1,res2,atom1._hbscore(atom2))

    if args.debug:
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
        if maxv > args.bpthres:
            bps.append(bnsTopLib.StructureWrapper.BPair(r1,pair,maxv))

    bpsref={}

    for bp in sorted(bps):
        print ("#BP ", bp, '(',bp.type,'):',','.join(bp.compsIds()), '(',str(bp.score),')')
        bpsref[bp.r1.residue.index]=bp
        if args.graphml:
            xml.addBond("bp",bp.r1, bp.r2)
        if args.json:
            jsondata.data['bpList'].append({'id': bp.bpid(),'type':bp.type,'score':float(bp.score), 'comps':bp.compsIds()})
    
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

    for bpstp in sorted(bpsteps):
        print ("#BPST ", bpstp, '(',bpstp.type,'):', ','.join(bpstp.compsIds()))
        bpstpref[bpstp.bp1.r1.residue.index]=bpstp
        if args.json:
            jsondata.data['bpStepList'].append({'id':bpstp.stepid(),'type':bpstp.type,'comps':bpstp.compsIds()})
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
        print ("#HF ", ','.join(seq))
        if args.json:
            jsondata.data['HelicalFrags'].append(seq)
    if args.graphml:
        xml.save(args.graphml)
    if args.json:
        jsondata.save(args.json)
        
if __name__ == "__main__":
    main()
