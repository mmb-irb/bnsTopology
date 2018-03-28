#! /usr/bin/python
#
# Simple parser to extact topology of NA from simulation snapshot
# Chain ids or TER not required
#
__author__ = "gelpi"
__date__ = "$29-ago-2017 16:14:26$"

import bnsTopLib

def main():

    args = bnsTopLib.cmdLine({'BPTHRESDEF': bnsTopLib.Topology.BPTHRESDEF, 'INTDIST':bnsTopLib.Topology.INTDIST, 'LIMIT':10000}).parse_args()
# ============================================================================
    loader = bnsTopLib.PDBLoader(args.pdb_path)
    
    args.useChains= args.useChains or loader.useChains      
    
    if not args.id and loader.id:
        args.id = loader.id
    args.format = loader.format
    print()
    st = loader.loadStructure()
    if args.useChains:
        print ('#INFO: Using chain ids from input data')
    else:
        print ('#INFO: Not using chain ids')
#Checking for chain ids in input TODO:make automatic
        if len(st) > 1:
            print ("#WARNING: Input PDB contains more than one chain ids, consider renumbering ")
# Checking max atoms
    if args.limit:
        if loader.numAts > args.limit:
            print ("#LIMIT atoms exceeded ("+str(loader.numAts)+")")
            sys.exit(1)
#==============================================================================
# Initializing topology object
    top = bnsTopLib.Topology(args)
#Detecting Covalent Pairs
    top.calcCovLinks(st)
# Building Chains
    top.chList = bnsTopLib.ResidueSet.ResidueSetList(top.covLinkPairs)
    
    print()
    print ('#INFO: Found %i chain(s)'% (top.chList.n))
    
    for s in top.chList.getSortedSets():
        print ("#CH", s)

# Residue list
    print ("\n#INFO: Residue Ids List")
    
    i=0
    for s in top.chList.getSortedSets():
        print ('#CH', s.__str__(1))
        i=i+1
# Contacts & interfaces
    if args.contacts or args.interface:
        top.calcContacts()

        if args.contacts:
            print()
            print ("#INFO: Getting interchain contacts")
            for ch1 in top.conts:
                for ch2 in top.conts[ch1]:
                    for c in top.conts[ch1][ch2]:
                        print ("#CT", '%i-%i %-10s %-10s %6.4f' % (ch1.id, ch2.id, c[0].atid(True), c[1].atid(True),c[2]))
        
        if args.interface:
            print()
            print ("#INFO: Interface residues at "+ str(args.intdist) +"A")
            for ch1 in top.intList:
                for ch2 in top.intList[ch1]:
                    for s in top.intList[ch1][ch2].getSortedSets():
                        print ('#INTRES %i-%i %s' % (ch1.id, ch2.id,','.join(s.getResidueIdList())))
                
    if not top.checkExistsNA():
        print()
        print ("#WARNING: No WC atom(s) found, skipping NA topology")
#Base Pairs
    else:
        print()
        print ("#INFO: Base pairs found")
        top.calcBasePairs()
        for bp in sorted(top.bps):
            print ("#BP ", bp.__str__(1))
    # Bpair steps from neighbour bps, relays on residue renumbering. TODO Check connectivity
        print()
        print ("#INFO: Base Pair steps")
        top.calcBPSteps()
        for bpstp in sorted(top.bpsteps):
            print ("#BPST", bpstp.__str__(1))

# Continuous helical segments from stretches of overlapping bsteps
        print()
        print ("#INFO: Helical segments")
        top.calcHelicalFrags()
        for frs in top.HFSeqs:
            print ("#HF", ','.join(frs))
        print ()
        print ("#INFO: Residues not in helical segments")
        nhf =[]
        for r in top.notInHF:
            nhf.append(r.resid(1))
        print ("#NHF",','.join(nhf))
            
    if args.json:
        top.json().save(args.json)
        print()
        print ("#INFO: JSON data written on "+args.json)

    if args.graphml:
# Graphml output
        xml = bnsTopLib.GraphmlWriter()
        for s in top.chList.getSortedSets():
            for res in s.items:
                xml.addResidue(i,res)
        for r in sorted(top.covLinkPairs, key=lambda i: i[0].residue.index):
            xml.addBond('ch',r[0],r[1])
        if top.NAOk:
            for bp in sorted(top.bps):
                xml.addBond("bp",bp.r1, bp.r2)
        xml.save(args.graphml)
        print ("#INFO: GRAPH data written on "+args.graphml)
        
if __name__ == "__main__":
    main()
