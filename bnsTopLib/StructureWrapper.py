#
# Convenient wrappers over BioPython classes
# Hint: Include here interaction energies
#

__author__ = "gelpi"
__date__ = "$27-oct-2017 13:59:14$"

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

class Residue():
    oneLetterResidueCode = {
        'DA' :'A', 'DC' :'C', 'DG' :'G', 'DT' :'T',
        'A'  :'A', 'C'  :'C', 'G'  :'G', 'U'  :'U',
        'DA3':'A', 'DC3':'C', 'DG3':'G', 'DT3':'T',
        'A3' :'A', 'C3' :'C', 'G3' :'G', 'U3' :'U',
        'DA5':'A', 'DC5':'C', 'DG5':'G', 'DT5':'T',
        'A5' :'A', 'C5' :'C', 'G5' :'G', 'U5' :'U',
        'MRA':'A',
        'ALA':'A', 'CYS':'C', 'ASP':'D', 'GLU':'E', 'PHE':'F', 'GLY':'G',
        'HIS':'H', 'HID':'H', 'HIE':'H', 'ILE':'I', 'LYS':'K', 'LEU':'L',
        'MET':'M', 'ASN':'N', 'PRO':'P', 'GLN':'Q', 'ARG':'R', 'SER':'S',
        'THR':'T', 'VAL':'V', 'TRP':'W', 'TYR':'Y'
    }

    def __init__(self, r, useChains=False):
        self.residue = r
        self.useChains = useChains
        self.chain = r.get_parent().id
        self.resNum = int(self.residue.id[1])

    def resid(self, compact=False):
        if self.useChains:
            ch = ":"+self.chain
        else:
            ch = ''
        if compact:
            return self._getOneLetterResidueCode() + ch + str(self.resNum)
        else:
            return self.residue.get_resname() + ch + ':'+ str(self.resNum)

    def bnsid(self):
        return self.chain+str(self.resNum)+"-"+self._getOneLetterResidueCode()

    def _getOneLetterResidueCode(self):
        resid = self.residue.get_resname().rstrip().lstrip()
        if not resid in Residue.oneLetterResidueCode:
            return 'X'
        else:
            return Residue.oneLetterResidueCode[resid]

    def __hash__(self):
        return hash(self.resid())

    def __eq__(self, other):
        return self.resid() == other.resid()

    def __lt__(self, other):
        if self.chain != other.chain:
            return self.chain < other.chain
        else:
            return self.resNum < other.resNum

    def __str__(self):
        return self.resid()
    
    def __index__(self):
        if self.useChains:
            return ord(self.chain)*1000 + int(self.resNum)
        else:
            return self.resNum
    

class Atom():
    def __init__ (self, at, useChains=False):
         self.at=at
         self.useChains=useChains

    def atid(self, compact=False):
        return self.resid(compact)+"."+self.at.id

    def resid(self, compact=False):
        return Residue(self.at.get_parent(),self.useChains).resid(compact)

    def resNum(self):
        return Residue(self.at.get_parent(),self.useChains).resNum()

    def attype(self):
        return Residue(self.at.get_parent(),self.useChains)._getOneLetterResidueCode()+'.'+self.at.id

    def __lt__(self,other):
        return self.at.get_serial_number()<other.at.get_serial_number()

    def __str__(self):
        return self.atid()

    def _hbscore(self,other):
        d = self.at - other.at
        return 2.6875 - 0.625*d

class BPair():
    def __init__(self,r1,r2,score):
        self.r1=r1
        self.r2=r2
        types =[self.r1._getOneLetterResidueCode(), self.r2._getOneLetterResidueCode()]
        self.type=''.join(sorted(types))
        self.score=score

    def bpid(self):
        return self.r1.chain + str(self.r1.resNum) + "-" \
            + self.r1._getOneLetterResidueCode() \
            + self.r2._getOneLetterResidueCode()

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
        bps = [
            self.bp1.r1._getOneLetterResidueCode() + self.bp2.r1._getOneLetterResidueCode(),
            self.bp2.r2._getOneLetterResidueCode() + self.bp1.r2._getOneLetterResidueCode()
        ]
        self.type= ''.join(sorted(bps))

    def stepid(self):
        return str(self.bp1.r1.resNum) + "-" \
            + self.bp1.r1._getOneLetterResidueCode() \
            + self.bp2.r1._getOneLetterResidueCode() \
            + self.bp2.r2._getOneLetterResidueCode() \
            + self.bp1.r2._getOneLetterResidueCode()

    def comps(self):
        return [self.bp1.bpid(),self.bp2.bpid()]

    def resNum(self):
        return self.bp1.r1.resNum

    def __eq__(self,other):
        return self.bp1==other.bp1 and self.bp2 == other.bp2

    def __lt__(self,other):
        return self.bp1<other.bp1

    def __str__(self):
      return self.stepid()

    def __hash__(self):
        return hash(self.stepid())
