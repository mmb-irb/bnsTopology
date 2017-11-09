# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.


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

