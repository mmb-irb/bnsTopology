# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.


class ChainList():
    def __init__(self):
        self.chains=[]
        self.n=0

    def find(self,item):
        i = 0
        while i < self.n and item not in self.chains[i].items:
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
        return sorted(self.chains,key=lambda s: s.ini)

class Chain():
    def __init__(self):
        self.ini=999999
        self.fin=0
        self.iniId=''
        self.finId=''
        self.items = set()

    def add(self,r):
        self.items.add(r)
        rn=r.__index__()
        if rn < self.ini:
            self.ini = rn
            self.inir = r.chain+str(r.resNum)
        if rn > self.fin:
            self.fin = rn
            self.finr = r.chain+str(r.resNum)
        
    def union(self,other):
        self.items = self.items.union(other.items)
        if other.ini < self.ini:
            self.ini  = other.ini
            self.inir = other.inir
        if other.fin > self.fin:
            self.fin = other.fin
            self.finr = other.finr

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
            seql.append(seq[i].chain+str(seq[i].resNum)+"-"+seq[i]._getOneLetterResidueCode())
        return seql

    def _getResidues(self):
        seq={}
        for r in self.items:
            seq[r.__index__()] = r
        #print (seq)
        return seq
    
    def __str__(self):
        return str(self.inir) + "-" + str(self.finr)+ ":" + self.getSequence()
    

class BPChain (Chain):
    def add(self,bpst):
        self.items.add(bpst)
        rn=bpst.bp1.r1.residue.index
        if rn < self.ini:
            self.ini = rn
            self.inir = bpst
        if rn > self.fin:
            self.fin = rn
            self.finr = bpst
    
    def getSequence(self):
        seq=self._getResidues()
        ss=''
        for i in sorted(seq.keys()):
            print (seq[i].bp1)
            ss=ss+seq[i].bp1.bpid()+","
        return ss
    
    def _getResidues(self):
        seq={}
        for bpst in self.items:
            seq[bpst.bp1.r1.residue.index] = bpst
        return seq
    
    def __str__(self):
        return str(self.inir.bp1.r1.resid(1)) + "-" + str(self.finr.bp1.r1.resid(1))+ ":  " + self.getSequence()
    
