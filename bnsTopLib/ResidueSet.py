
import re

class ResidueSetList():
    def __init__(self, pairsList=None, debug=False):
        self.sets=[]
        self.n=0
        if (pairsList):
            self.addData(pairsList,debug)

    def addData(self, pairsList, debug=False):
        for [res1,res2] in pairsList:        
            i = self.find(res1)
            j = self.find(res2)
            if i == -1 and j == -1:
                s = ResidueSet()
                s.add(res1)
                s.add(res2)
                self.append(s)
            elif i != -1 and j != -1 and i != j:
                self.sets[i].union(self.sets[j])
                self.delete(j)
            elif j == -1:
                self.sets[i].add(res2)
            elif i == -1:
                self.sets[j].add(res1)
            if debug:
                for s in self.getSortedSets():
                    print ("#DEBUG:" ,s)
    
    def find(self,item):
        i = 0
        while i < self.n and item not in self.sets[i].items:
            i = i + 1
        if i == self.n:
            return -1
        else:
            return i

    def append(self,item):
        self.sets.append(item)
        self.n = len(self.sets)

    def delete(self,i):
        del self.sets[i]
        self.n = len(self.sets)

    def getSortedSets(self):
        return sorted(self.sets,key=lambda s: s.ini)

class ResidueSet():
    def __init__(self):
        self.ini=999999
        self.fin=0
        self.iniId=''
        self.finId=''
        self.items = set()
        self.type=''
        self.id=''

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
    
    def isProtein(self):
        seq = self.getSequence()
        seq = re.sub('[ACTGUX]','',seq)
        if seq:
            return True
        else:
            return False

    def _getResidues(self):
        seq={}
        for r in self.items:
            seq[r.__index__()] = r
        #print (seq)
        return seq
    
    def __str__(self):
        return str(self.id) + ": " +str(self.inir) +  "-" + str(self.finr)+ "("+self.type+"):" + self.getSequence()
    

class BPSSetList (ResidueSetList):
    def addData(self, pairsList, debug=False):
        for [res1,res2] in pairsList:        
            i = self.find(res1)
            j = self.find(res2)
            if i == -1 and j == -1:
                s = BPSSet()
                s.add(res1)
                s.add(res2)
                self.append(s)
            elif i != -1 and j != -1 and i != j:
                self.sets[i].union(self.sets[j])
                self.delete(j)
            elif j == -1:
                self.sets[i].add(res2)
            elif i == -1:
                self.sets[j].add(res1)
            if debug:
                for s in self.getSortedSets():
                    print ("#DEBUG:" ,s)
    
class BPSSet (ResidueSet):
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
            ss=ss+seq[i].bp1.bpid()+","
        return ss
    
    def _getResidues(self):
        seq={}
        for bpst in self.items:
            seq[bpst.bp1.r1.residue.index] = bpst
        return seq
    
    def __str__(self):
        return str(self.inir.bp1.r1.resid(1)) + "-" + str(self.finr.bp1.r1.resid(1))+ ":  " + self.getSequence()
    
