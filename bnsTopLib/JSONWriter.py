
import json

class JSONWriter:
    def __init__(self):
        self.data = {
            'id':'',
            'useChains':'',
            'NOfChains':0, 
            'chains':[],
            'covLinks' : [],
            'contacts' : [],
            'interfaces': [],
            'bpList': [],
            'bpStepList' : [],
            'HelicalFrags': []
        }
         
    def __str__(self):
        return json.JSONEncoder(sort_keys=True, indent=1).encode(self.data)
        
    def save(self,file):
        jsout = open(file,"w+")
        jsout.write(self.__str__())
        jsout.close()


