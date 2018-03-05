
import json

class JSONWriter:
    def __init__(self):
        self.data = {}
        
    def insert(self,k,v):
        self.data[k]=v
    
    def write(self):
        return json.JSONEncoder(sort_keys=True, indent=1).encode(self.data)
        
    def save(self,file):
        jsout = open(file,"w+")
        jsout.write(self.write())
        jsout.close()


