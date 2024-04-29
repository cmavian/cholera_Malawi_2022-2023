"""
Created on Tue 9 April 2019 

@authors: david rasmussen
@edited: eduan wilkinson
"""
import baltic as bt
import pandas as pd

treeFile = "annotated_tree.nexus"
outFile = 'annottated_tree_events.csv'
myTree=bt.loadNewick(treeFile, absoluteTime=False)
myTree.setAbsoluteTime(2023.128767) # need to set this to time of last sampled tip

myTree.traverse_tree() ## required to set heights
myTree.treeStats() ## report stats about tree

changes = 0
times = []
origins = []
destinations = []
for k in myTree.Objects: ## iterate over a flat list of branches
    
    "Assign node UNKNOWN district if not give"
    if 'district' in k.traits:
        district = k.traits['district']
    else:
        district = 'UNKNOWN'
        k.traits['district'] = district 
    
    "Find parent district if given"
    if k.parent.traits:
        parent_district = k.parent.traits['district']
    else:
        parent_district = 'UNKNOWN'
        
    if district != parent_district:
        changes += 1
        times = times + [k.absoluteTime]
        origins = origins + [parent_district]
        destinations = destinations + [district]
        
        
print("Total number of state changes: " + str(changes))

df = pd.DataFrame({'EventTime' : times})
df['Origin'] = origins
df['Destination'] = destinations

df.to_csv(outFile)  
    

