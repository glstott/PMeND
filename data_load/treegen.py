import pandas as pd
import dendropy
import numpy as np
import sys
import phylopandas as ph

# generate a dataframe of edges.
path = str(sys.argv[1])
if path[-4:] == ".nex":
    tr1 = ph.read_nexus_tree(path)
else:
    tr1 = ph.read_newick(path)
print(path)

# Find list of unique parents. For each, generate a unique long id from concatenation of children.
#   This will likely run out of space. I need to come up with a better solution.
tr1.id = tr1.id.apply(lambda x: x.split("|")[0] if "|" in x else x)
parents = tr1.parent[1:].apply(int).sort_values(ascending=False).unique()
for parent1 in parents:
    parent = str(parent1)
    if parent:
        parentid = "|".join(tr1[tr1.parent == parent].id.sort_values().apply(lambda x: x if not "/" in x else x.split('/')[2]))

        tr1.loc[tr1.parent == parent, 'parent'] = parentid
        tr1.loc[tr1.id == parent, 'id'] = parentid

# Shorten the names of the leaves and save source info. In the future, we should move the name changes to outside this script in case someone wants to use this pipeline for non-SC2 sequences or sequences not from GISAID.
tr1.loc[tr1.type == 'leaf', 'id'] = tr1.loc[tr1.type == 'leaf', 'id'].apply(lambda x: x.split("/")[2])
tr1['source'] = path
print(tr1.head())
print("-------------------------")

# output to file
tr1.to_csv(str(sys.argv[2]), index=False)