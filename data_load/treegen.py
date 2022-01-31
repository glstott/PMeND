import pandas as pd
import dendropy
import numpy as np

# loop through directory with Newick Files, then generate a dataframe of edges.
path = "/scratch/gs69042/bahllab/"
combined_df = pd.DataFrame()
for file in os.listdir(path):
    if '.nwk' in file:
        tr1 = ph.read_newick(path+file)
    else:
        continue
    print(file)
    parents = tr1.parent[1:].apply(int).sort_values(ascending=False).unique()
    for parent1 in parents:
        parent = str(parent1)
        if parent:
            parentid = "|".join(tr1[tr1.parent == parent].id.sort_values().apply(lambda x: x if "|" in x else x.split('/')[2]))

            tr1.loc[tr1.parent == parent, 'parent'] = parentid
            tr1.loc[tr1.id == parent, 'id'] = parentid
    tr1.loc[tr1.type == 'leaf', 'id'] = tr1.loc[tr1.type == 'leaf', 'id'].apply(lambda x: x.split("/")[2])
    tr1['source'] = file
    print(tr1.head())
    if len(combined_df) > 1:
        combined_df = pd.concat([combined_df, tr1])
    else:
        combined_df = tr1
    print("-------------------------")

# output to file
combined_df.to_csv('neo_ready.csv', index=False)