from epiweeks import Week, Year
import pandas as pd
from Bio import SeqIO
import sys
file = str(sys.argv[1])
file2 = str(sys.argv[2])
prefix = file[:-13]
df = pd.read_csv(file, sep='\t')
mode='cumulative'
df.loc[len(df), ['strain', 'date']] = ["hCoV-19/Wuhan/WIV04/2019|EPI_ISL_402124|2019-12-30", "2019-12-30"]

df['epiweek'] = df.date.map(lambda x: Week.fromstring(x).cdcformat())
for epiweek in df.epiweek.unique():
    print(epiweek)
    if mode == 'cumulative':
        temp_df = df.loc[df.epiweek <= epiweek, :]
    else:
        temp_df = df.loc[df.epiweek == epiweek, :]
    if temp_df.strain.count() < 100:
        continue
    temp_df.loc[:, ['strain', 'date']].to_csv(prefix + "_" +str(epiweek)+".dt.tsv", index=False, sep='\t', header=False)
    temp_df.to_csv(prefix + "_" +str(epiweek)+".metadata.tsv", index=False, sep='\t' )
    with open(prefix + "_" + str(epiweek)+".fasta", 'a') as out_handle:
        for record in SeqIO.parse(file2, 'fasta'):
            if record.id in temp_df.strain.unique():
                SeqIO.write(record, out_handle, "fasta")
    print(temp_df.count()[0])
