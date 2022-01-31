import pandas as pd
import matplotlib.pyplot as plt
from uszipcode import SearchEngine
import mpu
#for extensive list of zipcodes, set simple_zipcode =False
search = SearchEngine(simple_zipcode=False)

df = pd.read_csv('ziplist.csv')
df['center'] = df.zipcode.apply(lambda x: (search.by_zipcode(x).lat, search.by_zipcode(x).lng))
df['lat'] = df.center.apply(lambda x: x[0])
df['long'] = df.center.apply(lambda x: x[1])
df.loc[:, ['name', 'zipcode', 'lat', 'long']].to_csv('latlong.csv', index=False)