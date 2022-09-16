
# coding: utf-8

# In[1]:


import pandas as pd
import copy
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pandas import DataFrame, Series
from matplotlib.patches import Polygon
import argparse


# In[ ]:


parser = argparse.ArgumentParser(description='Generate genotype over time plot')

parser.add_argument("-n", dest="filename", type=str,
                 help="name of simulation")

parser.add_argument("-l", dest="loci", action="store", type=int, help="number of loci")

parser.add_argument("-le", dest="epiloci", action="store", type=int,
                 help="specify the number of epistatic loci")


parser.add_argument("-r", dest="replicates", action="store", default=1, type=str,
                 help="specify replicates of the simulation")

parser.add_argument("-b", dest="burn", action="store", default=10000,
                 help="specify replicates of the simulation")

parser.add_argument("-d", dest="dispersal", action="store", default=500,
                 help="specify replicates of the simulation")

parser.add_argument("-q", action="store_false", dest="verbose", default=True,
                 help="don't print status messages to stdout")


args = parser.parse_args()


# In[2]:


reps=args.replicates
gens=(args.burn+args.dispersal)
path="./" +args.filename+"/quanti_genotype/"
loci= int(args.loci)
eloci=int(args.epiloci)

# In[3]:


repStr= str(reps)
repLen=len(repStr)

genStr=str(gens)
genLen=len(genStr)


# In[44]:


#Make dataframe with
#replicate      Generation         DMI genotype frequency
#plot by facet of Replicate so each line is a replicate
column_names = ["rep", "gen", "freq"]

df = pd.DataFrame(columns = column_names)
for i in range(1, int(reps)+1, 1):
    for n in range(10000, int(gens)+1, 1):
        #table=pd.read_table(path+args.filename+"_g"+n.zfill(genLen)+"_r"+i.zfill(repLen)+".dat", header=None, skiprows=(2*loci)+1+eloci, delim_whitespace = True)
        table=pd.read_table(path+args.filename+"_g"+str(n).zfill(genLen)+"_r"+str(i).zfill(repLen)+".dat", header=None, skiprows=(2*loci)+1+eloci, delim_whitespace = True)
        table=table.loc[table[0] == 3]
        Dfreq=0 #dmi freq
        Nfreq=0 #normal freq
        for p in range(1, (loci*2)+1):
            table=table.drop([p], axis=1)


        table[(loci*2)+eloci+1]='0'

        table=table.astype(str)

        for p in range((loci*2)+1,1+eloci+(loci*2)):
            #table[loci+1+eloci]=table[loci+1+eloci]+isSubstring('2', table[p])

            table[(loci*2)+1+eloci]+=np.where(table[p].str.contains('2'), '2', '1')

        freq=((table[(loci*2)+1+eloci].str.contains('22')).sum())/len(table.index)
        row = {'rep': i, 'cgen': n, 'freq': freq}
        df = df.append(row, ignore_index = True)
#freq of combined probability for each generation
#at end
#offpsring in patch 3
#time where migration starts what does patch 1 and 2 look like
#does DMI frequency go up or down

#Calculate genotype frequency for each generation over time


# In[ ]:


df.to_csv("genoTime_"+args.filename+".csv")
