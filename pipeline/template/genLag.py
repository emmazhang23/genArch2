
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
from optparse import OptionParser
import argparse


# In[ ]:


parser = argparse.ArgumentParser(description='Generate Lag File')

parser.add_argument("-n", dest="filename", 
                 help="name of simulation")
parser.add_argument("-l", dest="loci", action="store",
                 help="specify the number of loci")
parser.add_argument("-r", dest="replicates", action="store", default=1,
                 help="specify replicates of the simulation")

args = parser.parse_args()


# In[2]:


def lagphase(data):
    lag = []
    rep_counter = 1
    row = 0
    while rep_counter <= reps:
        if rep_counter == data.iloc[row,0]:
            for line in range(520):
                if data.iloc[line+row,4] == 1000:
                    lag.append(data.iloc[line+row,1]-10000)
                    rep_counter = rep_counter+1
                    break
                elif line == 519:
                    lag.append(500)
                    rep_counter = rep_counter+1
                    break
        elif rep_counter != data.iloc[row,0]:
            while rep_counter != data.iloc[row,0]:
                row = row+1
    return lag


# In[7]:


loci = args.loci
reps = args.replicates

data = pd.read_csv("./"+args.filename+"/stats/"+args.filename+"_stats.txt", delim_whitespace = True)
master = pd.DataFrame({'rep': range(reps)})
master.iloc[:,0] = master.iloc[:,0] + 1

#### Compute and output lag phase
inv_lag = lagphase(data)
master = pd.concat([master, pd.DataFrame({''.join(['lag phase duration']): inv_lag})], axis=1)


# In[8]:

master.loc["Total"] = master.loc[master['lag phase duration']==500].count()

master.to_csv('./'+args.filename+'_lag.csv')

