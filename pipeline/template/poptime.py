# coding: utf-8
# %%

# %%
#generate csv with population in each patch over time after dispersal

import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import argparse

# %%


my_parser = argparse.ArgumentParser(description='Generate Pop over Time csv')
my_parser.add_argument('-l',
                       metavar='loci',
                       dest="loci",
                       type=str,
                       help='the number of loci')
my_parser.add_argument('-s',
                       metavar='selection', dest="selection",
                       help='selection')
my_parser.add_argument('-d',
                       metavar='decrease', dest="decrease",
                       help='decrease to f')

args = my_parser.parse_args()


# %%

loci = args.loci
decrease = args.decrease
selection = args.selection

filename=str(loci)+'l_'+str(decrease)+'df_'+str(selection)+'s'


# %%

path = './'+filename+'/stats/'

df=pd.read_csv(path+filename+'_stats.txt', delimiter='\t')
df=df.iloc[:, [0,1,4]]

df = df.rename(columns={"replicate   ": "replicate", "generation  ": "generation", "adlt.nbInd_p3": "population"})
#extract only dispersal phase generations
df= df[df['generation'] > 10000] 


# %%

#save to csv
df.to_csv(filename+"_popsizes.csv")
#-----------------------------------------------------------------
#make plot

