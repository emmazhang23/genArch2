
# coding: utf-8

# In[239]:


import pandas as pd
import numpy as np
import argparse
import itertools
from re import search
import random as random
import re as re


# In[292]:

#script to generate the epistasis file that specifies DMI genotypes and the fitness effect
#Allow for specification of filename, number of DMI loci, and effect size in command line

my_parser = argparse.ArgumentParser(description='Generate Epistasis File')
my_parser.add_argument('-l',
                       metavar='loci',
                       dest="loci",
                       type=int,
                       help='the number of loci')
my_parser.add_argument('-n',
                       metavar='filename', dest="filename",
                       help='filename')
my_parser.add_argument('-d',
                       metavar='decrease to fitness',
                       dest="decrease",
                       type=float,
                       help='decrease to fitness')

args = my_parser.parse_args()


# In[273]:


loci= args.loci


# In[274]:

#open new file with same naming convention with _epi.txt extension 
f= open(args.filename +"_epi.txt","w")


# In[275]:

#write file headers in quantinemo convention
f.write("[FILE_INFO] {\n")
f.write("col_genotype 1\n")
f.write("col_fitness_factor 2\n")
f.write("}\n")
f.write("#genotype fitness\n")


# In[276]:

#Prupose of script:
#generate table of simplified mutation (all the combinations of 2 in pairs) and their fitness values (if dmi, decrease)
#generated table should be "loci" number of consecutive numbers, 2 options per loci (because 2 alleles)
#generate table of all permutation's (simplified permutations) fitnesses
#go through each permutation and add up all partial matches

# In[277]:

#helper function to find index of character in string
def find(s, ch):
    return [i for i, ltr in enumerate(s) if ltr == ch]


# In[278]:

#helper function to determine whether Array A is a subarray of array B
def isSubArray(A, B, n, m):
    i = 0; j = 0;

    while (i < n and j < m):

        if (A[i] == B[j]):

            i += 1;
            j += 1;

            if (j == m):
                return True;

        else:
            i = i - j + 1;
            j = 0;

    return False;


# In[279]:

#Create a column in df with all permutations of ancestral allele (1) and novel allele (2) as a string of "loci" 
simple = ["1", "2"]
simp = pd.DataFrame(list(itertools.product(simple, repeat = loci)))
simp[loci]=simp[0]


for k in range(1,loci,1):
    simp[loci]+=simp[k]
simp[loci+1]=0
simp[loci+2]=1


# In[280]:

series = []
for index, row in simp.iterrows():
    l=0
    length=simp.loc[index, loci].count("2")
    array=np.empty(length)

    n = find(row[loci], "2")
    series.append(n)


# In[281]:



simp[loci+1]=series


##pairwise interaction fitness set

for index, row in simp.iterrows():
    if (str(22) in row[loci]) and (row[loci].count("2")==2):
        ##where dmi intensity is calculated, will change second number (1-x) for matrix
        simp.loc[index, loci+2]= 1-args.decrease


##multiple pairwise fitness
##probably need to check if this works for 3 loci
for index, row in simp.iterrows():
    if (row[loci].count("2")>2):
        fit=0
        for index1, row1 in simp.iterrows():
            if index != index1:
                if (isSubArray(row[loci+1], row1[loci+1], len(row[loci+1]), len(row1[loci+1]))):
                    fit+=(1-(1-(row1[loci+2])))

            simp.loc[index, loci+2]=1-fit
        if (simp.loc[index, loci+2]<0):
            simp.loc[index, loci+2]=0


# In[282]:


simp_dict=dict(zip(simp[loci], simp[loci+2]))


# In[283]:


lists = ["0101", "0102", "0202"]
perms = pd.DataFrame(list(itertools.product(lists, repeat = loci)))


# In[284]:


perms[loci]="0"
for index, row in perms.iterrows():
    if (str(row[0]).count("2"))>0:
        row[loci] = "2"
    else:
        row[loci]="1"

for index, row in perms.iterrows():
    for x in range(1, loci):
        if (str(row[x]).count("2"))>0:
            row[loci] += "2"
        else:
            row[loci]+="1"


# In[285]:


perms[loci+1]= perms[loci].map(simp_dict)


# In[286]:


for index, row in perms.iterrows():
    line="{"
    for x in range(loci):
        line+=row[x]
        if x != (loci-1):
            line+=" "
    line+="} "+str(row[loci+1])+"\n"
    f.write(line)


# In[287]:


f.close()
