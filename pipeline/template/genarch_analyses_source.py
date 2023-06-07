
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

#%matplotlib inline


# In[ ]:

parser = OptionParser()
## options to enter in the command line
parser.add_option("-n", "--name", dest="filename",  action="store", type="string", 
                 help="name of simulation")
parser.add_option("-p", "--path", dest="path",  action="store", type="string", 
                 help="path to data directory")
parser.add_option("-l", "--loci", dest="loci",  action="store", type="int", 
                 help="number of loci")
parser.add_option("-r", "--reps", dest="reps",  action="store", type="int", 
                 help="simulation replicates")
parser.add_option("-f", "--minfrequency", dest="minfreq",  action="store", type="float", 
                 help="allele frequency cutoff")
parser.add_option("-q", "--quiet", action="store_false", dest="verbose", default=True,
                 help="don't print status messages to stdout")

(options, args) = parser.parse_args()


# In[2]:

#### generates a genetic architecture dataframe for a trait in a patch
def genarch(pop,generation,trait,loci,N):
    genos = pd.read_table(str(options.path)+str(options.filename)+'/quanti_genotype/'+str(options.filename)+'_g'+str(int(generation))+'_r'+str(rep_counter).zfill(3)+'.dat', header=None, skiprows=(2*loci)+1, delim_whitespace = True)
    alleles = pd.read_table(str(options.path)+str(options.filename)+'/allelic_values_t'+str(trait)+'_r'+str(rep_counter).zfill(3)+'.txt', header=None, skiprows=16, delim_whitespace = True)                    
    zfill = len(str(loci))
    row = 0
    genotypes = pd.DataFrame({'allele#': list(range(1,2*N))})
    while row < len(genos.index):
        if genos.ix[row,0] != pop:
            row = row+1
        elif genos.ix[row,0] == pop:
            popstart = row
            for locus in range(loci):
                row = popstart
                alpha = []
                for indv in range(N):
                    genotype = str(genos.ix[row,(locus+1+((trait-1)*loci))])
                    alpha.append(alleles.ix[int(genotype[:(int(len(genotype)/2))])-1,2])
                    alpha.append(alleles.ix[int(genotype[(int(len(genotype)/2)):int(len(genotype))])-1,2])
                    row = row+1 
                genotypes = pd.concat([genotypes, pd.DataFrame({''.join(['locus_',str(locus+1).zfill(zfill)]): alpha})], axis=1)
            break
    location = []
    alphas = []
    freqs = []
    for locus in range(loci):
        al = genotypes.ix[:,locus+1].unique()
        for allele in range(len(al)):
            location.append(locus+1)
            alphas.append(al[allele])
            freqs.append(genotypes[genotypes.ix[:,locus+1] == al[allele]].count()[genotypes.columns[locus+1]]/(2*N))
    arch = pd.DataFrame({'locus#': location, 'alpha': alphas, 'frequency': freqs})
    return arch


class population(object):
    def __init__(self,pop,generation,loci,N):
        """
        Initiate a population object with genetic architecture dataframes for the two quantitative traits
        """
        self.t1arch = genarch(pop,generation,1,loci,N)
        self.t2arch = genarch(pop,generation,2,loci,N)

def sourcepop(pop1,pop2,pop3,minfreq):
    """
    determine the source population of each adaptive allele in population 3 with frequency greater than 
    a minimum frequency and add source population ID to the pop3 architecture dataframe for each trait
    """
    #### for trait 1 
    alphas = []
    freqs = []
    source = []
    pop1_alphas = []
    pop1_loci = []
    pop2_alphas = []
    pop2_loci = []
    for allele in range(len(pop1.t1arch.ix[:,0])):
        if pop1.t1arch.ix[allele,0] > 0:
            pop1_alphas.append(pop1.t1arch.ix[allele,0])
            pop1_loci.append(pop1.t1arch.ix[allele,2])
    for allele in range(len(pop2.t1arch.ix[:,0])):
        if pop2.t1arch.ix[allele,0] > 0:
            pop2_alphas.append(pop2.t1arch.ix[allele,0])
            pop2_loci.append(pop2.t1arch.ix[allele,2])    
    for allele in range(len(pop3.t1arch.ix[:,0])):
        found = []
        if pop3.t1arch.ix[allele,1] > minfreq and pop3.t1arch.ix[allele,0] > 0:
            for donor in range(len(pop1_alphas)):
                if pop3.t1arch.ix[allele,0] == pop1_alphas[donor] and pop3.t1arch.ix[allele,2] == pop1_loci[donor]:
                    found.append(1)
                    break
            for donor in range(len(pop2_alphas)):
                if pop3.t1arch.ix[allele,0] == pop2_alphas[donor] and pop3.t1arch.ix[allele,2] == pop2_loci[donor]:
                    found.append(2)
                    break
            alphas.append(pop3.t1arch.ix[allele,0])
            freqs.append(pop3.t1arch.ix[allele,1])
            source.append(str(found))
    pop3.t1_source = pd.DataFrame({'alpha': alphas, 'frequency': freqs, 'source': source})
    pop3.t1_source1 = source.count('[1]')
    pop3.t1_source2 = source.count('[2]')
    pop3.t1_ambiguous = source.count('[1, 2]')
    #### for trait 2
    alphas = []
    freqs = []
    source = []
    pop1_alphas = []
    pop1_loci = []
    pop2_alphas = []
    pop2_loci = []
    for allele in range(len(pop1.t2arch.ix[:,0])):
        if pop1.t2arch.ix[allele,0] < 0:
            pop1_alphas.append(pop1.t2arch.ix[allele,0])
            pop1_loci.append(pop1.t2arch.ix[allele,2])
    for allele in range(len(pop2.t2arch.ix[:,0])):
        if pop2.t2arch.ix[allele,0] < 0:
            pop2_alphas.append(pop2.t2arch.ix[allele,0])
            pop2_loci.append(pop2.t2arch.ix[allele,2])    
    for allele in range(len(pop3.t2arch.ix[:,0])):
        found = []
        if pop3.t2arch.ix[allele,1] > minfreq and pop3.t2arch.ix[allele,0] < 0:
            for donor in range(len(pop1_alphas)):
                if pop3.t2arch.ix[allele,0] == pop1_alphas[donor] and pop3.t2arch.ix[allele,2] == pop1_loci[donor]:
                    found.append(1)
                    break
            for donor in range(len(pop2_alphas)):
                if pop3.t2arch.ix[allele,0] == pop2_alphas[donor] and pop3.t2arch.ix[allele,2] == pop2_loci[donor]:
                    found.append(2)
                    break
            alphas.append(pop3.t2arch.ix[allele,0])
            freqs.append(pop3.t2arch.ix[allele,1])
            source.append(str(found))                                                   
    pop3.t2_source = pd.DataFrame({'alpha': alphas, 'frequency': freqs, 'source': source})
    pop3.t2_source1 = source.count('[1]')
    pop3.t2_source2 = source.count('[2]')
    pop3.t2_ambiguous = source.count('[1, 2]') 


# In[ ]:

pops = 3
loci = options.loci
reps = options.reps
minfreq = options.minfreq

#### initiate master dataframe
data = pd.read_table(str(options.path)+str(options.filename)+'/stats/'+str(options.filename)+'_stats.txt', delim_whitespace = True)
master = pd.DataFrame({'rep': range(reps)})
master.ix[:,0] = master.ix[:,0] + 1
lagphase = []

#### intiate lists for trait architecture metrics
t1_pop1_weight_mean = []
t1_pop1_weight_sum = []
t1_pop2_weight_mean = []
t1_pop2_weight_sum = []
t2_pop1_weight_mean = []
t2_pop1_weight_sum = []
t2_pop2_weight_mean = []
t2_pop2_weight_sum = []
t1_pop1_alleles = []
t1_pop2_alleles = []
t1_ambig_alleles = []
t2_pop1_alleles = []
t2_pop2_alleles = []
t2_ambig_alleles = []

#### calculate trait architecture metrics for each simulation replicate and add to the master dataframe
rep_counter = 1
row = 0
while rep_counter <= reps:
    if rep_counter == data.ix[row,0]:
        for pop in range(pops):
            pop = pop+1
            if pop == 1:
                for line in range(520):
                    generation = data.ix[line+row,1]
                    if generation == 10001:
                        pop1 = population(pop,generation,loci,1000)
                        break
            elif pop == 2:
                for line in range(520):
                    generation = data.ix[line+row,1]
                    if generation == 10001:
                        pop2 = population(pop,generation,loci,1000)
                        break
            elif pop == 3:
                for line in range(520):
                    generation = data.ix[line+row,1]
                    if data.ix[line+row,4] == 1000:
                        pop3 = population(pop,generation,loci,1000)
                        sourcepop(pop1,pop2,pop3,minfreq)
                        t1_s1_weights = []
                        t1_s2_weights = []
                        for allele in range(len(pop3.t1_source.index)):
                            if pop3.t1_source.ix[allele,2] == '[1]':
                                t1_s1_weights.append(pop3.t1_source.ix[allele,0]*pop3.t1_source.ix[allele,1])
                            elif pop3.t1_source.ix[allele,2] == '[2]':
                                t1_s2_weights.append(pop3.t1_source.ix[allele,0]*pop3.t1_source.ix[allele,1])
                        t2_s1_weights = []
                        t2_s2_weights = []
                        for allele in range(len(pop3.t2_source.index)):
                            if pop3.t2_source.ix[allele,2] == '[1]':
                                t2_s1_weights.append(pop3.t2_source.ix[allele,0]*pop3.t2_source.ix[allele,1])
                            elif pop3.t2_source.ix[allele,2] == '[2]':
                                t2_s2_weights.append(pop3.t2_source.ix[allele,0]*pop3.t2_source.ix[allele,1])
                        if t1_s1_weights:
                            t1_pop1_weight_mean.append(np.mean(t1_s1_weights))
                            t1_pop1_weight_sum.append(np.sum(t1_s1_weights))
                        elif not t1_s1_weights:
                            t1_pop1_weight_mean.append('nan')
                            t1_pop1_weight_sum.append('nan')
                        if t1_s2_weights:
                            t1_pop2_weight_mean.append(np.mean(t1_s2_weights))
                            t1_pop2_weight_sum.append(np.sum(t1_s2_weights))
                        elif not t1_s2_weights:
                            t1_pop2_weight_mean.append('nan')
                            t1_pop2_weight_sum.append('nan')
                        if t2_s1_weights:
                            t2_pop1_weight_mean.append(np.mean(t2_s1_weights))
                            t2_pop1_weight_sum.append(np.sum(t2_s1_weights))
                        elif not t2_s1_weights:
                            t2_pop1_weight_mean.append('nan')
                            t2_pop1_weight_sum.append('nan')
                        if t2_s2_weights:
                            t2_pop2_weight_mean.append(np.mean(t2_s2_weights))
                            t2_pop2_weight_sum.append(np.sum(t2_s2_weights))
                        elif not t2_s2_weights:
                            t2_pop2_weight_mean.append('nan')
                            t2_pop2_weight_sum.append('nan')
                        t1_pop1_alleles.append(pop3.t1_source1)
                        t1_pop2_alleles.append(pop3.t1_source2)
                        t1_ambig_alleles.append(pop3.t1_ambiguous)
                        t2_pop1_alleles.append(pop3.t2_source1)
                        t2_pop2_alleles.append(pop3.t2_source2)
                        t2_ambig_alleles.append(pop3.t2_ambiguous)
                        break
                    elif generation == 10500:
                        t1_pop1_alleles.append('nan')
                        t1_pop2_alleles.append('nan')
                        t1_ambig_alleles.append('nan')
                        t2_pop1_alleles.append('nan')
                        t2_pop2_alleles.append('nan')
                        t2_ambig_alleles.append('nan') 
                        break          
        rep_counter = rep_counter+1
    elif rep_counter != data.ix[row,0]:
        while rep_counter != data.ix[row,0]:
            row = row+1

master = pd.concat([master, pd.DataFrame({''.join(['t1_pop1_weight_mean']): t1_pop1_weight_mean})], axis=1)
master = pd.concat([master, pd.DataFrame({''.join(['t1_pop1_weight_sum']): t1_pop1_weight_sum})], axis=1)
master = pd.concat([master, pd.DataFrame({''.join(['t1_pop2_weight_mean']): t1_pop2_weight_mean})], axis=1)
master = pd.concat([master, pd.DataFrame({''.join(['t1_pop2_weight_sum']): t1_pop2_weight_sum})], axis=1)
master = pd.concat([master, pd.DataFrame({''.join(['t2_pop1_weight_mean']): t2_pop1_weight_mean})], axis=1)
master = pd.concat([master, pd.DataFrame({''.join(['t2_pop1_weight_sum']): t2_pop1_weight_sum})], axis=1)
master = pd.concat([master, pd.DataFrame({''.join(['t2_pop2_weight_mean']): t2_pop2_weight_mean})], axis=1)
master = pd.concat([master, pd.DataFrame({''.join(['t2_pop2_weight_sum']): t2_pop2_weight_sum})], axis=1)
master = pd.concat([master, pd.DataFrame({''.join(['t1_pop1_alleles']): t1_pop1_alleles})], axis=1)
master = pd.concat([master, pd.DataFrame({''.join(['t1_pop2_alleles']): t1_pop2_alleles})], axis=1)
master = pd.concat([master, pd.DataFrame({''.join(['t1_ambig_alleles']): t1_ambig_alleles})], axis=1)
master = pd.concat([master, pd.DataFrame({''.join(['t2_pop1_alleles']): t2_pop1_alleles})], axis=1)
master = pd.concat([master, pd.DataFrame({''.join(['t2_pop2_alleles']): t2_pop2_alleles})], axis=1)
master = pd.concat([master, pd.DataFrame({''.join(['t2_ambig_alleles']): t2_ambig_alleles})], axis=1)


# In[ ]:

master.to_csv(str(options.filename)+'_source.csv')

