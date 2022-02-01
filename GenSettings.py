
# coding: utf-8

# In[ ]:


import argparse


# In[ ]:


parser = argparse.ArgumentParser(description='Generate Settings File')

parser.add_argument("-n", dest="filename", 
                 help="name of simulation")

parser.add_argument("-t", dest="template", 
                 help="name of template file")

parser.add_argument("-le", dest="epiloci", action="store",
                 help="specify the number of epistatic loci")

parser.add_argument("-l", dest="loci", action="store",
                 help="specify the number of loci")

parser.add_argument("-u", dest="mutation_rate", action="store", default=0.0001,
                 help="specify mutation rate")

parser.add_argument("-ue", dest="mutation_rate_e", action="store", default=0.0001,
                 help="specify epistatic mutation rate")

parser.add_argument("-r", dest="replicates", action="store", default=1,
                 help="specify replicates of the simulation")

parser.add_argument("-b", dest="burn", action="store", default=10000,
                 help="specify replicates of the simulation")

parser.add_argument("-d", dest="dispersal", action="store", default=500,
                 help="specify replicates of the simulation")

parser.add_argument("-s", dest="selection", action="store", default=10,
                 help="specify selection intensity of the epistasis trait")

parser.add_argument("-q", action="store_false", dest="verbose", default=True,
                 help="don't print status messages to stdout")


args = parser.parse_args()


# In[ ]:


with open(str(args.template+".ini"), 'r') as f:
    filedata = f.read()
    f.close()


# In[ ]:


## calculate mutational variance
effect = 0.5/int(args.loci)

effecte=.5/int(args.epiloci)

## calculate dispersal probabilities


args.dispersal
# In[95]:

## replace template ini settings with custom settings

filedata = filedata.replace('SELECTION', str(args.selection))

filedata = filedata.replace('GENERATIONS', str(args.burn + args.dispersal))
filedata = filedata.replace('N', str(args.filename))

filedata = filedata.replace('REPS', str(args.replicates))
filedata = filedata.replace('LOCI', str(args.loci))
filedata = filedata.replace('EL', str(args.epiloci))

filedata = filedata.replace('EFFECT', str(effect))
filedata = filedata.replace('DISP', str(args.burn + 1))

filedata = filedata.replace('MU', str(args.mutation_rate))



# In[ ]:


filename = str(args.filename)+'.ini'
with open(filename, 'w') as w:
    w.write(filedata)
    w.close()

