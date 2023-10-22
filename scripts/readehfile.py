#!/usr/bin/python3
# Pavel Loskot, January 2022
# functions to process event history files

from time import time
import os.path
from sys import argv
from collections import Counter
import matplotlib.pyplot as plt

# simulation folders
def get_simfolder(a):
    sims = {
        '1':'./BCR2/BCR_eh.gdat',
        '2':'./EGFR2/egfr_net_eh.gdat',
        '3':'./FceRI_gamma2/fceri_gamma2_eh.gdat',
        '4':'./FceRI_ji2/FceRI_ji2_eh.gdat',
        '5':'./Multi-site2/multisite2_eh.gdat',
        '6':'./Multi-state2/multistate_eh.gdat',
        '7':'./chemotaxisAdaptation2/AN_eh.gdat',
        '8':'./chemotaxisAdaptation2x/ANx_eh.gdat',
        '9':'./tlbr2/tlbr_eh.gdat'
    }
    if not a:    
        for k,v in sims.items():
            print('{}: {}'.format(k,v))
        a = input('simulation? ')
    else:
        a = a[0]
    return sims[str(a)]

# obtain model name from filename
def get_modelname(f):
    models = {
        'BCR2':'BCR',
        'EGFR2':'EGFR',
        'FceRI_gamma2':'FceRI-gamma',
        'FceRI_ji2':'FceRI-ji',
        'Multi-site2':'Multi-sites',
        'Multi-state2':'Multi-states',
        'chemotaxisAdaptation2':'Chemotaxis',
        'chemotaxisAdaptation2x':'Chemotaxis-ext',
        'tlbr2':'TLBR'
    }
    return models[f.rsplit("/",2)[1]]

# simpler exit
class _Exit:
    def __repr__(self):
        raise SystemExit
exit= _Exit()

# print global parameters
def print_params(G, sep=None):
    if sep is not None: print(sep)
    for p in dir(G):
        if p.startswith("_"): continue
        if isinstance(G.__dict__[p],int|float) and G.__dict__[p]>1e3:
            print(f'{p} = {G.__dict__[p]:g}')
        else:
            print(f'{p} = {G.__dict__[p]}')
    if sep is not None: print(sep)            

# elapsed time
def elapsed_time(start):
    t = time()-start
    print('Elapsed time: {:.0f} min {:.2f} sec'.format(t//60,t%60))

# export plot into file
def expfig(fname):
    plt.savefig(fname, bbox_inches='tight', dpi=100)
    
# read csv file line by line skipping comments
def getcsvlines(Params, lastrow=None):
    if not os.path.isfile(Params.fin):
        print('File {s} does not exist.'.format(Params.fin))
        exit()
    with open(Params.fin,'r') as f:
        if Params.verbose:
            print(getcsvlines.__name__,': Reading file',Params.fin)
        n = 0
        if lastrow is not None:
            for line in f:
                if line.startswith(lastrow): break                
                if line.startswith('#'): continue
                n += 1
                if n>Params.maxlines: break
                yield line.strip().split(',')
        else:
            for line in f:
                if line.startswith('#'): continue
                n += 1
                if n>Params.maxlines: break
                yield line.strip().split(',')
        if Params.verbose:
            print(getcsvlines.__name__,': Finished reading file',Params.fin)
        return None

# count the total number of reactions in a file
def nrxninfile(Params):
    if Params.verbose:
        print(nrxninfile.__name__,': Counting records in',Params.fin)
    if Params.block:
        Nrxn = {}
        for l in getcsvlines(Params, lastrow='#end'):
            Nrxn[int(l[0])] = Nrxn.get(int(l[0]),0) + 1
    else:
        Nrxn = 0
        for l in getcsvlines(Params, lastrow='#end'):
            Nrxn += 1
    return Nrxn

# count the total number of reactions in a list of tuples
def nrxninlist(tt, Params):
    if Params.verbose:
        print(nrxninlist.__name__,': Counting records in a tuple list')
    if Params.block:
        Nrxn = {}
        for l in list(tt):
            Nrxn[int(l[0])] = Nrxn.get(int(l[0]),0) + 1
    else:
        Nrxn = 0
        for l in list(tt):
            Nrxn += 1
    return Nrxn

# create dictionary of reactions {reaction ID : name}
def get_reactions(Params):
    if Params.verbose:
        print(get_reactions.__name__,': Identifying unique reactions')
    rxnSet = set()
    for l in getcsvlines(Params,lastrow='#end'): 
        rxnSet.add((l[5],int(l[4])))
    return rxnSet

# obtain tuples in each block
def get_tuples(Params):
    if Params.verbose: print(get_tuples.__name__,': Extracting tuples')
    b = 0 # block index
    t = [] 
    for l in getcsvlines(Params, lastrow='#end'):
        if not l: return
        if Params.block and int(l[0])>b: # new block arrived
            b += 1
            if b>Params.maxblocks: return            
            t = []
        t.append(int(l[4]))
        if len(t)==Params.N:
            yield (b,tuple(t)) if Params.block else tuple(t)
            t.pop(0)

# obtain tuples in each block with given size
def get_tuples1(fname, N, blocksize, maxlines=1e9, maxblocks=1e9, verbose=False):
    if verbose: print(get_tuples1.__name__,
                          ': Extracting tuples in blocks of size', blocksize)
    b = 0 # block index
    c = 0 # in-block tuple counter
    t = []
    for l in getcsvlines(fname,lastrow='#end', maxlines=maxlines, verbose=verbose):
        if not l: break
        if c>=blocksize: # start new block
            c = 0
            b += 1
            if b>maxblocks: return
            t = []
        t.append(int(l[4]))
        if len(t)==N:
            c += 1                        
            yield (b,tuple(t))
            t.pop(0)
            
# extract time differences in each block
def get_dtimes(Params):
    if Params.verbose: print(get_times.__name__,': Extracting time instances')
    if Params.block:
        b = 0 # block index
        dtimes = {}
    else:
        dtimes = []
    prevtime = None
    for l in getcsvlines(Params, lastrow='#end'):
        if not l: break
        if not prevtime:
            prevtime = float(l[1])
            continue
        if Params.block and int(l[0])>b: # new block arrived
            b += 1
            if b>Params.maxblocks: return
            dtimes[b] = []
        if Params.block:
            dtimes[b].append(float(l[1])-prevtime)
        else:
            dtimes.append(float(l[1])-prevtime)
        prevtime = float(l[1])
    return dtimes
            
# obtain histogram for clusters
def get_hist(tt, Params):
    if Params.verbose:
        print(get_hist.__name__,': Calculating histogram of tuples')
    return dict(Counter(e for e in tt))

# clustering of histogram tuples
def cluster_hist(hin, Params):
    if Params.verbose: print(cluster_hist.__name__,
                      ': Clustering tuples based on their frequency counts') 
    hout = {}
    c = -1
    v_prev = 1e10
    for k in [k for k in sorted(hin, key=hin.get, reverse=True)]:
        if Params.thr < (v_prev-hin[k])*100/v_prev:
                c += 1
                hout[c] = {}
        v_prev = hin[k]
        hout[c][k] = hin[k]
    return hout

# evaluate cluster statistics
def clust_stats_print(hc, Params):
    if Params.verbose:
        print(clust_stats.__name__,': Calculating cluster statistics')
    s1 = s2 = 0
    print('Tuples:  Unique   Frq-min Frq-max  Total')
    print('cluster---------------------------------')    
    for c in hc:
        g = [e for e in hc[c].values()]
        print("{:7} {:7} {:7} {:7} {:8}".format(c,len(g),min(g),max(g),sum(g)))
        s1 += len(g)
        s2 += sum(g)
    print('----------------------------------------')            
    print("tot:{:3} {:7} {:24}".format(c+1,s1,s2))
    return s2

# evaluate cluster statistics: summary
def clust_stats_print1(hc, Params):
    if Params.verbose:
        print(clust_stats.__name__,': Calculating cluster statistics')
    stat = {}
    for c in hc:
        g = [e for e in hc[c].values()]
        stat[c] = (len(g),min(g),max(g),sum(g))
    return stat

# find shared reactions
def clust_unions(hc, verbose=False):
    if verbose: print(clust_unions.__name,': Calculating common reactions')
    rxn = {}
    for c in hc:
        s = []
        for e in hc[c]: s.append(set(e))
        rxn[c] = set.union(*s)
    return rxn

        
