#!/usr/bin/python3 -i
# Pavel Loskot, January 2022
# find unique reactions and tuples and their counts in whole file

from readehfile import *

start_time = time()

# global parameters
class Params:
    N = 3
    block = True # do not change
    verbose = False
    maxlines = 1e9
    maxblocks = 100
    fin = get_simfolder(argv[1:])
print_params(Params, sep='--')
    
# get all reactions
#print(' ... identifying unique reactions')
RxnSet = get_reactions(Params)
# for r in RxnSet: print(r[0],':',r[1])
print('number of reaction types:',len(RxnSet))

# count tuples and reactions
#print(' ... counting tuples')
tt = get_tuples(Params)
print('There are {} {}-tuples'.format(len(list(tt)),Params.N))
#print()

#print(' ... counting reactions in',Params.fin)
n = nrxninfile(Params)
print('number of reactions total:',n)
tt = get_tuples(Params)
rxn = list(nrxninlist(tt,Params).values())
#print(' ... and per block:')
#print(rxn)
#print()

print('number of blocks:',len(rxn))
print('reactions per block:',
          min(rxn),'-',max(rxn),',',round(sum(rxn)/len(rxn)))

# elapsed time    
elapsed_time(start_time)

