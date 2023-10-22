#!/usr/bin/python3 
# Pavel Loskot, January 2022
# cluster tuples

from readehfile import *

start_time = time()

# global parameters
class Params:
    N = int(argv[2])
    block = True # do not change
    verbose = False
    maxlines = 1e9
    maxblocks = 100
    thr = 30    
    fin = get_simfolder(argv[1:])
print_params(Params, sep='--')

# create and statistically evaluate tuples
tt = get_tuples(Params)
if Params.block:
    # separate blocks    
    blocks = set((e[0] for e in tt))
    ttb = {}
    for b in blocks: ttb[b] = []
    tt = get_tuples(Params)
    for t in tt: ttb[t[0]].append(t[1])
    histc = {}
    clusters = set()
    d1,d2,d3,d4= {},{},{},{}
    a0 = (1e9,1e9)
    a1 = (0, 0)
    for b in blocks:
        print("\rBlock",b,end='')
        # histogram and clustering
        histc[b] = cluster_hist(get_hist(ttb[b], Params),Params)
        # show clusters info
        #stat = clust_stats_print(histc[b], Params)
        #continue
        stat = clust_stats_print1(histc[b], Params)
        clusters = clusters.union(set(stat.keys()))
        #print()
        for c in clusters:
            s0 = stat.get(c, (1e9,1e9,1e9,1e9))
            s1 = stat.get(c, (-1,-1,-1,-1))
            d1[c] = (min(d1.get(c,a0)[0], s0[0]), max(d1.get(c,a1)[1], s1[0]))
            d2[c] = (min(d2.get(c,a0)[0], s0[1]), max(d2.get(c,a1)[1], s1[1]))
            d3[c] = (min(d3.get(c,a0)[0], s0[2]), max(d3.get(c,a1)[1], s1[2]))
            d4[c] = (min(d4.get(c,a0)[0], s0[3]), max(d4.get(c,a1)[1], s1[3]))
    print()
    for c in clusters:
        print(c,'&',Params.N,'&','{}--{} & {}--{} & {}--{} & {}--{} \\\\'.format(
            *d1[c],*d2[c],*d3[c],*d4[c]))
    #print("\multicolumn{4}{c}{%d}" % N)
    #print(*list(clusters),sep=" & ")
    #for c in clusters: print("{}--{} &".format(d1[c][0],d1[c][1]),end='')
    #print()
    #for c in clusters: print("{}--{} &".format(d2[c][0],d2[c][1]),end='')
    #print()
    #for c in clusters: print("{}--{} &".format(d3[c][0],d3[c][1]),end='')
    #print()
    #for c in clusters: print("{}--{} &".format(d4[c][0],d4[c][1]),end='')
    #print()    
    print()
    
else:
    # no blocks
    tt = get_tuples(Params)
    histc = {}
    stot = 0
    # histogram and clustering
    histc = cluster_hist(get_hist(tt,Params),Params)
    # show clusters info
    clust_stats_print(histc, Params)

# elapsed time    
elapsed_time(start_time)
        

