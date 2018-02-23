import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
import os
import pprint
import time
import copy
pp = pprint.PrettyPrinter()
#matplotlib.rcParams['figure.figsize'] = [11.0, 8.0]


# Uncomment to enable matplotlib-D3.js integration -- enables zooming and panning of figures.
#import mpld3
#mpld3.enable_notebook()

def plotter(instance, method, gamma, gc_vals, main_logdir, max_x = None, symlog=True,savefig=False):
    ### Get the data
    cwd = os.getcwd()
    name = instance+"_"+method+"_"+gamma+"_"
    logdir = main_logdir+"/"+gamma+"/"+method+"/"
    def getdata(gc_vals):
        os.chdir(logdir)
        d = defaultdict(list)
        for gc in gc_vals:
            data = open(name+str(gc)+".log", "r")
            for line in data.readlines():
                vals = line.split(" ")
                while "\n" in vals:
                    vals.remove("\n")
                vals = [int(v) for v in vals]
                #vals.sort()
                if vals:
                    d[gc].append(vals)
            data.close()
        os.chdir(cwd)
        return d
    
    start=time.time()
    
    d = getdata(gc_vals)
    d_keys = d.keys()
    d_keys.sort()
    
    dd = copy.deepcopy(d)
    ddd = copy.deepcopy(d)
    
    # Unpack each d[gc]
    for gc in d_keys:
        acc, acc2 = [], []
        for i in d[gc]:
            acc += i
            acc2 += [min(i)]*len(i)
        dd[gc] = acc
        ddd[gc] = acc2
    
    
    if max_x is None:
        max_iter = max(len(dd[gc]) for gc in dd)
    else:
        max_iter = max_x
    
    max_obj = max(max(dd[gc]) for gc in dd)
    min_obj = min(min(dd[gc]) for gc in dd)
    xs = np.arange(0, max_iter)
    ys = []
    for gc in d_keys:
        y = [[], []]
        for i in range(max_iter):
            if i < len(dd[gc]):
                y[0].append(min(ddd[gc][:i+1]))
                y[1].append(dd[gc][i])
            else:
                y[0].append(min(dd[gc]))
                y[1].append(min(dd[gc]))
        ys.append(y)
    
    cmap = plt.get_cmap('jet')


    for (gc,y) in zip(d_keys,ys):
        color = cmap((d_keys.index(gc))/float(len(d_keys)))
        # adding minus one because we always a) write initial obj, b) final obj
        total_iters = len(dd[gc]) - 1
        plt.fill_between(xs, y[0], y[1], color=color, alpha=0.3, label="g1="+str(gc)+" AugILP: " + str(total_iters))
        plt.plot(xs, y[0], color=color, alpha=0.7)

    plt.xlabel('iteration')
    plt.ylabel('objective')
    plt.ylim(min_obj-1,max_obj+100)
    if symlog:
        plt.gca().set_xscale('log')


    #plt.title("instance: " + instance + " method: " + method + " gamma: " + gamma)

    plt.legend()
    fig = plt.gcf()
    if savefig:
        fig.savefig(name+"clip_"+str(max_x)+"_inner.pdf", format='pdf', dpi=1000,bbox_inches="tight")
    plt.show()
