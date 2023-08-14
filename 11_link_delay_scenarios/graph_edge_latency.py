import os
import re
import sys
import numpy as np
from matplotlib import pyplot as plt
import csv
import matplotlib.ticker as mtick

def get_data(filename):
    data = []
    with open(filename, newline='') as csvfile:
        reader = csv.reader(csvfile)
        for i,row in enumerate(reader):
            data.append(row)
    return data

def get_parameters(data):
    V = len(data)/2
    return V

def get_edge_data(data):
    # each row will have the number of the vertex
    # plus several triples of numbers (v,luv,lvu)
    # after that, infinitiy and -1
    defined = lambda x : (x!="Inf") and (x!="-1")
    l_out = []
    l_in  = []

    for i,row in enumerate(data):
        out = (i%2==0)
        new_l = []
        N = int((len(row)-1)/3)
        for x in range(N):
            dx = 1 if out else 2
            luv = row[1+3*x+dx]
            if not defined(luv):
                break
            new_l.append(float(luv))
        if out: l_out.append(new_l)
        else: l_in.append(new_l)
    return l_out,l_in

def convert_to_cdf(data):
    l_out,l_in = get_edge_data(data)
    L = [latency for outgoing in l_out for latency in outgoing]
    xx = sorted(L)
    yy = np.linspace(start=1./len(xx),stop=100.,num=len(xx),endpoint=True)
    return xx,yy

def num(f):
    res = 0
    try:
        res = int(re.findall(r'\d+',f)[-1])
    except: pass
    return res

clamp = lambda x : min(1,max(0,x))

def get_each_graph_file_in_order(exp):
    #get the graphs folder
    graphs_folder = [f for f in os.scandir(exp) if "graphs" in f.name][0]
    # search within
    random_graph_file = os.path.join(exp,"graphs/randomized_connections_startup/graph.csv")
    subset_graphs       = [(0,random_graph_file)]
    edgepriority_graphs = [(0,random_graph_file)]

    join = lambda f : os.path.join(os.path.join(exp,"graphs"),os.path.join(f.name,"graph.csv"))

    subset_graphs.extend(sorted([
        (num(f.name),join(f)) for f in os.scandir(graphs_folder) if "subset_connections" in f.name
    ]))

    edgepriority_graphs.extend(sorted([
        (num(f.name),join(f)) for f in os.scandir(graphs_folder) if "edge_priority" in f.name
    ]))


    return {
        "subset" : subset_graphs,
        "edgepriority" : edgepriority_graphs
    }


### 

fig,ax = plt.subplots(nrows=2,figsize=(7,7))
cmap = plt.get_cmap('jet')
experiments = [f for (i,f) in enumerate(sys.argv) if i!=0]
alpha = clamp(0.3 + 1./len(experiments))

for n,experiment in enumerate(experiments):
    G = get_each_graph_file_in_order(experiment)

    subset = G["subset"]
    edgepriority = G["edgepriority"]

    for i,(blocks,g) in enumerate(subset):
        color = cmap(i/len(subset))
        data = get_data(g)
        xx,yy = convert_to_cdf(data)
        label = "$|B|$ = {b:<4,}".format(b=blocks)

        if n==0:
            ax[0].plot(xx,yy,color=color,alpha=alpha,label=label)
        else:
            ax[0].plot(xx,yy,color=color,alpha=alpha)

    for i,(blocks,g) in enumerate(edgepriority):
        color = cmap(i/len(edgepriority))
        data = get_data(g)
        xx,yy = convert_to_cdf(data)
        label = "$|B|$ = {b:<4,}".format(b=blocks)

        if n==0:
            ax[1].plot(xx,yy,color=color,alpha=alpha,label=label)
        else:
            ax[1].plot(xx,yy,color=color,alpha=alpha)


for _ax in ax:
    _ax.yaxis.set_major_formatter(mtick.PercentFormatter())
    _ax.grid()
    _ax.set_xlabel("Latencia par a par [ms]")
    _ax.set_ylabel("Conexiones")
    _ax.legend(loc="lower right",fontsize=6)

ax[0].set_title("CDF de latencias de conexiones para protocolo Subset\n para experimento con valor de latencia m=10")
ax[1].set_title("CDF de latencias de conexiones para protocolo EdgePriority")


plt.show()