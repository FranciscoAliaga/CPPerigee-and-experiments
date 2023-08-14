import os
import sys
import re
import numpy as np
from matplotlib import pyplot as plt
import csv

import matplotlib.ticker as mtick
def file_number(filename):
    res = None
    try:
        res = int(re.findall(r'\d+',filename)[-1])
    except: pass
    return res

clamp = lambda x : max(min(x,1.),0.)

def get_curve_data(filename):
    data = []
    with open(filename, newline='') as csvfile:
        reader = csv.reader(csvfile)
        for i,row in enumerate(reader):
            if i==0: continue
            data.append([float(r) for r in  row])
    data = np.array(data)
    return data


def get_csvfiles(dir='.'):
    files    = [os.path.join(dir,f.name) for f in os.scandir(dir) if os.path.isfile(f)]
    csvfiles = [f for f in files if f.endswith('.csv')]
    return csvfiles

def get_curves_in_order(dir = ''):
    files = get_csvfiles(dir)

    edge_priority_percentile_curves = sorted([
        (file_number(f),f) for f in files if "edge_priority_percentile" in f
    ])

    edge_priority_edge_weight_curves = sorted([
        (file_number(f),f) for f in files if "edge_priority_edge" in f
    ])

    subset_percentile_curves = sorted([
        (file_number(f),f) for f in files if "subset_percentile" in f
    ])

    subset_edge_weight_curves = sorted([
        (file_number(f),f) for f in files if "subset_edge" in f
    ])

    _,edge_priority_percentile_curves = zip(*edge_priority_percentile_curves)
    _,edge_priority_edge_weight_curves = zip(*edge_priority_edge_weight_curves)
    _,subset_percentile_curves = zip(*subset_percentile_curves)
    _,subset_edge_weight_curves = zip(*subset_edge_weight_curves)

    e0 = [f for f in files if "edge_weight_start.csv" in f][0]
    p0 = [f for f in files if "percentile_90_start.csv" in f][0]

    edge_priority_edge_weight = [e0]
    edge_priority_percentiles = [p0]
    subset_edge_weight = [e0]
    subset_percentiles = [p0]

    edge_priority_edge_weight.extend(edge_priority_edge_weight_curves)
    edge_priority_percentiles.extend(edge_priority_percentile_curves)
    subset_edge_weight.extend(subset_edge_weight_curves)
    subset_percentiles.extend(subset_percentile_curves)

    res = {
        "edgepriority edgeweight":  edge_priority_edge_weight  ,
        "subset edgeweight":        subset_edge_weight  ,
        "edgepriority percentiles": edge_priority_percentiles  ,
        "subset percentiles":       subset_percentiles,
    }
    return res

# Chosen Colormap
cmap = plt.get_cmap('jet')
dirs = [d for x,d in enumerate(sys.argv) if x!=0]

alpha = clamp(0.3 + 1./len(dirs))

# percentiles #########################################################
fig,ax = plt.subplots(nrows=2,ncols=1,figsize=(9,8))

for n,d in enumerate(dirs):
    F = get_curves_in_order(d)

    subset_percentiles = F["subset percentiles"]
    edge_priority_percentiles = F["edgepriority percentiles"]

    #subset
    for i,f in enumerate(subset_percentiles):
        color = cmap(i / len(subset_percentiles))
        data = get_curve_data(f)
        xx = sorted(data[:,1])
        yy = np.linspace(start=1./len(xx),stop=100.,num=len(xx),endpoint=True)

        filenumber = file_number(f)
        if filenumber==90: filenumber=0
        label = "$|B|$ = {b:<4,}".format(b=filenumber)

        # only plot the labels for the first ones
        if n==0:
            ax[0].plot(xx,yy,color=color,label=label,alpha=alpha)
        else:
            ax[0].plot(xx,yy,color=color,alpha=alpha)

    for i,f in enumerate(edge_priority_percentiles):
        color = cmap(i / len(edge_priority_percentiles))
        data = get_curve_data(f)
        xx = sorted(data[:,1])
        yy = np.linspace(start=1./len(xx),stop=100.,num=len(xx),endpoint=True)

        filenumber = file_number(f)
        if filenumber==90: filenumber=0
        label = "$|B|$ = {b:<4,}".format(b=filenumber)

        # only plot the labels for the first ones
        if n==0:
            ax[1].plot(xx,yy,color=color,label=label,alpha=alpha)
        else:
            ax[1].plot(xx,yy,color=color,alpha=alpha)

for _ax in ax:
    _ax.yaxis.set_major_formatter(mtick.PercentFormatter())
    _ax.set_xlabel("Latencia [ms]")
    _ax.set_ylabel("Nodos")
    _ax.grid()
    _ax.legend(loc="lower right",fontsize=7)

ax[0].set_title("Curva $\lambda^{90\%}$ para protocolo Subset")
ax[1].set_title("Curva $\lambda^{90\%}$ para protocolo Edge Priority")

# edge weights ########################################################

fig2,ax2 = plt.subplots(figsize=(7,4))
#
for n,d in enumerate(dirs):
    F = get_curves_in_order(d)
    subset_edge_weights = F["subset edgeweight"]
    edgepriority_edge_weights = F["edgepriority edgeweight"]

    f_clean = lambda f : file_number(f) if file_number(f)%1000 == 0 else 0

    ##subset
    bb,ww = zip(*[(f_clean(f),np.mean(get_curve_data(f)[:,1])/1000.) for f in subset_edge_weights])
    # label only on first n
    if n==0:
        ax2.plot(bb,ww,color = "blue",label="Subset",marker="*",alpha=alpha)
    else:
        ax2.plot(bb,ww,color = "blue",marker="*",alpha=alpha)
    #
    bb,ww = zip(*[(f_clean(f),np.mean(get_curve_data(f)[:,1])/1000.) for f in edgepriority_edge_weights])
    if n==0:
        ax2.plot(bb,ww,color="red",label="EdgePriority",marker="*",alpha=alpha)
    else:
        ax2.plot(bb,ww,color="red",marker="*",alpha=alpha)
#
ax2.set_title("Métrica de latencia global v/s número de bloques")
ax2.set_xlabel("Número de bloques |B|")
ax2.set_ylabel("Latencia [s]")
ax2.legend()
ax2.grid()


#
plt.show()


