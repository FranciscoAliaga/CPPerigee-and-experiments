import matplotlib
from matplotlib import pyplot as plt
from itertools import product
import numpy as np
import os
import string
import csv

dirs = os.scandir("Experimentos/3_scaling_on_size_fixed_parameters")
#CMAP = matplotlib.cm.get_cmap('Pastel1')
CMAP = matplotlib.cm.get_cmap('hsv')

PLOT_x_vals = []
PLOT_N = 1

PLOT_COLORS = iter(CMAP(i/10.) for i in range(20))
def plot_curve(experiment_dir,ax):
    global PLOT_x_vals
    global PLOT_N
    global PLOT_COLORS
    ## get .csv file
    files = os.scandir(experiment_dir.path)
    csv_f = None
    for csv_f in (f for f in files if f.name.endswith('.csv')):
        #try:
        #    csv_f = next((f for f in files if f.name.endswith('.csv')))
        #except Exception:
        #    pass
        #if(csv_f==None):
        #    return
        ## open the csv file and retrieve second column
        values = []
        with open(csv_f.path) as csvfile:
            reader = csv.reader(csvfile)
            for i,row in enumerate(reader):
                if i==0: continue
                values.append(row[1])
        values = np.array(sorted(values),dtype=float)
    
        C =  "red" #next(PLOT_COLORS)
        ax.hist(values,20,density=True,label = "" , alpha = 0.8,color = C)
        #ax.hist(values,30,density=True,label = "experimento {n}".format(n=PLOT_N) , alpha = 0.8,color = C)
        #if(np.mean(values)<3300): return
        #yy = values
        #xx = np.linspace(0.,1.,num=len(yy))
        #ax.plot(yy,xx,label="experimento{n}".format(n=PLOT_N))
        PLOT_x_vals.append((np.mean(values),C))
        PLOT_N +=1

def regularize_plot(ax):
    global PLOT_x_vals
    y_lim = ax.get_ylim()[1]
    for x,c in PLOT_x_vals:
        ax.vlines(x,0.,y_lim * 1.25,color=c)
    ax.set_ylim(0.,y_lim*1.5*1.618)

#seeds = [
#    '0',
#    '21574',
#    '32452',
#    '74567',
#    '12312'
#]

seeds = [
    '4',
]

experimentos = {s:[] for s in seeds}
experimentos_totales = []

for s,d in product(seeds,dirs):
    #if (d.name.endswith(s)):
    if (d.name.startswith(s)):
        experimentos[s].append(d)
        experimentos_totales.append(d)

_,experimentos_totales = zip(*sorted(zip([x.name for x in experimentos_totales],experimentos_totales)))

fig,ax = plt.subplots(figsize=(10,10))

for exp in experimentos_totales:
    plot_curve(exp,ax)

regularize_plot(ax)
ax.legend()


plt.show()