import matplotlib
from matplotlib import pyplot as plt
from itertools import product
import numpy as np
import os
import string
import csv
CMAP = matplotlib.cm.get_cmap('hsv')

def read_csv(csv_f):
    values = []
    with open(csv_f.path) as csvfile:
        reader = csv.reader(csvfile)
        for i,row in enumerate(reader):
            if i==0: continue
            values.append(row[1])
    return np.array(values,dtype=np.float64)

def make_data(directories,targets,condition = (lambda f : True)):
    D = {}
    for d,t in product(dirs,targets):
        if (d.is_dir) and (d.name.startswith(t[0])):
            files_inside = [read_csv(f) for f in os.scandir(d.path) if f.name.endswith(".csv") and condition(f)]
            D[t[1]] = files_inside
    return D

dirs = os.scandir("Experimentos/CPPerigee-and-experiments/41_different_term_rounds/")

tar = [
    ('4',"primero","red"),
    ('5',"segundo","green"),
    ('6',"tercero","blue"),
]

cond = lambda f: ("percentile" in f.name) and (not ("bound" in f.name))
datos = make_data(dirs,tar,cond)

fig,ax = plt.subplots(figsize=(10,10))

PLOT_x_vals = []
for t in tar:
    XX = datos[t[1]]
    for xx in XX:
        ax.hist(xx,20,alpha=0.2,color=t[2])
        PLOT_x_vals.append((np.mean(xx),t[2]))

def regularize_plot(ax):
    global PLOT_x_vals
    y_lim = ax.get_ylim()[1]
    for x,c in PLOT_x_vals:
        ax.vlines(x,0.,y_lim * 1.25,color=c)
    ax.set_ylim(0.,y_lim*1.5*1.618)

regularize_plot(ax)
    

ax.legend()
plt.show()