import sys
import numpy as np
from matplotlib import pyplot as plt
import csv

import matplotlib.ticker as mtick

def get_curve_data(filename):
    data = []
    with open(filename, newline='') as csvfile:
        reader = csv.reader(csvfile)
        for i,row in enumerate(reader):
            if i==0: continue
            data.append([float(r) for r in  row])
    data = np.array(data)
    return data

def get_number_rows(data):
    return data.shape[0]


# color ramp
curves_plotted = 1 ; cmap = plt.get_cmap('hsv') ; color_parameter = 3

def plot_curve_data(data,fig=None,ax=None, force_N : int = 0, first_row : bool = False):
    N = data.shape[1]
    V = data.shape[0]
    if (force_N!=0):
        V = force_N
    num = data.shape[0]
    y_axis = np.linspace(start=1.0,stop=100.,num=num,endpoint=True)
    if fig==None:
        fig, ax = plt.subplots()

    global curves_plotted ; global cmap; global color_parameter
    #color = lambda : cmap(curves_plotted/color_parameter)
    color = lambda x : "darkblue" if (x>950) else "red"
    for i in range(N):
        if i==0:
            if not first_row: continue
        xx = sorted(1000.*data[:,i])
        
        ax.plot(xx,y_axis,color=color(np.min(xx)),alpha = 0.2)
        # plot each tenth percentile
        py = [10*i for i in range(0,11)]
        px = [np.percentile(xx,i) for i in py]
        #ax.scatter(px,py,color='gray',marker="|",alpha=1.)
        # plot each average over percentile
        avepx = [np.average([x for x in xx if low < x and x < high]) for low,high in zip(px, px[1:]) ]
        avepy = [ (low+high)*0.5 for low,high in zip(py, py[1:]) ]
        ax.scatter(avepx,avepy,color=color(np.min(xx)),alpha=0.1,s=100.,marker="|")

        # plot hlines
        min_x = np.min(xx); max_x = np.max(xx)
        for y in py:
            ax.hlines(y,min_x,max_x,alpha=0.05,color="black")
    #curves_plotted+=1
    return fig,ax

fig,ax = plt.subplots(figsize=(7,3))
#fig,ax = (None,None)
files = [f for i,f in enumerate(sys.argv) if i!=0]
V = 0

for i,f in enumerate(files):
    d = get_curve_data(f)
    if (i==0):
        N = get_number_rows(d)
    fig,ax = plot_curve_data(d,fig,ax,V)


ax.set_title("Curvas $\lambda^{90\%}$ para redes aleatorias (azul) y para redes entrenadas \n con $|B|=4000$ bloques bajo protocolo Edge Priority (rojo).")
ax.set_xlabel("Latencia [ms]")
ax.set_ylabel("Nodos")
ax.yaxis.set_major_formatter(mtick.PercentFormatter())





plt.show()


