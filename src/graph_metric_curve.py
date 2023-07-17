import sys
import numpy as np
from matplotlib import pyplot as plt
import csv

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
curves_plotted = 1 ; cmap = plt.get_cmap('cool') ; color_parameter = 2

def plot_curve_data(data,fig=None,ax=None, force_N : int = 0, first_row : bool = False):
    N = data.shape[1]
    V = data.shape[0]
    if (force_N!=0):
        V = force_N
    x_axis = np.linspace(start=1,stop=V,num=data.shape[0],endpoint=True)

    if fig==None:
        fig, ax = plt.subplots()

    global curves_plotted ; global cmap; global color_parameter
    color = lambda : cmap(1 - color_parameter/(curves_plotted+color_parameter-1))
    for i in range(N):
        if i==0:
            if not first_row: continue
        ax.plot(x_axis,sorted(data[:,i]),color=color()) ; curves_plotted+=1
    return fig,ax

fig,ax=(None,None)
files = [f for i,f in enumerate(sys.argv) if i!=0]
V = 0


for i,f in enumerate(files):
    d = get_curve_data(f)
    if (i==0):
        N = get_number_rows(d)
    fig,ax = plot_curve_data(d,fig,ax,V)

plt.show()


