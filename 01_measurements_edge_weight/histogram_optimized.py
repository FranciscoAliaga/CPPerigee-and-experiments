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
        res = int(re.findall(r'\d+',filename)[0])
    except: pass
    return res

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
    files    = [f.name for f in os.scandir(dir) if os.path.isfile(f)]
    csvfiles = [os.path.join(dir,f) for f in files if f.endswith('.csv')]
    return csvfiles

# Chosen Colormap
cmap = plt.get_cmap('hsv')

dirs = [f for i,f in enumerate(sys.argv) if i!=0]

# percentiles #########################################################
fig,ax = plt.subplots(figsize=(6,4))

alpha = 0.4

for n,d in enumerate(dirs):
    archivo = get_csvfiles(d)[0]
    data = get_curve_data(archivo)
    xx = data[:,1]

    desv_estandar = np.std(xx)
    label = "$\sigma = {s:.1f}$ [s]".format(s=desv_estandar)
    color = cmap(n/len(dirs))

    ax.hist(xx,bins=30,alpha=alpha,label=label,color=color)

#_ax.yaxis.set_major_formatter(mtick.PercentFormatter())
ax.set_xlabel("Latencia [s]")
ax.set_ylabel("Nodos")
ax.set_title("Latencia global para red bajo protocolo Edge Priority.")
ax.grid()
ax.legend(loc="upper right",fontsize=9)

#
plt.show()


