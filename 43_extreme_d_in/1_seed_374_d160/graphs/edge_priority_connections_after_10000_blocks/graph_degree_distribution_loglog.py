import sys
import numpy as np
from matplotlib import pyplot as plt
import csv

def plot_loghist(x, nbins,cutoff=8.):
  dexp = pow(np.max(x),1./nbins)
  N = int(np.ceil(np.log(np.max(x))/np.log(dexp)))
  bins = [ pow(dexp,k) for k in range(0,N+1) ]
  hist, bins = np.histogram(x, bins=bins)
  pairs = [(x,y) for y,x in zip(hist,bins) if y!=0 and x>cutoff]
  xx,yy = zip(*pairs)
  plt.plot(xx,yy)
  plt.xscale("log")
  plt.yscale("log")
  return yy,xx

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

def get_degree_data(data):
    d_out = []
    d_in  = []
    # each row will have the number of the vertex
    # plus several triples of numbers (v,luv,lvu)
    # after that, infinitiy and -1
    # then d_out(v) = (number of non infinite entries - 1)/3 
    defined = lambda x : (x!="Inf") and (x!="-1")

    for i,row in enumerate(data):
        d = (len([val for val in row if defined(val)])-1)/3 
        if i%2==0: # outgoing
            d_out.append(d)
        else:      # incoming
            d_in.append(d)
    return d_out, d_in

data = get_data(sys.argv[1])
d_out,d_in = get_degree_data(data)

nbins =  20
bins = plot_loghist(d_in,nbins)

for d,U in zip(bins[1],bins[0]): print("d: ",d," -> ",U)

plt.show()