import sys
import numpy as np
from matplotlib import pyplot as plt
import csv

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

data = get_data(sys.argv[1])
xx,yy = convert_to_cdf(data)
plt.plot(xx,yy)
plt.show()