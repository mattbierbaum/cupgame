import numpy as np
import scipy as sp
import pylab as pl

cuplist = []

"""
For example:
./cupgame-track 1 0.451171875 0.392578125 2 huh
./cupgame-track 1 0.451171875 0.393578125 2 huh1
./cupgame-track 1 0.452171875 0.393578125 2 huh2
./cupgame-track 1 0.452171875 0.394578125 2 huh3
./cupgame-track 1 0.452371875 0.394578125 2 huh4
./cupgame-track 1 0.452371875 0.391578125 2 huh5
./cupgame-track 1 0.452471875 0.391578125 2 huh6
./cupgame-track 1 0.452671875 0.391578125 2 huh7
./cupgame-track 1 0.456971875 0.390978125 2 huh8
"""

def cupnear(x, y):
    tx = np.round(y / np.sqrt(3))
    ty = np.round((x-tx)/2)
    cups = []
    for i in xrange(-1, 2):
        for j in xrange(-1, 2):
            cups.append((2*(i+ty) + (j+tx), np.sqrt(3)*(j+tx)))
    return cups

def plot_hex_grid(h, r, filename):
    tmp = np.fromfile(filename, dtype='float64')
    z = tmp.reshape(tmp.shape[0]/3, 3)

    circs = set([ cup for pt in z for cup in cupnear(pt[0], pt[1]) if cup not in cuplist])
    cuplist.extend(circs)
    for circ in circs:
        pl.gca().add_artist(pl.Circle(circ, h+r, color='k', alpha=0.3))
        pl.gca().add_artist(pl.Circle(circ, h-r, color='w', alpha=1.0))

    pl.plot(z[:,0], z[:,1], linewidth=1)
    pl.gca().set_aspect('equal', 'datalim') 
    pl.xticks([])
    pl.yticks([])
