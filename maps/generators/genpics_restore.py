import numpy as np
import scipy as sp
import scipy.ndimage as nd
from scipy.misc import *
from subprocess import check_call
import os
import Image
import matplotlib as mpl
mpl.use('Agg')
import pylab as pl

from matplotlib.colors import Normalize
from matplotlib import cm

join = os.path.join
outdir = '/var/www/runat.me/tmp/maps/energy'

steps = 30

for i, height in enumerate(np.arange(0.5, 1.0, 0.005)):
    tfile = '/media/scratch/cuptemp_energy'
    N = 1024 

    x0 = 0.4*np.cos(30./180 * np.pi)
    size = 1 - x0

    args = ['./cupgame-nowait-mapgen', N, 2.0, height, 
            x0, size + x0, 
            0., size + 0., tfile]
    
    check_call([str(a) for a in args])
    tt = np.fromfile(tfile, dtype='int32').reshape(N,N)
    tt = np.ma.masked_array(tt, tt==0)
    out = cm.jet(Normalize(vmin=np.log(2), vmax=np.log(160))(np.log(tt+1)))

    imsave(join(outdir, 'restore_%04d.png' % i), out)
