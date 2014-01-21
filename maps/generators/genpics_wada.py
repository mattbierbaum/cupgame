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
outdir = '/var/www/runat.me/tmp/maps/wada'

steps = 30
size0 = 1.0

targetx = 0.7800808571829128
targety = 0.2944330689467324

for zoom in xrange(155):
    tfile = '/media/scratch/cuptemp_wada'
    N = 1024
    
    width = size0 * 0.8**zoom

    print zoom, width, targetx - width/2, targety - width/2

    args = ['./cupgame-nowait-mapgen', N, 2.0, 1.00,
            targetx - width/2, targetx + width/2,
            targety - width/2, targety + width/2, tfile]
    
    check_call([str(a) for a in args])
    tt = np.fromfile(tfile, dtype='int32').reshape(N,N)
    out = cm.flag(Normalize(vmin=-3000, vmax=3000)(tt))
    out[tt == 0] = [1.0,1.0,1.0,1.0]

    imsave(join(outdir, 'wada_%04d.png' % zoom), out)
