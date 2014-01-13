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
outdir = '/var/www/runat.me/tmp/maps/zoom'

steps = 30
size0 = 1.0

targetx = 0.7844019803357969
targety = 0.2640594016350469

for zoom in xrange(150):
    tfile = '/media/scratch/cuptemp_zoom'
    N = 128
    
    width = size0 * 0.9**zoom

    print zoom, width, targetx - width/2, targety - width/2

    args = ['./cupgame-nowait-mapgen', N, 2.0, 0.90, 
            targetx - width/2, targetx + width/2,
            targety - width/2, targety + width/2, tfile]
    
    check_call([str(a) for a in args])
    tt = np.fromfile(tfile, dtype='int32').reshape(N,N)
    tt = np.ma.masked_array(tt, tt==0)
    out = cm.jet(Normalize(vmin=np.log(2), vmax=np.log(160))(np.log(tt+1)))
    out = cm.jet(Normalize()(np.log(tt+1)))

    imsave(join(outdir, 'zoom_%04d.png' % zoom), out)
