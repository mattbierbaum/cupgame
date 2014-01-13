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
outdir = '/var/www/runat.me/tmp/maps/'

for zoom in xrange(4, 10):
    try:
        check_call(['mkdir', join(outdir, str(zoom))])
    except:
        pass

    for x in xrange(1 << zoom):
        try:
            check_call(['mkdir', join(outdir, join(str(zoom), str(x)))])
        except:
            pass

        for y in xrange(1 << zoom):
            tfile = '/media/scratch/cuptemp'
            zz = (1 << zoom)
            N = 256

            x0 = 0.4*np.cos(30./180 * np.pi)
            size = 1 - x0

            args = ['./cupgame-nowait-mapgen', N, 2.0, 0.90, 
                    size*x/zz + x0, size*(x+1)/zz + x0, 
                    size*y/zz + 0., size*(y+1)/zz + 0., tfile]
    
            #args = ['./cupgame-nowait-mapgen', N, 2.0, 0.90, 
            #        2.0*x/zz-1.0, 2.0*(x+1)/zz-1.0, 
            #        2.0*y/zz-1.0, 2.0*(y+1)/zz-1.0, tfile]
    
            check_call([str(a) for a in args])
            tt = np.fromfile(tfile, dtype='int32').reshape(N,N)
            tt = np.ma.masked_array(tt, tt==0)
            out = cm.jet(Normalize(vmin=0, vmax=np.log(360))(np.log(tt+1)))

            imsave(join(outdir, join(join(str(zoom), str(x)), str(y)+'.png')), out)
