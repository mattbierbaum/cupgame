from numpy import *
from scipy import *
from scipy.misc import *
from subprocess import check_call
import os
import Image

z = fromfile("./cupgame_large.npy", dtype='uint8').reshape(16384, 16384, 4)

for zoom in xrange(6):
    try:
        check_call(['mkdir', str(zoom)])
    except:
        pass

    for x in xrange(1 << zoom):
        try:
            check_call(['mkdir', os.path.join(str(zoom), str(x))])
        except:
            pass

        for y in xrange(1 << zoom):
            ff = 256*(1 << 6 - zoom)
            yy = (1 << 6 - zoom - 1)
            tt = z[ff*y : ff*(y+1) : yy, ff*x : ff*(x+1) : yy,:]
            a = Image.fromarray(tt)
            a.save(os.path.join(os.path.join(str(zoom), str(x)), str(y)+'.png'))
