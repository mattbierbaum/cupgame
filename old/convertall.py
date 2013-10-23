import matplotlib as mpl
mpl.use('Agg')  # use this for no display running
import pylab as pl
import numpy as np
from scipy import misc
import glob, re

def plotRawPicture(field, name, cmap=None):
    mapper = cmap or mpl.cm.jet 
    q = mapper(mpl.colors.Normalize()(np.log(field+1)))
    misc.imsave(name, q)

for f in glob.glob("./temp*.npy"):
    index = re.search("temp_([0-9]*).npy", f)
    a = np.fromfile(f, dtype='int32').reshape(512,512)
    aa = np.ma.masked_array(a, a==0)
    plotRawPicture(aa, "./picture_%04d.png" % int(index.groups()[0]), cmap=pl.cm.ocean_r)
