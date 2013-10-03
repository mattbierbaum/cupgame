import numpy as np
import scipy as sp
import pylab as pl
import matplotlib as mpl
import mpl_toolkits.mplot3d as mp3d
import mpl_toolkits.mplot3d.axes3d as a3
from scipy import weave

import sys, glob
sys.path.append(glob.glob("../build/lib*")[0])
import cupgamelib

def plotSquareGrid(traj, h, r, fig=None):
    u,v = np.linspace(0, 2*np.pi, 30), np.linspace(0, 2*np.pi, 30)
    ux, vx = np.meshgrid(u,v)
    x,y,z = (h + r*np.cos(vx))*np.cos(ux), (h + r*np.cos(vx))*np.sin(ux), r*np.sin(vx)

    fig = fig or pl.figure()
    axs = a3.Axes3D(fig)
    
    segs = ([ (xt,yt,zt) for xt,yt,zt in traj ],)
    lines = a3.art3d.Line3DCollection(segs)

    axs.add_collection3d(lines)
    for i in xrange(-2,4,2):
        for j in xrange(-2,4,2):
            plot = axs.plot_surface(x+i, y+j, z, rstride = 1, cstride = 1, cmap = pl.get_cmap('jet'),
                                   linewidth = 0, antialiased = False)
    return fig

def plotPoints(x,y,z):
    fig = pl.figure()
    axs = a3.Axes3D(fig)
    axs.scatter(x,y,z)
    pl.show()

def testscan(height=1.0, m=1.0, N=500, r=0.3, h=0.7, d=0.9, disp=0):
    x0 = np.array([0, 0, height]).astype('double')
    v0 = np.array([0, 0, -1e-1]).astype('double')
    out = np.zeros((N,N))
    angles = np.mgrid[-m:m:N*1j, -m:m:N*1j].T

    segs, xs, ys, zs = [], [], [], []
    for i in xrange(N):
        if disp:
            print i
        for j in xrange(N):
            xt = x0+np.array([angles[i,j,0], angles[i,j,1], 0])
            pos = xt.copy(); vel = v0.copy()

            tz = np.array([0]).astype('double')
            tz[0] = cupgamelib.singleCollisions(pos, vel, h, r)
            xs.append(pos[0])
            ys.append(pos[1])
            zs.append(pos[2])
   
    plotPoints(xs,ys,zs)

def corbits(height=1.0, m=1.0, N=500, r=0.3, h=0.7, d=0.9, disp=0):
    x0 = np.array([0, 0, height]).astype('double')
    v0 = np.array([0, 0, -1e-1]).astype('double')
    bounces = np.zeros((N,N))
    angles = np.mgrid[0:m:N*1j, 0:m:N*1j].T
  
    for i in xrange(N):
        if disp == 1:
            print "Slice",i
        for j in xrange(N):
            offset = np.array([angles[i,j,0], angles[i,j,1], 0])
            if np.sqrt(offset[0]**2 + offset[1]**2)  - h > r:
                continue

            pos = x0+offset
            vel = v0

            bounces[i,j] = cupgamelib.trackCollisions(pos, vel, h, r, d, 10000)
  
    return bounces

def slowraise(prefix="rising"):
    start = 0.5
    end = 1.0
    N = 500
    for i in xrange(N):
        print i
        height = start + (end-start)/N*i
        b = corbits(d=height, disp=1)
        pl.figure(10)
        pl.imshow(np.log(b+1), interpolation='nearest', origin='lower')
        pl.xticks([])
        pl.yticks([])
        pl.savefig(prefix+"%04d.png" % i)

def trajectory(i,j, height=1.0, m=1.0, N=500, r=0.3, h=0.7, d=0.9, disp=0):
    x0 = np.array([0, 0, height]).astype('double')
    v0 = np.array([0, 0, -1e-1]).astype('double')
    angles = np.mgrid[0:m:N*1j, 0:m:N*1j].T

    with open("cupgamelib.c") as f:
        import hashlib
        clib = "\n".join(f.readlines())
        sha1 = hashlib.sha1(clib)
        uuid = sha1.hexdigest()

    offset = np.array([angles[i,j,0], angles[i,j,1], 0])
    pos = x0+offset
    vel = v0

    # set up return structures
    mlen = 10000
    tbounces = np.array([0])
    traj = np.zeros((3*mlen*(mlen>0)))
    clen = np.zeros((1,)).astype('int32')

    weave.inline("tbounces[0] = trackTrajectory(pos, vel, h, r, d, 10000, traj, clen, mlen); char id[] = \"%s\";" % uuid, 
            ["tbounces", "pos", "vel", "h", "d", "r", "traj", "clen", "mlen"], support_code=clib, libraries=["m"])
  
    traj = traj.reshape((mlen, 3))
    plotSquareGrid(traj ,h, r)
    #plotPoints(traj[:,0], traj[:,1], traj[:,2])
