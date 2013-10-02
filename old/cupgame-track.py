import numpy as np
import scipy as sp
import pylab as pl
from scipy.optimize import fminbound
import matplotlib as mpl
import mpl_toolkits.mplot3d as mp3d
import mpl_toolkits.mplot3d.axes3d as a3

def position(x0, v0, t):
    a = np.array([0,0,-1])
    return x0+v0*t+0.5*a*t**2

def velocity(v0, t):
    a = np.array([0,0,-1])
    return v0 + a*t

def normpos(r):
    tr = r.copy()
    tr[:2] -= 2*((tr[:2].astype('int')+(tr[:2]>0))/2)
    return tr

def distance_to_torus2(t, x0, v0, hsq, rsq):
    r = position(x0, v0, t)
    r = normpos(r)
    return abs((r.dot(r) + hsq - rsq)**2 - 4*hsq*(r[0]*r[0] + r[1]*r[1]))

def distance_to_torus(r, hsq, rsq):
    r = normpos(r)
    return (r.dot(r) + hsq - rsq)**2 - 4*hsq*(r[0]*r[0] + r[1]*r[1])

def inside_torus(r, hsq, rsq):
    return np.sign(distance_to_torus(r, hsq, rsq))

def rotate_about_axis(r, axis, theta):
    return r*np.cos(theta) + np.cross(axis, r)*np.sin(theta) + axis*(axis.dot(r))*(1-np.cos(theta))

def plotPoints(x,y,z):
    fig = pl.figure()
    axs = a3.Axes3D(fig)
    axs.scatter(x,y,z)
    pl.show()
    
def plotPorcupine(lines, x, y, z):
    segs = mp3d.art3d.Line3DCollection(lines)
    fig = pl.figure()
    axs = a3.Axes3D(fig)
    axs.scatter(x,y,z)
    axs.add_collection3d(segs)
    pl.show()

def orbits():
    m = 1.0
    N = 50
    r = 0.3
    h = 1-r

    hsq, rsq = h*h, r*r

    x0 = np.array([0, 0, 2])
    v0 = np.array([0, 0, 0])
    bounces = np.zeros((N,N))
    inout = np.zeros((N,N))
    angles = np.mgrid[0:m:N*1j, 0:m:N*1j].T
   
    segs, xs, ys, zs = [], [], [], []
    for i in xrange(0,50,5):
        for j in xrange(0,50,5):
            pos = x0+np.array([angles[i,j,0], angles[i,j,1], 0])
            vel = v0
            tbounces = 0
            while 1:
                # get the next collision
                px0,py0,pz0 = pos
                vx0,vy0,vz0 = vel
                dt, t, dist = 5e-2, 0, 1e10
                maxspeed = max([vz0, np.sqrt(vz0**2 + 2*pz0)])
                z = pz0

                sign = inside_torus(pos, hsq, rsq)
                prevsign = sign

                while z > 0 and prevsign == sign:
                    prevsign = sign
                    t = t + dt
                    x = px0 + vx0*t
                    y = py0 + vy0*t
                    z = pz0 + vz0*t - 0.5*t**2

                    sign = inside_torus(np.array([x,y,z]), hsq, rsq)
                    xs.append(x); ys.append(y); zs.append(z)

                if sign != prevsign:
                    #ts = np.arange(t-dt, t, 0.001)
                    #pl.plot(ts, [distance_to_torus2(tt, pos, vel, hsq, rsq) for tt in ts], '-')
                    topt, fopt, iter, funcalls = fminbound(distance_to_torus2, t-dt, t, 
                        args=(pos, vel, hsq, rsq), disp=0, full_output=1, xtol=1e-5)
                else:
                    topt = t
                    fopt = distance_to_torus(np.array([x,y,z]), hsq, rsq) 

                tcoll, closest_approach = topt, fopt 

                # figure out where it hit and what speed
                pos = position(pos, vel, tcoll)
                vel = velocity(vel, tcoll)
                xs.append(pos[0])
                ys.append(pos[1])
                zs.append(pos[2])

                if closest_approach < 2e-4:
                    tbounces += 1

                    # findind the normal
                    x,y,z = normpos(pos)
                    theta = np.arctan2(y,x)
                    yin = h*np.sin(theta)
                    xin = h*np.cos(theta)
                    dr = pos-np.array([xin, yin, 0])
                    dr = dr / np.sqrt(dr.dot(dr))
                    ax = dr
                 
                    segs.append((pos, pos-vel/10))
                    segs.append((pos, pos+ax/3))

                    # rotate the velocity about the normal
                    vel = -rotate_about_axis(vel, ax, np.pi)
                    segs.append((pos, pos+vel/10))

                    # go forward in time slightly
                    pos = position(pos, vel, 2e-2)
                    #segs.append((pos, vel/100))

                    if pos[2] < 0:
                        break
                else:
                    break
   
    plotPorcupine(segs,xs,ys,zs)
    #plotPoints(xs,ys,zs)
    #pl.figure()
    #pl.imshow(out)
    return bounces, inout

