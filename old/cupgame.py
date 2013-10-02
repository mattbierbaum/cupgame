import numpy as np
import scipy as sp
import pylab as pl
from scipy.optimize import fmin
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
    rt = r.copy()
    rt[:2] -= 2*((rt[:2].astype('int')+(rt[:2]>0))/2)
    return rt

def distance_to_torus2(t, x0, v0, hsq, rsq):
    r = position(x0, v0, t)
    r = normpos(r)
    return abs((r.dot(r) + hsq - rsq)**2 - 4*hsq*(r[0]*r[0] + r[1]*r[1]))

def distance_to_torus(r, hsq, rsq):
    r = normpos(r)
    return (r.dot(r) + hsq - rsq)**2 - 4*hsq*(r[0]*r[0] + r[1]*r[1])

def inside_torus(r, hsq, rsq):
    return np.sign(distance_to_torus(r, hsq, rsq))

def collision_time(x0, v0, hsq, rsq):
    px0,py0,pz0 = x0
    vx0,vy0,vz0 = v0
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
    return topt,fopt 

def collision_normal(pos, vel, h):
    x,y,z = normpos(pos)
    theta = np.arctan2(y,x)
    yin = h*np.sin(theta)
    xin = h*np.cos(theta)
   
    dr = pos-np.array([xin, yin, 0])
    dr = dr / np.sqrt(dr.dot(dr))
    return dr #(pos, pos+dr)

def rotate_about_axis(r, axis, theta):
    return r*np.cos(theta) + np.cross(axis, r)*np.sin(theta) + axis*(axis.dot(r))*(1-np.cos(theta))

def reflect_velocity(pos, vel, h):
    ax = collision_normal(pos, vel, h)
    return -rotate_about_axis(vel, ax, np.pi)

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

def testscan():
    m = 2.0
    N = 40
    r = 0.3
    h = 1-r

    hsq, rsq = h*h, r*r

    x0 = np.array([0, 0, 1])
    v0 = np.array([0, 0, 0])
    out = np.zeros((N,N))
    angles = np.mgrid[-m:m:N*1j, -m:m:N*1j].T
   
    segs, xs, ys, zs = [], [], [], []
    for i in xrange(N):
        for j in xrange(N):
            xt = x0+np.array([angles[i,j,0], angles[i,j,1], 0])
            tcoll, closest_approach = collision_time(x0=xt, v0=v0, hsq=hsq, rsq=rsq)

            if abs(closest_approach) < 1e-4:
                pos = position(xt, v0, tcoll)
                vel = velocity(v0, tcoll)
                #segs.append(collision_resultant(pos, vel, h))
                xs.append(pos[0])
                ys.append(pos[1])
                zs.append(pos[2])
                out[i,j] = pos[2]
   
    #plotPorcupine(segs,xs,ys,zs)
    plotPoints(xs,ys,zs)
    #pl.figure()
    #pl.imshow(out)
    return out

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
    for i in xrange(N):
        print i
        for j in xrange(i,N):
            pos = x0+np.array([angles[i,j,0], angles[i,j,1], 0])
            vel = v0
            tbounces = 0
            while 1:
                # get the next collision
                tcoll, closest_approach = collision_time(x0=pos, v0=vel, hsq=hsq, rsq=rsq)

                # figure out where it hit and what speed
                pos = position(pos, vel, tcoll)
                vel = velocity(vel, tcoll)

                if closest_approach < 1e-4:
                    tbounces += 1
                    #segs.append((pos, vel/100))
                    xs.append(pos[0])
                    ys.append(pos[1])
                    zs.append(pos[2])
                    vel = 0.9*reflect_velocity(pos, vel, h)
                    pos = position(pos, vel, 5e-2)
                    #segs.append((pos, vel/100))
                else:
                    pos = normpos(pos)
                    if pos.dot(pos) < hsq:
                        inout[i,j] = 1
                    else:
                        inout[i,j] = 0
                    break
                if tbounces % 500 == 0:
                    print tbounces
                if tbounces > 5000:
                    break
            bounces[i,j] = tbounces
   
    #plotPorcupine(segs,xs,ys,zs)
    plotPoints(xs,ys,zs)
    #pl.figure()
    #pl.imshow(out)
    return bounces, inout

