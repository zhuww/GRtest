import numpy as np
import matplotlib.pyplot as plt
from pylab import rand
from datatools.tempo import *
from astropy import coordinates as coord
from tools.Coordinate import *

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
#fig = plt.figure(figsize=(10,10))
#ax = fig.add_subplot(111, aspect='equal')

data = np.genfromtxt(open('log_arms.out', 'r'), dtype = [('Num', 'i'), ('n', 'i'), ('x', 'f8'),('y', 'f8')])[1:]
UniqNum = np.unique(data['Num'])

solardist = 8.3
for num in UniqNum:
    x = data['x'][data['Num'] == num]
    y = data['y'][data['Num'] == num]
    ax.plot(x,y,0,'b-')
    #ax.plot([0],[solardist], 0.17, 'yo', ms=15) #OLausen & Kaspi
    ax.plot([0],[solardist], 0., 'yo', ms=5) #OLausen & Kaspi
    ax.plot([0],[0], 0, 'k*', ms=15)
    #ax.text(0,0,'GC')
    ax.set_xlabel('x (kpc)')
    ax.set_ylabel('y (kpc)')
    #ax.set_xlim((-4.0,6.8))
    #ax.set_ylim((-0.2,10.6))

def getGpos(pf):
    ra = RA(pf.RAJ[0])
    dec = Dec(pf.DECJ[0])
    pos = coord.FK5Coordinates(str(ra) +' '+ str(dec))
    l = pos.galactic.l.radians
    b = pos.galactic.l.radians
    return l,b

pf = PARfile('1713.Oct.mcmc.par')

#print getGpos(pf)
gl, gb = getGpos(pf)
D = 1./float(pf.PX[0])
x = D * np.cos(gb) * np.sin(gl)
y = solardist - D * np.cos(gb) * np.cos(gl) 
angle = np.arctan(y/x)*180./np.pi


from matplotlib.patches import Ellipse, Circle
import mpl_toolkits.mplot3d.art3d as art3d
from matplotlib.collections import PatchCollection
#ells = [Ellipse(xy=[x[i], y[i]], width=(0.6+0.4*rand())*(1.-0.5), height=0.6*1, angle=90, fill=False, lw=2)  for i in range(len(psrs))]
e = Ellipse(xy = [x,y], width=(0.9), height=0.2, angle=90, fill=True, lw=1)
#ax.add_artist(ell)
#p = PatchCollection([e])
#ax.add_collection3d(p, zs=[0])
#p = Circle((1,1),1)
#ax.add_patch(p)
ax.add_patch(e)
art3d.pathpatch_2d_to_3d(e, z=0, zdir=(0,1,0))
e.set_alpha(1.)
e.set_edgecolor('r')
plt.show()


