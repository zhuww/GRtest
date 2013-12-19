import numpy as np
import matplotlib.pyplot as plt
from pylab import rand
from datatools.tempo import *
from astropy import coordinates as coord
from tools.Coordinate import *
from psrpbdot import M1

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
#fig = plt.figure(figsize=(10,10))
#ax = fig.add_subplot(111, aspect='equal')

G = Decimal(6.673e-11)
c = Decimal(2.99792458e10)
PI = Decimal(np.pi)
AU = Decimal(1.469e13)
Msun = Decimal(1.9882e33)
secperday = 24*3600


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

GCpos = coord.FK5Coordinates('17h45m40.04s -29d00m28.1s')
pos1713 = coord.FK5Coordinates(str(RA(pf.RAJ[0]))+' '+str(Dec(pf.DECJ[0])))
DRA = pos1713.ra.radians - GCpos.ra.radians
DDec = pos1713.dec.radians - GCpos.dec.radians
Omega = float(pf.PAASCNODE)/180.*np.pi
Theta_g = np.pi - np.arctan(np.tan(DRA)/np.sin(DDec))

#print Theta_g*180/np.pi, Omega*180/np.pi

z = np.sin(gb) * D
R0 = solardist
R1 = np.sqrt(R0**2 + (D*np.cos(gb))**2 -2 * R0 * D * np.cos(gb) * np.cos(gl))
lbd = np.arccos((R1**2 + D**2 + z**2 - R0**2)/(2*D*np.sqrt(R1**2 + z**2)))#*180/np.pi
print 'lambda: ', lbd, 'R_PSR: ', R1, 'z: ', z
rat = np.sqrt(1 - (np.cos(i)*np.cos(lbd) + np.sin(i)*np.sin(lbd)*np.sin(Theta_g - Omega))**2)
print 'ratio:',rat
Mtot = M1(pf) + pf.M2[0]
print 'M1+M2', Mtot
Kz = lambda z:(2.27*z + 3.68*(1-np.exp(-4.31*z)) ) * 1.e-9 #Galactic acceleration in z direction (cm/s^2)
Omega_G = 27.2 #km s^-1 kpc^-1
#R_G = 8.33 # +/-0.35 kpc
kpcinkm = 3.0857e16
Kr =  Omega_G**2 * R1 / kpcinkm * 1.e5 #Galactic acceleration in radio direction (cm/s^2)
KG = np.sqrt(Kr**2 + Kz(z)**2) 
print 'Galactic acceleration: ',Kz(z), Kr, KG
print 'Projected Galactic acceleration: ', KG*rat

print 'EF:', Decimal(0.009) * Decimal(0.5) * Decimal(KG*rat) * c**2 / G / Mtot / Msun /(2*PI/pf.PB[0]/secperday)**2

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
art3d.pathpatch_2d_to_3d(e, z=z, zdir=(0,1,0))
e.set_alpha(1.)
e.set_edgecolor('r')
plt.show()



