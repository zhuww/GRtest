import numpy as np
import matplotlib.pyplot as plt
from pylab import rand
from datatools.tempo import *
from astropy import coordinates as coord
from tools.Coordinate import *
from psrpbdot import M1

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
import numpy.linalg as linalg

from Arrow3D import Arrow3D

G = Decimal(6.673e-8)
c = Decimal(2.99792458e10)
PI = Decimal(pi)
AU = Decimal(1.469e13)
Msun = Decimal(1.9882e33)
secperday = 24*3600

def getGpos(pf):
    ra = RA(pf.RAJ[0])
    dec = Dec(pf.DECJ[0])
    pos = coord.FK5Coordinates(str(ra) +' '+ str(dec))
    l = pos.galactic.l.radians
    b = pos.galactic.b.radians
    return l,b

pf = PARfile('1713.Dec.mcmc.par')

from GalacticGeometry import *

T, GT = GetTransferMatrix(pf)
pos1713 = coord.FK5Coordinates(str(COORD.RA(pf.RAJ[0]))+' '+str(COORD.Dec(pf.DECJ[0])))
GCpos = coord.FK5Coordinates('17h45m40.04s -29d00m28.1s')


def Xi(theta):
    theta %= (2*np.pi)
    if theta >= 0. and theta < np.pi/2:
        return 1/np.sin(theta)
    elif theta >= np.pi/2 and theta <=1.5*np.pi:
        return 1
    else:
        return -1./np.sin(theta)

def alpha3(M1, M2, PB, F0, ECC, PMRA, PMDEC, PX, SINI, PAASCNODE, OM, w, xi):
    G = 6.673e-8
    c = 2.99792458e10
    PI = np.pi
    AU = 1.469e13
    Msun = 1.9882e33
    secperday = 24*3600
    secperyear = secperday*365
    kpc = 3.0857e21
    Mtot = M1+M2
    EF = lambda w:0.21* M1 * w * (PB**2) * c**2 *F0/24 / PI /G /Mtot/Msun

    wsolar = 369.e5 #See ref below (aaa+13 Planck Team: Aghanim, N. et al. 2013. Planck confirms this using a different method)
    wserr = 0.9e5 #Kogut et al. 1993, Fixsen et al 1996, Hinshaw et al. 2009
    lws, bws = 263.99/180.*np.pi, 48.26/180.*np.pi
    w_s = wsolar * (GT.I * np.matrix((np.cos(bws)*np.sin(lws),np.cos(bws)*np.cos(lws),np.sin(bws))).T)
    ws_NSEW = T * w_s

    D = kpc/PX
    
    wx = w * (np.matrix((-1., 0., 0.)).T)
    wy = PMRA*1.e-3/60./60.*np.pi/180./secperyear * D * (np.matrix((0.,-1.,0.)).T)
    wz = PMDEC*1.e-3/60./60.*np.pi/180./secperyear * D * (np.matrix((0.,0.,1.)).T)

    #print wx+wy+wz
    #print ws_NSEW
    w = wx + wy + wz + ws_NSEW

    incang = np.arcsin(SINI)
    Omgang = PAASCNODE/180.*np.pi
    A = -1./np.tan(incang)
    B = -1./np.tan(Omgang)
    C = 1.
    A_ref = np.matrix((0, -1.* np.sin(Omgang), np.cos(Omgang)))

    n_orb = np.matrix((A, B, C)).T
    n_orb= n_orb/linalg.norm(n_orb)

    w_proj = w - n_orb * (w.T*n_orb) 
    w_leg = linalg.norm(w_proj)
    w_dir = w_proj/w_leg
    w_ang = np.arccos(A_ref * w_dir)
    w_ecc = np.cos(w_ang - OM/180.*np.pi) * w_leg

    return ECC*xi/EF(linalg.norm(w_ecc)) 
    #return ECC/EF(linalg.norm(w_proj)) 

#M1, M2, PB, F0, ECC, w = 1.3, 0.3, float(pf.PB[0]*secperday), float(pf.F0[0]), float(pf.E[0]), 25.e5
#PMRA, PMDEC, PX = float(pf.PMRA[0]), float(pf.PMDEC[0]), float(pf.PX[0])
#print alpha3(M1, M2, PB, F0, ECC, PMRA, PMDEC, PX, w)

''' load in the MCMC results for Delta estimation'''
import cPickle as pickle
import sys
import numpy.random as npr
#from Coordinate import RA, Dec
secperday = 24*3600

dic = pickle.load(open('bestpar.p', 'rb'))
best = dic['BEST']
plist = dic['parameters'] 
MChain = pickle.load(open('MChain.p','rb'))
MarkovChain = MChain['Chain']
MCMCSize = len(MarkovChain)
pi = 3.141592653589793
G = 6.673e-11
Msun = 1.98892e30
c = 2.99792458e8
twopi = 6.283185307179586
fac = 1.536e-16 
if0 = plist.index('F0')
im2 = plist.index('M2')
ipb = plist.index('PB')
isini = plist.index('SINI')
ia = plist.index('A1')
ichisq = plist.index('chisq')
iecc = plist.index('E')
ipx = plist.index('PX')
ipmra = plist.index('PMRA')
ipmdec = plist.index('PMDEC')
iomega = plist.index('PAASCNODE')
iom = plist.index('OM')
M2 = np.array([float(p[im2])*Msun for p in MarkovChain])
PB = np.array([float(p[ipb])*secperday for p in MarkovChain])
SINI = np.array([float(p[isini]) for p in MarkovChain])
a = np.array([float(p[ia])*c for p in MarkovChain])
#print PB, M2, SINI, 
F0 = np.array([float(p[if0]) for p in MarkovChain])
M1 = (PB/2/pi*np.sqrt(G*(M2*SINI)**3/a**3)-M2)/Msun
M2 = M2/Msun
ECC = np.array([float(p[iecc]) for p in MarkovChain])
PX = np.array([float(p[ipx]) for p in MarkovChain])
PMRA = np.array([float(p[ipmra]) for p in MarkovChain])
PMDEC = np.array([float(p[ipmdec]) for p in MarkovChain])
PAASCNODE = np.array([float(p[iomega]) for p in MarkovChain])
OM = np.array([float(p[iom]) for p in MarkovChain])
chisq = [p[ichisq] for p in MarkovChain]
bestidx = chisq.index(min(chisq))
w = np.array([ y for y in npr.normal(0., 2900000., MCMCSize)])
xi = np.array([ Xi(y) for y in npr.uniform(0., 2*np.pi, MCMCSize)])
alpha = np.array([alpha3(M1[i], M2[i], PB[i], F0[i], ECC[i], PMRA[i], PMDEC[i], PX[i], SINI[i], PAASCNODE[i], OM[i], w[i], xi[i]) for i in range(len(w))])

dsize = alpha.size
alpha.sort()
print alpha[int(dsize*0.95)] , alpha[int(dsize*0.05)]

fig = plt.figure()
ax = fig.add_subplot(111)
ax.hist(alpha, bins=50, normed=1)
ax.semilogy(nonposy='clip')
plt.show()
