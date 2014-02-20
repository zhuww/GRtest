import numpy as np
import matplotlib.pyplot as plt
from pylab import rand

from psrpbdot import M1
from datatools.tempo import *
from astropy import coordinates as coord
from tools.Coordinate import *
import numpy.linalg as linalg

#from Arrow3D import Arrow3D

G = Decimal(6.673e-8)
c = Decimal(2.99792458e10)
PI = Decimal(np.pi)
AU = Decimal(1.469e13)
Msun = Decimal(1.9882e33)
secperday = 24*3600
solardist = 8.33

from GalacticGeometry import *
pf = PARfile('1713.Dec.mcmc.par')

gl, gb = getGpos(pf)
D = 1./float(pf.PX[0])
x = D * np.cos(gb) * np.sin(gl)
y = solardist - D * np.cos(gb) * np.cos(gl) 
angle = np.arctan(y/x)*180./np.pi

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

def Delta(PX, SINI, PAASCNODE, M1, M2, PB, ECC, OM, xi):
    G = 6.673e-8
    c = 2.99792458e10
    PI = np.pi
    AU = 1.469e13
    Msun = 1.9882e33
    secperday = 24*3600
    D = 1./PX
    DRA = pos1713.ra.radians - GCpos.ra.radians
    DDec = pos1713.dec.radians - GCpos.dec.radians
    Omega = PAASCNODE/180.*np.pi
    Theta_g = np.pi - np.arctan(np.tan(DRA)/np.sin(DDec))

    z = np.sin(gb) * D
    R0 = solardist
    R1 = np.sqrt(R0**2 + (D*np.cos(gb))**2 -2 * R0 * D * np.cos(gb) * np.cos(gl))
    coszeta = (R0**2 + R1**2 - D**2 + z**2)/R0/R1/2.
    zeta = np.arccos(coszeta)

    Mtot =  M1 + M2
    Kz = lambda z:(2.27*z + 3.68*(1-np.exp(-4.31*z)) ) * 1.e-9 #Galactic acceleration in z direction (cm/s^2)
    Omega_G = 27.2 #km s^-1 kpc^-1
    kpcinkm = 3.0857e16
    Kr =  Omega_G**2 * R1 / kpcinkm * 1.e5 #Galactic acceleration in radio direction (cm/s^2)

    def getKG(kr, zeta, z, sini, paascnode, om):
        g_r = kr * (GT.I * np.matrix((np.sin(0.-zeta),np.cos(0.-zeta),0)).T) #X/linalg.norm(X)
        #g_z = -1. * Kz(z) * Gc
        g_z = GT.I * Kz(z) * (np.matrix((0., 0., -1.)).T) 
        g = g_r + g_z
        g_NSEW = T * g
        #print g_NSEW

        incang = np.arcsin(sini)
        Omgang = paascnode/180.*np.pi
        A = -1./np.tan(incang)
        B = -1./np.tan(Omgang)
        C = 1.

        n_orb = np.matrix((A, B, C)).T
        n_orb= n_orb/linalg.norm(n_orb)
        #print n_orb

        g_proj = g_NSEW - n_orb * (g_NSEW.T*n_orb) 
        KG  = float(linalg.norm(g_proj))
        A_ref = np.matrix((0, -1.* np.sin(Omgang), np.cos(Omgang)))
        g_dir = g_proj/KG
        g_ang = np.arccos(A_ref * g_dir)
        g_ecc = np.cos(g_ang - om/180.*np.pi) * KG

        return float(np.abs(g_ecc))
    KGarray = []
    for i,sini in enumerate(SINI):
        KGarray.append(getKG(Kr[i], zeta[i], z[i], sini, PAASCNODE[i], OM[i]))
    KG = np.array(KGarray)
    #print 'Correct projected Galactic acceleration: ', KG
    return  ECC*xi/( 0.5 * KG * c**2 / G / Mtot / Msun /(2*PI/PB)**2 )


import numpy.random as npr

''' load in the MCMC results for Delta estimation'''
import cPickle as pickle
import sys
#from Coordinate import RA, Dec
secperday = 24*3600

dic = pickle.load(open('bestpar.p', 'rb'))
best = dic['BEST']
plist = dic['parameters'] 
MChain = pickle.load(open('MChain.p','rb'))
MarkovChain = MChain['Chain']
MCMCSize = len(MarkovChain)
pi = 3.141592653589793
#G = 6.673e-11
#Msun = 1.98892e30
#c = 2.99792458e8
twopi = 6.283185307179586
#fac = 1.536e-16 
G = 6.673e-8
c = 2.99792458e10
Msun = 1.9882e33
im2 = plist.index('M2')
ipb = plist.index('PB')
isini = plist.index('SINI')
ia = plist.index('A1')
ichisq = plist.index('chisq')
iecc = plist.index('E')
ipx = plist.index('PX')
iomega = plist.index('PAASCNODE')
iom = plist.index('OM')
M2 = np.array([float(p[im2])*Msun for p in MarkovChain])
PB = np.array([float(p[ipb])*secperday for p in MarkovChain])
SINI = np.array([float(p[isini]) for p in MarkovChain])
OM = np.array([float(p[iom]) for p in MarkovChain])
a = np.array([float(p[ia])*c for p in MarkovChain])
#print PB, M2, SINI, 
M1 = (PB/2/pi*np.sqrt(G*(M2*SINI)**3/a**3)-M2)/Msun
M2 = M2/Msun
ECC = np.array([float(p[iecc]) for p in MarkovChain])
PX = np.array([float(p[ipx]) for p in MarkovChain])
PAASCNODE = np.array([float(p[iomega]) for p in MarkovChain])
chisq = [p[ichisq] for p in MarkovChain]
bestidx = chisq.index(min(chisq))
xi = np.array([ Xi(y) for y in npr.uniform(0., 2*np.pi, MCMCSize)])

#print 'Ecc:', ECC
#print 'M1', M1, 'M2', M2

delta = Delta(PX, SINI, PAASCNODE, M1, M2, PB, ECC, OM, xi)

dsize = delta.size
delta.sort()
print delta[int(dsize*0.95)] , delta[int(dsize*0.05)]

fig = plt.figure()
ax = fig.add_subplot(111)
ax.hist(delta, bins=50, normed=1)
ax.semilogy(nonposy='clip')
plt.show()

