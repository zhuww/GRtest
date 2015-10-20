import numpy as np
import matplotlib.pyplot as plt
from pylab import rand
from cmath import *

from psrpbdot import M1
from datatools.tempo import *
from astropy import coordinates as coord
from tools.Coordinate import *
import numpy.linalg as linalg

#from Arrow3D import Arrow3D
secperday = 24*3600
solardist = 8.34

from GalacticGeometry import *
#pf = PARfile('1713.Dec.mcmc.par')
pf = PARfile('mcmcresult.par')
#pf = PARfile('1713.sns.par')

gl, gb = getGpos(pf)
D = 1./float(pf.PX[0])
x = D * np.cos(gb) * np.sin(gl)
y = solardist - D * np.cos(gb) * np.cos(gl) 
#angle = np.arctan(y/x)*180./np.pi

T, GT = GetTransferMatrix(pf)
pos1713 = coord.FK5Coordinates(str(COORD.RA(pf.RAJ[0]))+' '+str(COORD.Dec(pf.DECJ[0])))
GCpos = coord.FK5Coordinates('17h45m40.04s -29d00m28.1s')

Eerr = float(pf.E[1])
OMerr = float(pf.OM[1])/180.*np.pi

G = 6.673e-8
c = 2.99792458e10
PI = np.pi
AU = 1.469e13
Msun = 1.9882e33
Tsun = 4.925490947e-6
secperday = 24*3600
R0 = solardist
DRA = pos1713.ra.radians - GCpos.ra.radians
DDec = pos1713.dec.radians - GCpos.dec.radians
#Kz = lambda z:(2.27*z + 3.68*(1-np.exp(-4.31*z)) ) * 1.e-9 #Galactic acceleration in z direction (cm/s^2)
def Kz(z):
    absz = np.abs(z)
    sign = z/absz
    return sign * (2.27*absz + 3.68*(1-np.exp(-4.31*absz)) ) * 1.e-9


def getKG(kr, zeta, z, sini, paascnode, om):
    #g_r = kr * (GT.I * np.matrix((np.sin(0.-zeta),np.cos(0.-zeta),0)).T) #X/linalg.norm(X)
    g_r = GT.I * kr * ( np.matrix((np.cos(0.-zeta),np.sin(0.-zeta),0)).T) #X/linalg.norm(X)
    #g_z = -1. * Kz(z) * Gc
    g_z = GT.I * Kz(z) * (np.matrix((0., 0., -1.)).T) 
    g = g_r + g_z
    g_NSEW = T * g
    #print 'GT:', GT
    #print 'g:', g
    #print 'T:', T
    #print 'g_NSEW:', g_NSEW

    incang = np.arcsin(sini)
    Omgang = paascnode/180.*np.pi
    A = -1./np.tan(incang)
    B = -1./np.tan(Omgang)
    C = 1.

    n_orb = np.matrix((A, B, C)).T
    n_orb= n_orb/linalg.norm(n_orb)
    #print 'n_orb', n_orb.T

    g_proj = g_NSEW - n_orb * (g_NSEW.T*n_orb) 
    KG  = float(linalg.norm(g_proj))
    A_ref = np.matrix((0, -1.* np.sin(Omgang), np.cos(Omgang)))
    g_dir = g_proj/KG
    #print 'A_ref,', A_ref, 'g_dir', g_dir
    g_ang = np.arccos(A_ref * g_dir)
    #print 'g angle on orbit:', g_ang*180/np.pi, A_ref * g_dir
    #print 'omega: ', om
    #print 'theta:', float(g_ang - om/180.*np.pi) * 180./np.pi
    #g_ecc = np.cos(g_ang - om/180.*np.pi) * KG
    g_norm = linalg.norm(g)
    #print 'ratio:', KG/g_norm, 'KG:', KG, 'g:', g_norm
    return float(g_ang - om/180.*np.pi), KG

def EccArea(ECC, EF, THETA):
    global Eerr, OMerr
    #THETA[THETA<0]+=(np.pi*2)
    #print 'mean theta', np.mean(THETA)*180/np.pi
    Emax, Emin = ECC +Eerr*3 , ECC-Eerr*3
    Temax, Temin = THETA +OMerr*3, THETA-OMerr*3
    fourpnts = [[(x*np.cos(y)-EF[i]+x*np.sin(y)*1j) for x in (Emax[i], Emin[i]) for y in (Temax[i], Temin[i])] for i in range(len(ECC))]
    areas = []
    for v4 in fourpnts:
        abs4 = [abs(v) for v in v4]
        p4 = np.array([phase(v) for v in v4])
        if p4.min() < 0. and p4.max() - p4.min() > np.pi:
            p4[p4<0.] += 2.*np.pi
        if np.mean(abs4) < Eerr*3 and abs(np.mean(p4))<OMerr*3:
            #print "hei, it's happening"
            value = max((np.log10(min(max(abs4), 0.05)) - np.log10(max(min(abs4),1.e-6))),0)*np.pi*2
        else:
            value = max((np.log10(min(max(abs4), 0.05)) - np.log10(max(min(abs4),1.e-6))),0)*(max(p4)-min(p4))
        #print abs4, p4
        areas.append(value)
    areas = np.array(areas)
    #print 'EF, ECC, THETA, AREA, A4, p4', EF.mean(), ECC.mean(), THETA.mean(), areas.mean(), max(abs4),min(abs4), max(p4), min(p4)
    #print 'Eerr, OMerr', Eerr,OMerr
    return areas

def Pintegrant(PX, SINI, PAASCNODE, M1, M2, PB, ECC, OM, Delta):
    D = 1./PX
    Omega = PAASCNODE/180.*np.pi
    #Theta_g = np.pi - np.arctan(np.tan(DRA)/np.sin(DDec))
    z = np.sin(gb) * D
    R1 = np.sqrt(R0**2 + (D*np.cos(gb))**2 -2 * R0 * D * np.cos(gb) * np.cos(gl))
    coszeta = (R0**2 + R1**2 - D**2 + z**2)/R0/R1/2.
    zeta = np.arccos(coszeta)
    #print 'zeta: ',zeta

    Mtot =  M1 + M2
    Omega_G = 27.2 #km s^-1 kpc^-1
    kpcinkm = 3.0857e16
    Kr =  Omega_G**2 * R0**2 / R1 / kpcinkm * 1.e5 #Galactic acceleration in radio direction (cm/s^2)
    KGarray = []
    THETA = []
    for i,sini in enumerate(SINI):
        theta, kg= getKG(Kr[i], zeta[i], z[i], sini, PAASCNODE[i], OM[i])
        THETA.append(theta)
        #THETA.append(np.abs(theta-np.pi))
        KGarray.append(kg)
    THETA = np.array(THETA)
    THETA[THETA<0]+=(np.pi*2)
    KG = np.array(KGarray)
    #print 'Correct projected Galactic acceleration: ', KG
    EF = Delta * ( 0.5 * KG / Mtot /Tsun/c/(2*PI/PB)**2 )
    #EF = Delta * ( 0.5 * KG * c**2 / G / Mtot / Msun /(2*PI/PB)**2 )
    #return  ECC*xi/( 0.5 * KG * c**2 / G / Mtot / Msun /(2*PI/PB)**2 )
    #print 'Delta:', Delta,
    Areas = 0.
    Areas += sum(EccArea(ECC, EF, THETA))
    Areas += sum(EccArea(ECC, EF, THETA-np.pi))
    return Areas

''' load in the MCMC results for Delta estimation'''
import cPickle as pickle
import sys
#from Coordinate import RA, Dec
secperday = 24*3600

dic = pickle.load(open('bestpar.p', 'rb'))
best = dic['BEST']
plist = dic['parameters'] 
MChain = pickle.load(open('TinyMChain.p','rb'))
#MChain = pickle.load(open('SmallMChain.p','rb'))
MarkovChain = MChain['Chain']
MCMCSize = len(MarkovChain)
#MChain = pickle.load(open('MChain.p','rb'))
#MarkovChain = MChain['Chain']
#MCMCSize = len(MarkovChain)
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
#xi = np.array([ Xi(y) for y in npr.uniform(0., 2*np.pi, MCMCSize)])

#print 'Ecc:', ECC
#print 'M1', M1, 'M2', M2

#integ = Pintegrant(PX, SINI, PAASCNODE, M1, M2, PB, ECC, OM, 0.001)
Integ = lambda d: Pintegrant(PX, SINI, PAASCNODE, M1, M2, PB, ECC, OM, d)
#Integ(1.e-5)
#sys.exit(0)
#print integ

#dsize = delta.size
#delta.sort()
#print delta[int(dsize*0.95)] , delta[int(dsize*0.05)]

#delta = np.arange(1.e-5, 1.e-1, 2.e-5) #ingrid setting
delta = np.arange(5.e-5, 0.03, 5.e-5) #ingrid setting
#delta = np.arange(1.e-8, 1.e-4, 5.e-8)
res = np.array([Integ(d) for d in delta])

np.save(open('SEPresult.npy', 'wb'), res)

cdf = [res[0]]
for r in res[1:]:
    s = cdf[-1] + r
    cdf.append(s)
sumres = cdf[-1]
cdf = np.array(cdf)/sumres
for i,c in enumerate(cdf):
    if c > 0.95:
        print delta[i], c
        break

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(delta, res, '-')
#ax.logx()
#ax.hist(delta, bins=50, normed=1)
ax.semilogx(nonposy='clip')
ax.semilogy(nonposy='clip')
plt.show()

