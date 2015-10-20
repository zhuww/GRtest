import numpy as np
import matplotlib.pyplot as plt
from pylab import rand
from cmath import *

from psrpbdot import M1
from datatools.tempo import *
from astropy import coordinates as coord
from tools.Coordinate import *
import numpy.linalg as linalg
from GalacticGeometry import *

#from Arrow3D import Arrow3D
secperday = 24*3600
#solardist = 8.34 # old
solardist = 8.34 # Reid et al. 2014, (rmb+14)
G = 6.673e-8
c = 2.99792458e10
PI = np.pi
AU = 1.469e13
Msun = 1.9882e33
Tsun = 4.925490947e-6
secperday = 24*3600
R0 = solardist

"""load the parfile"""
#pf = PARfile('1713.Dec.mcmc.par')
#pf = PARfile('mcmcresult.par')
#pf = PARfile('1713.sns.par')
#pf = PARfile('1713.final.par')
pf = PARfile('Feb.T1.RN.par')

"""read some information from the parfile"""
gl, gb = getGpos(pf)
D = 1./float(pf.PX[0])
Eerr = float(pf.E[1])
OMerr = float(pf.OM[1])/180.*np.pi
"""Calculate the coordinate transfermation matrix"""
T, GT = GetTransferMatrix(pf)
"""Galactic acceleration for low z """
#Kz = lambda z:(2.27*z + 3.68*(1-np.exp(-4.31*z)) ) * 1.e-9 #Galactic acceleration in z direction (cm/s^2)
def Kz(z):
    absz = np.abs(z)
    sign = z/absz
    return sign * (2.27*absz + 3.68*(1-np.exp(-4.31*absz)) ) * 1.e-9


def getKG(kr, zeta, z, sini, paascnode, om):
    """calculate the projection of Galactic acceleration on the orbital plane
       return  "angle between g and periastron", "projected acceleration (km/s)"
    """
    g_r = GT.I * kr * ( np.matrix((np.cos(0.-zeta),np.sin(0.-zeta),0)).T) #X/linalg.norm(X)
    g_z = GT.I * Kz(z) * (np.matrix((0., 0., -1.)).T) 
    g = g_r + g_z
    g_NSEW = T * g

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
    g_ang = np.arccos(A_ref * g_dir)
    g_norm = linalg.norm(g)
    return float(g_ang - om/180.*np.pi), KG

def EccArea(ECC, EF, THETA):
    """ calculate the "area" in the four 3-sigma points for given ECC_observed, ECC_forced and the angle between them.
    """
    global Eerr, OMerr
    #THETA[THETA<0]+=(np.pi*2)
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
            value = max((np.log10(min(max(abs4), 0.05)) - np.log10(max(min(abs4),1.e-6))),0)*np.pi*2
        else:
            value = max((np.log10(min(max(abs4), 0.05)) - np.log10(max(min(abs4),1.e-6))),0)*(max(p4)-min(p4))
        areas.append(value)
    areas = np.array(areas)
    #print 'EF, ECC, THETA, AREA, A4, p4', EF.mean(), ECC.mean(), THETA.mean(), areas.mean(), max(abs4),min(abs4), max(p4), min(p4)
    #print 'Eerr, OMerr', Eerr,OMerr
    return areas

def Pintegrant(PX, SINI, PAASCNODE, M1, M2, PB, ECC, OM, Delta):
    """ Calculate the "probability" for given Delta and timing parameters
    """
    D = 1./PX
    Omega = PAASCNODE/180.*np.pi
    z = np.sin(gb) * D
    R1 = np.sqrt(R0**2 + (D*np.cos(gb))**2 -2 * R0 * D * np.cos(gb) * np.cos(gl))
    coszeta = (R0**2 + R1**2 - D**2 + z**2)/R0/R1/2.
    zeta = np.arccos(coszeta)
    #print R1.size, zeta.size, z.size

    Mtot =  M1 + M2
    #Omega_G = 27.2 #km s^-1 kpc^-1
    Omega_G = 30.57 #km s^-1 kpc^-1 +/- 5.1 Ref: Reid et al. 2014 (rmb+14)
    kpcinkm = 3.0857e16
    Kr =  Omega_G**2 * R0**2 / R1 / kpcinkm * 1.e5 #Galactic acceleration in radio direction (cm/s^2)
    KGarray = []
    THETA = []
    for i,sini in enumerate(SINI):
        theta, kg= getKG(Kr[i], zeta[i], z[i], sini, PAASCNODE[i], OM[i])
        #theta, kg= getKG(Kr, zeta, z, sini, PAASCNODE[i], OM[i])
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
    Areas = 0.
    Areas += sum(EccArea(ECC, EF, THETA))
    Areas += sum(EccArea(ECC, EF, THETA-np.pi))
    return Areas

''' load in the MCMC results for Delta estimation'''
import cPickle as pickle
import sys
secperday = 24*3600

#dic = pickle.load(open('bestpar.p', 'rb'))
#best = dic['BEST']
#plist = dic['parameters'] 
plist = open('pars.txt', 'r').readlines()
pl = [p.strip('\n') for p in plist]
plist = pl
#MChain = pickle.load(open('TinyMChain.p','rb'))
#MChain = pickle.load(open('SmallMChain.p','rb'))
#MarkovChain = MChain['Chain']
MarkovChain = np.loadtxt('chain_1.0.txt')
MCMCSize, cols = MarkovChain.shape
#MChain = pickle.load(open('MChain.p','rb'))
#MarkovChain = MChain['Chain']
#MCMCSize = len(MarkovChain)
pi = 3.141592653589793
twopi = 6.283185307179586
G = 6.673e-8
c = 2.99792458e10
Msun = 1.9882e33
im2 = plist.index('M2')
ipb = plist.index('PB')
#isini = plist.index('SINI')
ii = plist.index('KIN')
ia = plist.index('A1')
iecc = plist.index('ECC')
ipx = plist.index('PX')
#iomega = plist.index('PAASCNODE')
iomega = plist.index('KOM')
iom = plist.index('OM')
#ichisq = plist.index('chisq')
#M2 = np.array([float(p[im2])*Msun for p in MarkovChain])
#PB = np.array([float(p[ipb])*secperday for p in MarkovChain])
#SINI = np.sin(np.array([float(p[ii]) for p in MarkovChain]))
#OM = np.array([float(p[iom]) for p in MarkovChain])
#a = np.array([float(p[ia])*c for p in MarkovChain])
#M1 = (PB/2/pi*np.sqrt(G*(M2*SINI)**3/a**3)-M2)/Msun
#M2 = M2/Msun
#ECC = np.array([float(p[iecc]) for p in MarkovChain])
#PX = np.array([float(p[ipx]) for p in MarkovChain])
#PAASCNODE = np.array([float(p[iomega]) for p in MarkovChain])
#chisq = [p[ichisq] for p in MarkovChain]
#bestidx = chisq.index(min(chisq))
M2 = MarkovChain[:,im2]
PB = secperday * MarkovChain[:,ipb]
SINI = np.sin(MarkovChain[:,ii]/180.*np.pi)
OM = MarkovChain[:,iom]
a =  MarkovChain[:, ia]
#M1 = (PB/2/pi*np.sqrt(G*(M2*SINI)**3/a**3)-M2)#/Msun
M1 = PB/2/pi*(np.sqrt(Tsun*(M2*SINI)**3/a**3))-M2
M2 = M2#/Msun
ECC = MarkovChain[:,iecc]
PX = MarkovChain[:,ipx]
#PX = float(pf.PX)
PAASCNODE = MarkovChain[:,iomega]

#print sini.mean()
#print M1.mean()
#print M2.mean()
#sys.exit(0)

"""Calculate and plot the pdf """
Integ = lambda d: Pintegrant(PX, SINI, PAASCNODE, M1, M2, PB, ECC, OM, d)
delta = np.arange(5.e-5, 0.03, 5.e-5) #ingrid setting
res = np.array([Integ(d) for d in delta])

np.save(open('SEPjae.npy', 'wb'), res)

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

