from decimal import *
from math import *
import numpy as np
from commands import getoutput
from astropy import coordinates as coord
from astropy import constants as const
from tools.Coordinate import RA, Dec
from round import TexStyle as SF
from threadit import spamit

Tsun = Decimal('4.925490947')*Decimal('0.000001') #Tsun == GM/c^3 in seconds
PI = Decimal(np.pi)
secperday = 24*3600
secperyear = secperday*365.24218967
#AU = Decimal(1.469e13)
#c = Decimal(2.99792458e10)
c = const.c.cgs.value
kpc = const.kpc.cgs.value
#Msun = Decimal(0.9882e33)

def McMillanPot(L, B, d):
    return float(getoutput('./GalPotMcMillan2016/calcGalPdot.exe %s %s %s' % (L, B, d)))

def SurfaceDensity(Mod, R):
    return float(getoutput('./GalPotMcMillan2016/surfacedensity.exe %s %s' % (Mod, R)))/1.e6 #Msun/pc^2, while the original unit is Msun/kpc^2


def dPoP(L, B, d, Mod):
    res = getoutput('./GalPotMcMillan2016/compPdotMods.exe %s %s %s %s' % (L, B, d, Mod))
    hor, ver= res.split(' ')
    return float(hor), float(ver)

def Pbdot_Gal(pf):
    try:
        name = pf.PSRJ
    except:
        name = pf.PSR
    ra = RA(pf.RAJ[0])
    dec = Dec(pf.DECJ[0])
    pos = coord.SkyCoord(str(ra) +' '+ str(dec))
    l = pos.galactic.l.rad
    b = pos.galactic.b.rad
    L = pos.galactic.l.deg
    B = pos.galactic.b.deg

    Pb = float(pf.PB[0] * secperday)
    d = float(1/pf.PX[0])
    pf.DIST = d
    #print 'L, B, D:', L, B, d
    #output = getoutput('./GalPotMcMillan2016/calcGalPdot.exe %s %s %s' % (L, B, d))
    #print output
    px = np.random.normal(float(pf.PX[0]), float(pf.PX[1]), 10000)
    res = np.array(spamit(McMillanPot, [(L, B, 1/x) for x in px]))
    
    val = res.mean() * Pb
    fac2 = 8/240 #Omega_0
    fac3 = 0.16/8.34 #R0 Reid et al 2014
    fac1 = res.std()/res.mean()
    err = abs(sqrt(fac1**2 + fac2**2 + fac3**2) * val)
    return val, err


if __name__ == '__main__':
    from datatools.tempo import tempofit, tempo2fit, touchparfile, uniquename, PARfile
    pf = PARfile('1713.cut.par')
    print Pbdot_Gal(pf)
