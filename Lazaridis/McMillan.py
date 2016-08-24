from decimal import *
from math import *
import numpy as np
from commands import getoutput
from astropy import coordinates as coord
from astropy import constants as const
from tools.Coordinate import RA, Dec
from round import TexStyle as SF

Tsun = Decimal('4.925490947')*Decimal('0.000001') #Tsun == GM/c^3 in seconds
PI = Decimal(np.pi)
secperday = 24*3600
secperyear = secperday*365.24218967
#AU = Decimal(1.469e13)
#c = Decimal(2.99792458e10)
c = const.c.cgs.value
kpc = const.kpc.cgs.value
#Msun = Decimal(0.9882e33)

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
    print 'L, B, D:', L, B, d
    pf.DIST = d
    z_kpc = float(d)*(abs(sin(b)))#/kpc must be positive
    a_z = ((2.27)*z_kpc + (3.68)*(1 - exp((-4.31)*z_kpc)))*(1.e-9) #cm s^-2
    A_z = -1 * a_z *abs(sin(b))/c
    pf.A_z = A_z
    R0 = (8.34) #* kpc # Reid et al. 2014
    beta = float(d/R0) * cos(b) - cos(l)
    Omega0 = (240. * 1.e5) #240+/-8 km/s; Reid et al  2014
    A_x = -1/c * (cos(b)) * (Omega0**2/R0/kpc) * (cos(l) + beta/(sin(l)**2 + beta**2))
    pf.A_x = A_x
    fac1 = float(pf.PX[1]/pf.PX[0])
    fac2 = 8/240 #Omega_0
    fac3 = 0.16/8.34 #R0 Reid et al 2014
    val = float(Pb*(A_z + A_x))
    err1 = A_x * fac1 
    err2 = A_z * sqrt(fac1**2 + fac2**2 + fac3**2)
    err = sqrt(err1**2 + err2**2) * float(Pb)
    return val, err


if __name__ == '__main__':
    from datatools.tempo import tempofit, tempo2fit, touchparfile, uniquename, PARfile
    pf = PARfile('1713.cut.par')
    print Pbdot_Gal(pf)
