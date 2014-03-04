from decimal import *
from math import *
from tools.Coordinate import RA, Dec
from astropy import coordinates as coord
from astropy import constants as const

Tsun = Decimal('4.925490947')*Decimal('0.000001') #Tsun == GM/c^3 in seconds
PI = Decimal(pi)
secperday = 24*3600
AU = Decimal(1.469e13)
#c = Decimal(2.99792458e10)
c = Decimal(const.c.cgs.value)
kpc = Decimal(const.kpc.cgs.value)
#Msun = Decimal(0.9882e33)


def M1(pf):
    Tsun = Decimal('4.925490947')*Decimal('0.000001') #Tsun == GM/c^3 in seconds
    m2 = pf.M2[0]
    Pb = pf.PB[0]*secperday
    a = pf.A1[0]
    if pf.__dict__.has_key('KIN'):
        I = pf.KIN[0]/180*PI
        sini = Decimal(sin(float(I)))
    else:
        sini = pf.SINI[0]
    #result = sqrt(930.998*m2**3*Pb**2/a**3) - m2
    #return (Pb/2/PI*Decimal(sqrt(G*(m2*Decimal(str(sini)))**3/a**3))-m2)/Msun
    return Pb/2/PI*((Tsun*(m2*sini)**3/a**3)**Decimal(0.5))-m2


def Pbdot_GW(pf):
    Pb = pf.PB[0] * secperday
    M2 = pf.M2[0]
    q = M1(pf)/M2
    return - 192*PI/5 * (2*PI/Pb)**Decimal(5./3) * (Tsun*M2)**Decimal(5./3) * q / (q+1)**Decimal(1./3)


def Shlkovskii(pf):
    d = AU * 180 / PI * 3600 / pf.PX[0] * 1000 
    Pb = pf.PB[0] * secperday
    PMRA = pf.PMRA[0] / 3600000 * PI / 180 / secperday / Decimal(365.24218967)
    PMDEC = pf.PMDEC[0] / 3600000 * PI / 180 / secperday / Decimal(365.24218967) #* Decimal(str(sin(inc)))
    return float((PMRA**2 + PMDEC**2)*Pb*d/c)


def Pbdot_Gal(pf):
    name = pf.PSRJ
    ra = RA(pf.RAJ[0])
    dec = Dec(pf.DECJ[0])
    pos = coord.FK5Coordinates(str(ra) +' '+ str(dec))
    l = pos.galactic.l.radians
    b = pos.galactic.b.radians
    Pb = pf.PB[0] * secperday
    #d = AU * 180 / PI * 3600 / pf.PX[0] * 1000 
    d = 1/pf.PX[0]
    z_kpc = d*Decimal(sin(b))#/kpc
    a_z = (Decimal(2.27)*z_kpc + Decimal(3.68)*Decimal(1 - exp(Decimal(-4.31)*z_kpc)))*Decimal(1.e-9) #cm s^-2
    A_z = -1 * a_z *abs(Decimal(sin(b)))/c
    R0 = Decimal(8.34) #* kpc # Reid et al. 2014
    Omega0 = Decimal(240. * 1.e5) #240+/-8 km/s; Reid et al  2014
    beta = float(d/R0) * cos(b) - cos(l)
    A_x = -1/c * Decimal(cos(b)) * (Omega0**2/R0/kpc) * Decimal(cos(l) + beta/(sin(l)**2 + beta**2))
    return Pb*(A_z + A_x)


if __name__ == '__main__':
    from datatools.tempo import PARfile
    pf  = PARfile('./J1909.zww.par')
    #pf  = PARfile('./mcmcresult.par')
    pf.PBDOT[0] = pf.PBDOT[0]*Decimal('1.e-12')
    pf.PBDOT[1] = pf.PBDOT[1]*Decimal('1.e-12')

    print 'PBdot(shlkovskii):', Shlkovskii(pf)
    print 'PBdot(Galaxtic):' , Pbdot_Gal(pf)
    print 'PBdot(GW):', Pbdot_GW(pf)
    print "PBdot - PBdot(Galaxtic) - PBdot(shlkovskii) - PBdot(GW), PBdoterr*2:",
    print (pf.PBDOT[0] - Pbdot_Gal(pf) - Decimal(Shlkovskii(pf)) - Pbdot_GW(pf)), pf.PBDOT[1]*2
    print "PBdot , PBdot(Galaxtic) + PBdot(shlkovskii) + PBdot(GW):",
    print pf.PBDOT[0] , Pbdot_Gal(pf) + Decimal(Shlkovskii(pf)) +  Pbdot_GW(pf)
    print "PBdot_excess/PB --> Gdot/G:",
    from round import shortform as SF
    print SF(( float((pf.PBDOT[0] - Pbdot_Gal(pf) - Decimal(Shlkovskii(pf)) )/ pf.PB[0]/ secperday) , float( pf.PBDOT[1]*2 / pf.PB[0]/ secperday ) ))

    Pbdot_exc_1713 = (pf.PBDOT[0] - Pbdot_Gal(pf) - Decimal(Shlkovskii(pf)) )
