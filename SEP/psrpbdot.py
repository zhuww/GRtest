from decimal import *
from math import *
from tools.Coordinate import RA, Dec
from astropy import coordinates as coord
from astropy import constants as const

Tsun = Decimal('4.925490947')*Decimal('0.000001') #Tsun == GM/c^3 in seconds
PI = Decimal(pi)
secperday = 24*3600
secperyear = secperday*365.24218967
#AU = Decimal(1.469e13)
#c = Decimal(2.99792458e10)
c = const.c.cgs.value
kpc = const.kpc.cgs.value
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
    fac1 = float(pf.M2[1]/pf.M2[0])
    fac2 = float(pf.A1[1]/pf.A1[0])
    if pf.SINI == 'KIN':
        fac3 = abs(cos(float(pf.KIN[0])*pi/180))*float(pf.KIN[1])
    else:
        fac3 = float(pf.SINI[1]/pf.SINI[0])
    val = float(- 192*PI/5 * (2*PI/Pb)**Decimal(5./3) * (Tsun*M2)**Decimal(5./3) * q / (q+1)**Decimal(1./3))
    err = sqrt(5*fac1**2 + 3* fac2**2 + 3*fac3**2) * abs(val)
    return val, err


def Shlkovskii(pf):
    #d = AU * 180 / PI * 3600 / pf.PX[0] * 1000 
    d = float(1/pf.PX[0]) * kpc
    Pb = float(pf.PB[0]) * secperday
    PMRA = pf.PMRA[0] / 3600000 * PI / 180 / secperday / Decimal(365.24218967)
    PMDEC = pf.PMDEC[0] / 3600000 * PI / 180 / secperday / Decimal(365.24218967) #* Decimal(str(sin(inc)))
    val = float(PMRA**2 + PMDEC**2)*Pb*d/c
    fac1 = float(pf.PX[1]/pf.PX[0])
    fac2 = (sqrt(float(pf.PMRA[1]**2 + pf.PMDEC[1]**2)/float(pf.PMRA[0]**2+pf.PMDEC[0]**2)))
    err = sqrt(fac1**2 + fac2**2)*abs(val)
    return val, err


def Pbdot_Gal(pf):
    try:
        name = pf.PSRJ
    except:
        name = pf.PSR
    ra = RA(pf.RAJ[0])
    dec = Dec(pf.DECJ[0])
    pos = coord.FK5Coordinates(str(ra) +' '+ str(dec))
    l = pos.galactic.l.radians
    b = pos.galactic.b.radians
    Pb = float(pf.PB[0] * secperday)
    #d = AU * 180 / PI * 3600 / pf.PX[0] * 1000 
    d = float(1/pf.PX[0])
    pf.DIST = d
    z_kpc = float(d)*(abs(sin(b)))#/kpc must be positive
    a_z = ((2.27)*z_kpc + (3.68)*(1 - exp((-4.31)*z_kpc)))*(1.e-9) #cm s^-2
    #print 'a_z:', a_z
    A_z = -1 * a_z *abs(sin(b))/c
    pf.A_z = A_z
    R0 = (8.34) #* kpc # Reid et al. 2014
    beta = float(d/R0) * cos(b) - cos(l)
    Omega0 = (240. * 1.e5) #240+/-8 km/s; Reid et al  2014
    #print b,l, cos(b), cos(l), beta
    A_x = -1/c * (cos(b)) * (Omega0**2/R0/kpc) * (cos(l) + beta/(sin(l)**2 + beta**2))
    pf.A_x = A_x
    #print 'Ax, Az: ',A_x, A_z
    fac1 = float(pf.PX[1]/pf.PX[0])
    fac2 = 8/240 #Omega_0
    fac3 = 0.16/8.34 #R0 Reid et al 2014
    val = float(Pb*(A_z + A_x))
    err1 = A_x * fac1 
    err2 = A_z * sqrt(fac1**2 + fac2**2 + fac3**2)
    err = sqrt(err1**2 + err2**2) * float(Pb)
    return val, err



if __name__ == '__main__':
    from datatools.tempo import PARfile
    #pf  = PARfile('1909-3744.par')
    #pf  = PARfile('./J1909.zww.par')
    #pf  = PARfile('./J1713+0747.par')
    #pf  = PARfile('./1713.Sep.par')
    #pf  = PARfile('./mcmcresult.par')
    #pf  = PARfile('./1713.Apr.par')
    pf  = PARfile('./Oct.T2.par')
    if pf.PBDOT[0] > Decimal('1.e-10'):#tempo 1 units
        pf.PBDOT[0] *= Decimal('1.e-12')
        pf.PBDOT[1] *= Decimal('1.e-12')

    from round import shortform as SF
    Shl = Shlkovskii(pf)
    Gal = Pbdot_Gal(pf)
    GW = Pbdot_GW(pf)
    if pf.__dict__.has_key('PSR'):
        print pf.PSR
    else:
        print pf.PSRJ
    print 'PBdot(shlkovskii):', SF(Shl)
    print 'PBdot(Galaxtic):' , SF(Gal)
    print 'PBdot(GW):', SF(GW)
    pb_exc = (float(pf.PBDOT[0]) - Gal[0] - Shl[0] - GW[0]), sqrt(float(pf.PBDOT[1])**2 + Gal[1]**2 + Shl[1]**2 + GW[1]**2)
    print "PBdot - PBdot(Galaxtic) - PBdot(shlkovskii) - PBdot(GW), PBdot_exc_err:",
    print SF(pb_exc)   
    print float(pf.PBDOT[0]) - Gal[0] - Shl[0] - GW[0], sqrt(float(pf.PBDOT[1])**2 + Gal[1]**2 + Shl[1]**2 + GW[1]**2)
    print "PBdot , PBdot(Galaxtic) + PBdot(shlkovskii) + PBdot(GW):",
    print pf.PBDOT[0] , Shl[0]+Gal[0]+GW[0]
    print "-PBdot_excess/PB/2 --> Gdot/G (95%):",

    print SF((-1*pb_exc[0]/float(pf.PB[0])/secperday/2*secperyear, pb_exc[1]/float(pf.PB[0])/secperday*secperyear))
    print ((pb_exc[0]/float(pf.PB[0])/secperday*secperyear, pb_exc[1]/float(pf.PB[0])/secperday*secperyear*2))
    #print SF(( float((pf.PBDOT[0] - Pbdot_Gal(pf) - Decimal(Shlkovskii(pf)) )/ pf.PB[0]/ secperday) , float( pf.PBDOT[1]*2 / pf.PB[0]/ secperday ) ))

    #Pbdot_exc_1713 = (pf.PBDOT[0] - Pbdot_Gal(pf) - Decimal(Shlkovskii(pf)) )
