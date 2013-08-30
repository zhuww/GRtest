#mass loss term and tidal effects
from math import *
c = 2.99792458e10
PI = pi
AU = 1.469e13
Msun = 1.9882e33
secperday = 24*3600
G = 6.67259e-8 #G in Gaussian unit

Edot_1713 = 3.5e33
Pb_1713  = 67.825168825 * secperday
Mtotal = (1.4+0.3)*Msun

Mdot = Edot_1713/c**2
print 'Pbdot due to Mass loss: ', 2*Mdot/Mtotal *Pb_1713 

k = 0.2
q = 1.4/0.285 #mass ratio
mc = 0.285 #companion mass
sini = sin(72./180*PI)
x = 32.34#*c

Rsun = 69550000000. #solar radius in cm
R_c = 0.018 * Rsun #WD radius (wikipedia, non-degenerate Fermi gas)
Omega_c = sqrt(G*mc*Msun/R_c**3)
print 'Omega_c', Omega_c
ts = 218.81184390161924649/4.0834687472948314223e-16 #tc of the pulsar
print 'ts: ', ts

Pbdot_tidal = k*Omega_c/3./PI/q/(q+1)*(R_c*Pb_1713*sini/x/c)**2 / ts
print 'Pbdot due totidal effect: ', Pbdot_tidal
print 'Pbdot/Pb: ', Pbdot_tidal/Pb_1713
print 'Pbdot/Pb*year: ', Pbdot_tidal/Pb_1713*365*secperday

#Numbers for 1738:
q = 8.1
mc = 0.18
R_c = 0.037 * Rsun
ts = 4.1e9*365*secperday
Pb= 0.35479 * secperday
x = 0.3434
sini = sin(32.6/180*PI)
Omega_c = sqrt(G*mc*Msun/R_c**3)
#Omega_c = 0.038
print 'Omega_c', Omega_c
print 'J1738 ts: ', ts
Pbdot_tidal = k*Omega_c/3./PI/q/(q+1)*(R_c*Pb*sini/x/c)**2 / ts
print 'J1738 Pbdot due totidal effect: ', Pbdot_tidal
