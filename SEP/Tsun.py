import numpy as np
from decimal import Decimal


G = Decimal(6.67259e-8)
c = Decimal(2.99792458e10)
PI = Decimal(np.pi)
AU = Decimal(1.469e13)
Msun = Decimal(1.9891e33)
secperday = 24*3600
solardist = 8.33

print G*Msun/c**3 * 1000000
