import numpy as np
import matplotlib.pyplot as plt
from pylab import rand
from cmath import *

res = np.load('alpha3.npy')
#res = np.load('SEPresult.npy')

delta = np.arange(5.e-22, 1.e-19, 5.e-22) #ingrid setting
#delta = np.arange(5.e-5, 0.03, 5.e-5) #ingrid setting


fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(delta, res, '-')
#ax.logx()
#ax.hist(delta, bins=50, normed=1)
ax.semilogx(nonposy='clip')
ax.semilogy(nonposy='clip')

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

plt.show()
