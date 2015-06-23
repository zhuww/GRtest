from datatools.tempo import *
from pylab import *

md = model('J0023+0923-tempo1-SEPtesting.par')
tf = TOAfile('J0023+0923-tempo1-SEPtesting.tim')

md.tempofit(tf)

md.plot('mjd', 'res')

show()
