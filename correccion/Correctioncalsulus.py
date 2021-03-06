# IPython log file

import asciitable,sys
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad
from scipy.interpolate import UnivariateSpline
from scipy.interpolate import InterpolatedUnivariateSpline
import os
import re
from sklearn import linear_model, datasets
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure, show, rc, grid, savefig
from matplotlib.ticker import MultipleLocator
import matplotlib.pyplot as plt
import math
import csv
import numpy as np
import sys
from pylab import *
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import asciitable
waveref=np.arange(250,15000,1)
intesity=np.arange(250.0,15000.0,1)*0
LPSDD=asciitable.read('Halo_PAR30_domestic_use_softhwith_color_Sylvania_brand_75W_120V_15_.dat')
UCM=asciitable.read('SPECTRUM0019.dat')
plot(LPSDD.col2,LPSDD.col1)
show()
plot(UCM.col2,UCM.col1)
show()
sLPSDD = InterpolatedUnivariateSpline(LPSDD.col2, LPSDD.col1, k=3)
s_LPSDD=sLPSDD(waveref)
s_LPSDD
plot(waveref,s_LPSDD)
show()
sLPSDD = InterpolatedUnivariateSpline(LPSDD.col1, LPSDD.col2, k=3)
s_LPSDD=sLPSDD(waveref)
plot(waveref,s_LPSDD)
show()
sLPSDD = InterpolatedUnivariateSpline(LPSDD.col2, LPSDD.col1, k=3)
s_LPSDD=sLPSDD(waveref)
plot(waveref,s_LPSDD)
show()
sLPSDD = InterpolatedUnivariateSpline(LPSDD.col2, LPSDD.col1, k=3)
s_LPSDD=sLPSDD(waveref)
plot(waveref,s_LPSDD)
show()
plot(LPSDD.col2,LPSDD.col1)
show()
waveref=np.arange(3800,15000,1)
s_LPSDD=sLPSDD(waveref)
plot(waveref,s_LPSDD)
show()
waveref=np.arange(3800,7300,1)
s_LPSDD=sLPSDD(waveref)
plot(waveref,s_LPSDD)
show()
sUCM = InterpolatedUnivariateSpline(UCM.col2, UCM.col1, k=3)
s_UCM=sUCM(waveref)
plot(waveref,s_UCM)
show()
plot(waveref,s_LPSDD/s_UCM)
show()

correctionJAZ=(s_LPSDD/s_UCM)/np.min(s_LPSDD/s_UCM)
z = np.polyfit(waveref, correctionJAZ, 15)
p = np.poly1d(z)
plot(waveref,(s_LPSDD/s_UCM)/np.min(s_LPSDD/s_UCM))
plot(waveref,p(waveref))
show()
asciitable.write({'y': p(waveref),'x': waveref },'CorrectionJAZ.csv',Writer=asciitable.NoHeader ,names=['y', 'x'])
get_ipython().magic(u'logstart')
exit()
