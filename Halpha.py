# IPython log file
from pylab import *
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

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Verdana']
plt.rcParams['font.size'] = 20
plt.rcParams['lines.linewidth'] = 4.
plt.rcParams['axes.labelsize'] = 'large'
plt.rcParams['grid.linewidth'] = 1.0
plt.rcParams['grid.linestyle'] = ':'
plt.rcParams['xtick.minor.size']=4
plt.rcParams['xtick.major.size']=8
plt.rcParams['figure.figsize'] = 9,9
plt.rcParams['figure.subplot.bottom'] = 0.15

LPS=asciitable.read('LPS_Low_Pressure_Sodium_1797K_LICA.dat')
PC_ambar=asciitable.read('LED_Lytepro7_Architectural.dat')
ambar=asciitable.read('LED_595nm_LICA.dat')
HPS=asciitable.read('HPS_High_Pressure_Sodium_1859K_SFH_Toledo.dat')
HPS2=asciitable.read('HPS_Philips_noHG.dat')

gca()
gca().add_patch(Rectangle((6400,0),600,1,facecolor="yellow"))
plt.text(7000,0.5,r'H$\alpha$ region',ha='right', va='bottom')
gca().add_patch(Rectangle((4800,0),500,1,facecolor="orange"))
plt.text(5300,0.5,r'H$\beta$ region',ha='right', va='bottom')
plt.plot(PC_ambar.col2,PC_ambar.col1/PC_ambar.col1[635],label='PC-Ambar')
plt.plot(LPS.col2,LPS.col1/np.max(LPS.col1),label='LPS')
plt.plot(HPS.col2,HPS.col1/np.max(HPS.col1[915]),label='HPS')
#plt.plot(HPS.col2,(HPS2.col1-HPS2.col1[-1])/np.max(HPS2.col1[915]))
plt.plot(ambar.col2,ambar.col1/ambar.col1[915],label='Ambar')
plt.xlim(4500,7100)
plt.ylim(0,1)
plt.legend()
xlabel(r"Wavelenght($\AA$)")
ylabel("Intensity (relative intensity)")
plt.savefig('Halpha.png')
plt.show()
