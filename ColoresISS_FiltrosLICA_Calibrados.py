#!/usr/bin/python
# coding: utf-8
# Canales DMSP/OLS y VIIRS
# Espectros de lamparas
# Respuesta de Nikon D3s de LICA
# Colores esperados en ISS pictures

import matplotlib.pyplot as plt
from matplotlib.pyplot import figure, show, rc, grid, savefig
from matplotlib.ticker import MultipleLocator
import math
import csv
import numpy as np
import sys
from pylab import *
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import asciitable

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
#plt.rcParams['ytick.labelsize'] = 10
JAZ=asciitable.read("JAZ_99.csv")
ISS=asciitable.read('zonas_medidas.csv',delimiter=';')

RN=ISS.R1/ISS.Npix
RND=RN-RN[0]
GN=ISS.G2/ISS.Npix
GND=GN-GN[0]
BN=ISS.B4/ISS.Npix
BND=BN-BN[0]
GR=GND/(RND)
BG=BND/GND

JAZ['(05)']=(((10**((JAZ['B_D3s'])/-2.5))*0.87977)/10**((JAZ['G_D3s']/-2.5)))*0.8800116943392164
JAZ['(07)']=(10**((JAZ['G_D3s'])/-2.5)/(10**((JAZ['R_D3s']/-2.5))*0.768953))*1.4
plt.plot(BG,GR,'ko')
#plt.text(BG,GR,ISS.name,ha='left', va='bottom') 
plt.plot(JAZ['(05)'],JAZ['(07)'],'g*')
plt.text(JAZ['(05)'][0],JAZ['(07)'][0],JAZ.name[0],ha='left', va='bottom')
plt.text(JAZ['(05)'][1],JAZ['(07)'][1],JAZ.name[1],ha='left', va='bottom')
plt.text(JAZ['(05)'][2],JAZ['(07)'][2],JAZ.name[2],ha='left', va='bottom')
plt.text(JAZ['(05)'][3],JAZ['(07)'][3],JAZ.name[3],ha='left', va='bottom')
plt.text(JAZ['(05)'][4],JAZ['(07)'][4],JAZ.name[4],ha='left', va='bottom')
plt.text(JAZ['(05)'][5],JAZ['(07)'][5],JAZ.name[5],ha='left', va='bottom')
plt.text(JAZ['(05)'][6],JAZ['(07)'][6],JAZ.name[6],ha='left', va='bottom')
plt.text(JAZ['(05)'][7],JAZ['(07)'][7],JAZ.name[7],ha='left', va='bottom')
plt.text(JAZ['(05)'][8],JAZ['(07)'][8],JAZ.name[8],ha='left', va='bottom')
plt.text(JAZ['(05)'][9],JAZ['(07)'][9],JAZ.name[9],ha='left', va='bottom')
plt.xlabel("B/G")
plt.ylabel("G/R")
plt.title("ISS colors for typical lamps")
plt.show()
