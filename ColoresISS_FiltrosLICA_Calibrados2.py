#!/usr/bin/python
# coding: utf-8
# Canales DMSP/OLS y VIIRS
# Espectros de lamparas
# Respuesta de Nikon D3s de LICA
# Colores esperados en ISS pictures

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
JAZ=asciitable.read("JAZ_100.txt")
ISS=asciitable.read('zonas_medidas2.csv',delimiter='\t')

RN=ISS.R1/ISS.Npix
RND=RN-RN[-1]
GN=ISS.G2/ISS.Npix
GND=GN-GN[-1]
BN=ISS.B4/ISS.Npix
BND=BN-BN[-1]
GR=GND/(RND)
BG=BND/GND

JAZ['(05)']=(((10**((JAZ['B_D3s'])/-2.5)))/10**((JAZ['G_D3s']/-2.5)))
JAZ['(07)']=(10**((JAZ['G_D3s'])/-2.5)/(10**((JAZ['R_D3s']/-2.5))))
j=BG>0
i=BG<0
plt.plot(BG[j],GR[j],'ko')


#plt.text(BG,GR,ISS.name,ha='left', va='bottom') 
plt.text(BG[0],GR[0],ISS.Name[0],ha='right', va='bottom')
plt.text(BG[1],GR[1],ISS.Name[1],ha='left', va='bottom')
plt.text(BG[2],GR[2],ISS.Name[2],ha='right', va='bottom')
plt.text(BG[3],GR[3],ISS.Name[3],ha='right', va='bottom')
plt.text(0.03,GR[4],ISS.Name[4],ha='right', va='bottom')
plt.text(BG[5],GR[5]-0.045,ISS.Name[5],ha='right', va='bottom')
plt.text(BG[6],GR[6],ISS.Name[6],ha='right', va='bottom')
plt.text(BG[7],GR[7]-0.04,ISS.Name[7],ha='center', va='baseline')
plt.text(BG[8],GR[8],ISS.Name[8],ha='right', va='bottom')
plt.text(BG[9],GR[9],ISS.Name[9],ha='right', va='bottom')
plt.text(0.03,GR[10],ISS.Name[10],ha='right', va='top')
plt.text(0.03,GR[11],ISS.Name[11],ha='right', va='baseline')
plt.text(0.03,GR[12],ISS.Name[12],ha='right', va='center')
#plt.text(BG[13],GR[13],ISS.Name[13],ha='right', va='bottom')
plt.text(0.05,GR[13]+0.01,ISS.Name[13],ha='left', va='center')
#plt.text(JAZ['(05)'][9],JAZ['(07)'][11],JAZ.name[11],ha='left', va='bottom')

plt.plot(JAZ['(05)'],JAZ['(07)'],'g*',ms=9)
plt.text(JAZ['(05)'][0],JAZ['(07)'][0],JAZ.name[0],ha='left', va='bottom',color='green')
plt.text(JAZ['(05)'][1],JAZ['(07)'][1]+0.02,JAZ.name[1],ha='right', va='bottom',color='green')
plt.text(JAZ['(05)'][2],JAZ['(07)'][2],JAZ.name[2],ha='left', va='bottom',color='green')
plt.text(JAZ['(05)'][3],JAZ['(07)'][3],JAZ.name[3],ha='left', va='bottom',color='green')
plt.text(JAZ['(05)'][4],JAZ['(07)'][4],JAZ.name[4],ha='left', va='bottom',color='green')
plt.text(JAZ['(05)'][5],JAZ['(07)'][5]-0.02,JAZ.name[5],ha='left', va='top',color='green')
plt.text(JAZ['(05)'][6],JAZ['(07)'][6],JAZ.name[6],ha='left', va='bottom',color='green')
plt.text(JAZ['(05)'][7],JAZ['(07)'][7],JAZ.name[7],ha='right', va='bottom',color='green')
plt.text(JAZ['(05)'][8],JAZ['(07)'][8],JAZ.name[8],ha='right', va='bottom',color='green')
plt.text(JAZ['(05)'][9],JAZ['(07)'][9],JAZ.name[9],ha='left', va='bottom',color='green')
plt.text(JAZ['(05)'][10],JAZ['(07)'][10],JAZ.name[10],ha='left', va='bottom',color='green')
plt.text(JAZ['(05)'][11],JAZ['(07)'][11],JAZ.name[11],ha='left', va='bottom',color='green')
plt.text(JAZ['(05)'][12],JAZ['(07)'][12],JAZ.name[12],ha='left', va='bottom',color='green')
plt.text(JAZ['(05)'][13],JAZ['(07)'][13],JAZ.name[13],ha='right', va='bottom',color='green')
plt.text(JAZ['(05)'][14],JAZ['(07)'][14],JAZ.name[14],ha='left', va='bottom',color='green')
plt.text(JAZ['(05)'][15],JAZ['(07)'][15],JAZ.name[15],ha='left', va='bottom',color='green')
plt.text(JAZ['(05)'][16],JAZ['(07)'][16],JAZ.name[16],ha='right', va='top',color='green')
plt.text(JAZ['(05)'][17],JAZ['(07)'][17],JAZ.name[17],ha='left', va='top',color='green')
plt.text(0.26,JAZ['(07)'][19]+0.01,JAZ.name[19],ha='left', va='bottom',color='green')
plt.text(JAZ['(05)'][20],JAZ['(07)'][20],JAZ.name[20],ha='left', va='bottom',color='green')
plt.text(JAZ['(05)'][21],JAZ['(07)'][21],JAZ.name[21],ha='left', va='bottom',color='green')
plt.text(JAZ['(05)'][22],JAZ['(07)'][22],JAZ.name[22],ha='left', va='bottom',color='green')
plt.text(JAZ['(05)'][23],JAZ['(07)'][23],JAZ.name[23],ha='right', va='bottom',color='green')
plt.text(JAZ['(05)'][24],JAZ['(07)'][24],JAZ.name[24],ha='right', va='top',color='green')
plt.text(JAZ['(05)'][25],JAZ['(07)'][25],JAZ.name[25],ha='left', va='top',color='green')

arrowright_verts = [[0.,0.], [2., 0],[0.,0.],[1,1],[0,0],[1,-1]]
scatter(BG[i]/BG[i]*0.035,GR[i], s=100,  color='red',verts=arrowright_verts,marker=None)
plt.xlabel("B/G")
plt.ylabel("G/R")
plt.title("ISS colors for typical lamps")
plt.xlim(0,0.7)
plt.ylim(0,1.4)
plt.savefig('colores_espectros_LaPalma.png')
plt.show()
