#!/usr/bin/env python

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

def ransacfit(X,y):
    X=X.reshape(len(X),1)
    model_ransac = linear_model.RANSACRegressor(linear_model.LinearRegression())
    model_ransac.fit(X, y)
    inlier_mask = model_ransac.inlier_mask_
    outlier_mask = np.logical_not(inlier_mask)
    model = linear_model.LinearRegression()
    model.fit(X, y)
    line_X = np.arange(-5, 5)
    line_y = model.predict(line_X[:, np.newaxis])
    line_y_ransac = model_ransac.predict(line_X[:, np.newaxis])
    return model_ransac,line_y_ransac,line_X

def normaliza(X): #Normalize to the maximum
    max=np.max(X.intrel)
    X.intrel=X.intrel/max
    return X

def srtoarcsec2(X): #Change of units
    X=X* 4.254517e+10
    return X

def srtoarcsec2(X): #Change of units
    X=X* 2.3504431e-11
    return X

def ergstoW(X): #Change of units
      X=X* 1e-07
      return X
def Wtoergs(X): #Change of units
      X=X/1e-07
      return X

def cm_2tom_2(X): #Change of units
      X=X*10000
      return X

def m_2tocm_2(X): #Change of units
      X=X/10000
      return X

def Angstrontonm(X): #Change of units
      X=X/10
      return X

def nmtoAngstron(X): #Change of units
      X=X*10
      return X

def ABtoFlux(AB,filtro): #AB magnitudes to flux #[W cm -2 arsec-2 amstrong-1]
    w=np.sum(filtro.intrel)/float(filtro.size)
    Flux=10^((AB-5*np.log10(w)-2.401)/-2.5)
    return Flux
def FluxtoAB(Flux,filtro): #Flux [W cm -2 arsec-2 amstrong-1] to AB
    w=np.sum(filtro.intrel)/float(filtro.size)
    AB = -2.5*np.log10(Flux)-5*np.log10(w)+2.401 #[W cm -2 arsec-2 amstrong-1]
    return AB


#
#
#   Candela definition: http://physics.nist.gov/cuu/Units/candela.html
#
#
####################################
c=299792458.0
frec=540*10**12
lamnda=c/frec
Inten= 1.0/683#[W.str]
waveref=np.arange(250,15000,1)
intesity=np.arange(250.0,15000.0,1)*0
ncd=list(waveref).index(int(lamnda*1E10))
intesity[ncd]=1.0/683

#intensityAB=Wtoergs(X)
cd=np.recarray((len(waveref),),dtype=[('wave', float), ('flux', float)])
cd.wave=np.array(waveref)
cd.flux=np.array(intesity)

####################################


####################
#
# AB magnitude definition #erg s-1 cm-2 Amstrong-1 
#
# The AB magnitude system is defined as (Oke and Gunn 1983)
# where F_nu is the flux in erg cm-2 s-1 Hz-1
#
# AB = V +AB_nu = - 2.5 *logF_nu  - 48.60
#
#   Here it is used in #erg s-1 cm-2 Amstrong-1 
#
####################
waveref=np.arange(250,15000,1)
NPMAX=10000   
WL_VEGA_AB=[]
FLUX_VEGA_AB=[]
for X in range(NPMAX):
    WL_VEGA1=1000.0+(float(X)/(NPMAX-1)*(30000.0-1000.0))
    WL_VEGA_AB.append(WL_VEGA1)
    FLUX_VEGA1=0.1088/(WL_VEGA1*WL_VEGA1)
    FLUX_VEGA_AB.append(FLUX_VEGA1)




AB=np.recarray((len(FLUX_VEGA_AB),),dtype=[('wave', float), ('flux', float)])
AB.wave=np.array(WL_VEGA_AB)
AB.flux=np.array(FLUX_VEGA_AB) #erg s-1 cm-2 Amstrong-1

####################
#
# AB magnitude definition END
#
####################

####################
#
# ST magnitude definition 
#
####################

WL_VEGA_ST=[]
FLUX_VEGA_ST=[]
for X in range(NPMAX):
    WL_VEGA1=1000.0+(float(X)/(NPMAX-1)*(30000.0-1000.0))
    WL_VEGA_ST.append(WL_VEGA1)
    FLUX_VEGA1=3.63E-9
    FLUX_VEGA_ST.append(FLUX_VEGA1)

####################
#
# ST magnitude definition END
#
####################

sun_lambda=asciitable.read('SUN_CIE.txt') #Cie D65 Sun spectrum

#filename='outfile.dat'
#filename=sys.argv[1]
#spectro=asciitable.read(filename)

####################
#
# Filters Start
#
####################
B_D3=asciitable.read("B_D3.csv") #Nikon D3
G_D3=asciitable.read("G_D3.csv")
R_D3=asciitable.read("R_D3.csv")
B_D3s=asciitable.read("B_D3s.csv") #Nikon D3s
G_D3s=asciitable.read("G_D3s.csv")
R_D3s=asciitable.read("R_D3s.csv")
U_johnson=asciitable.read("U_johnson.csv") #Filters
B_johnson=asciitable.read("B_johnson.csv")
V_johnson=asciitable.read("V_johnson.csv")
R_johnson=asciitable.read("R_johnson.csv")
msas=asciitable.read("msas.csv")
scoto=asciitable.read("scoto.csv")
photo=asciitable.read("photo.csv")
psas=asciitable.read("psas.csv")
SQM_LICA=asciitable.read("SQM_LICA.csv")
VIIRS_2013=asciitable.read("VIIRS_2013.csv")
DMSP_1999=asciitable.read("DMSP_1999.csv")
Halpha=asciitable.read("Halpha.csv")
Hbeta=asciitable.read("Hbeta.csv")

####################
#
# Filters end
#
####################

##########################
#
#
# Fiting planes from: http://stackoverflow.com/questions/15959411/fit-points-to-a-plane-algorithms-how-to-iterpret-results
#
###############################33


import numpy as np
import scipy.optimize

def fitPLaneLTSQ(XYZ):
    # Fits a plane to a point cloud, 
    # Where Z = aX + bY + c        ----Eqn #1
    # Rearanging Eqn1: aX + bY -Z +c =0
    # Gives normal (a,b,-1)
    # Normal = (a,b,-1)
    [rows,cols] = XYZ.shape
    G = np.ones((rows,3))
    G[:,0] = XYZ[:,0]  #X
    G[:,1] = XYZ[:,1]  #Y
    Z = XYZ[:,2]
    (a,b,c),resid,rank,s = np.linalg.lstsq(G,Z) 
    normal = (a,b,-1)
    nn = np.linalg.norm(normal)
    normal = normal / nn
    return normal


def fitPlaneSVD(XYZ):
    [rows,cols] = XYZ.shape
    # Set up constraint equations of the form  AB = 0,
    # where B is a column vector of the plane coefficients
    # in the form b(1)*X + b(2)*Y +b(3)*Z + b(4) = 0.
    p = (np.ones((rows,1)))
    AB = np.hstack([XYZ,p])
    [u, d, v] = np.linalg.svd(AB,0)        
    B = v[3,:];                    # Solution is last column of v.
    nn = np.linalg.norm(B[0:3])
    B = B / nn
    return B[0:3]


def fitPlaneEigen(XYZ):
    # Works, in this case but don't understand!
    average=sum(XYZ)/XYZ.shape[0]
    covariant=np.cov(XYZ - average)
    eigenvalues,eigenvectors = np.linalg.eig(covariant)
    want_max = eigenvectors[:,eigenvalues.argmax()]
    (c,a,b) = want_max[3:6]    # Do not understand! Why 3:6? Why (c,a,b)?
    normal = np.array([a,b,c])
    nn = np.linalg.norm(normal)
    return normal / nn  

def fitPlaneSolve(XYZ):
    X = XYZ[:,0]
    Y = XYZ[:,1]
    Z = XYZ[:,2] 
    npts = len(X)
    A = np.array([ [sum(X*X), sum(X*Y), sum(X)],
                   [sum(X*Y), sum(Y*Y), sum(Y)],
                   [sum(X),   sum(Y), npts] ])
    B = np.array([ [sum(X*Z), sum(Y*Z), sum(Z)] ])
    normal = np.linalg.solve(A,B.T)
    nn = np.linalg.norm(normal)
    normal = normal / nn
    return normal.ravel()

def fitPlaneOptimize(XYZ):
    def residiuals(parameter,f,x,y):
        return [(f[i] - model(parameter,x[i],y[i])) for i in range(len(f))]


    def model(parameter, x, y):
        a, b, c = parameter
        return a*x + b*y + c

    X = XYZ[:,0]
    Y = XYZ[:,1]
    Z = XYZ[:,2]
    p0 = [1., 1.,1.] # initial guess
    result = scipy.optimize.leastsq(residiuals, p0, args=(Z,X,Y))[0]
    normal = result[0:3]
    nn = np.linalg.norm(normal)
    normal = normal / nn
    return normal


#####################################3
#
# Fiting planes End
#
#######################333
def normaliza(X):
    max=np.max(X.intrel)
    X.intrel=X.intrel/max
    return X

#AB_ref = np.array(WL_VEGA_AB, dtype=[('wave', float)])
def mag(filtroX,espectro,referencia): #Calculate the sintetic magnitude and also the Intensity.
    if 0:
     if filtroX.wave[-1]>espectro.col2[-1]:
        dif=filtroX.wave[-1]-espectro.col2[-1]
        n=espectro.size
        new=n+int(dif)+2
        espectro=np.resize(espectro,new)
        espectro[n:-1]['col2']=np.linspace(espectro['col2'][n-1],filtroX.wave[-1],num=espectro[n:-1]['col2'].size)  
    dark=0.0
    if 0:
     if sum(espectro.col1[-10:-1])>1:
       dark=(sum(espectro.col1[-10:-1])/10)
    espectro.col1=espectro.col1-dark
    filtro1=espectro.col1<0
    #espectro.col1[filtro1]=0
    s = InterpolatedUnivariateSpline(espectro['col2'], espectro['col1'], k=3)
    AB_r =  InterpolatedUnivariateSpline(referencia.wave, referencia.flux, k=3)
    filtro = InterpolatedUnivariateSpline(filtroX.wave, filtroX.intrel,k=3)

    xs = np.arange(filtroX.wave[0], filtroX.wave[-1], 1)
    s_1 = s(xs)
    f_1 = filtro(xs)
    AB_1= AB_r(xs)
    S1=s_1*f_1
    S2=AB_1*f_1
    S1_1=InterpolatedUnivariateSpline(xs,S1,k=3)
    S1_2=InterpolatedUnivariateSpline(xs,S2,k=3)
    Up = quad(S1_1, xs[0], xs[-1])
    Dn = quad(S1_2, xs[0], xs[-1])
    m=-2.5*np.log10(Up[0]/Dn[0])
    I=10**(m/-2.5)
    return m,I

def IndiceIsoLux(filtroX,espectro,referencia,photo): #Calculate Index as is definide in http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0067798
    s = InterpolatedUnivariateSpline(espectro.col2, espectro.col1, k=3)
    ref =  InterpolatedUnivariateSpline(referencia.wave, referencia.flux, k=3)
    filtro = InterpolatedUnivariateSpline(filtroX.wave, filtroX.intrel,k=3)
    photo2 = InterpolatedUnivariateSpline(photo.wave, photo.intrel,k=3)
    xs = np.arange(photo.wave[0], photo.wave[-1], 50)
    s_1 = s(xs)
    f_1 = filtro(xs)
    AB_1= ref(xs)
    ph_1=photo2(xs)
    S1n=s_1*ph_1
    AB_1n=AB_1*ph_1
    NormS=np.sum(S1n)
    NormAB=np.sum(AB_1n)
    s_1x=s_1/NormS
    AB_1x=AB_1/NormAB
    s_1f=np.sum(s_1x*f_1)
    AB_1f=np.sum(AB_1x*f_1)
    result=np.sum(s_1f/AB_1f)

    
    return result




def calculamagnitudes(filename):
    spectro=asciitable.read(filename)
    m_Uj,I_Uj=mag(U_johnson,spectro,AB)
    m_Bj,I_Bj=mag(B_johnson,spectro,AB)
    m_Vj,I_Vj=mag(V_johnson,spectro,AB)
    m_Rj,I_Vj=mag(R_johnson,spectro,AB)
    m_Bn,I_Bn=mag(B_D3,spectro,AB)
    m_Gn,I_Gn=mag(G_D3,spectro,AB)
    m_Rn,I_Rn=mag(R_D3,spectro,AB)
    m_SQM,I_Rn=mag(SQM_LICA,spectro,AB)
    m_VIIRS,I_VIIRS=mag(VIIRS_2013,spectro,AB)
    m_DMSP,I_DMSP=mag(DMSP_1999,spectro,AB)
    m_photo,I_photo=mag(photo,spectro,AB)
    m_scoto,I_scoto=mag(scoto,spectro,AB)
    m_Halpha,I_Halpha=mag(Halpha,spectro,AB)
    m_Hbeta,I_Hbeta=mag(Hbeta,spectro,AB)
    Uj_Vj=m_Uj-m_Vj
    Bj_Vj=m_Bj-m_Vj
    Rj_Vj=m_Rj-m_Vj
    Bn_Vj=m_Bn-m_Vj
    Gn_Vj=m_Gn-m_Vj
    Rn_Vj=m_Rn-m_Vj
    SQM_Vj=m_SQM-m_Vj
    Bn_Gn=I_Bn/I_Gn
    Gn_Rn=I_Gn/I_Rn
    i_sli=IndiceIsoLux(scoto,spectro,sun_lambda,photo)
    i_msi=IndiceIsoLux(msas,spectro,sun_lambda,photo)
    i_ipi=IndiceIsoLux(psas,spectro,sun_lambda,photo)
    if 0:
        fig, (ax0,ax1,ax2) = plt.subplots(nrows=3)

        ax0.plot(WL_VEGA_AB, FLUX_VEGA_AB)
        ax0.plot(WL_VEGA_ST, FLUX_VEGA_ST)
        ax0.set_title(filename)
        ax1.plot(spectro.col2, spectro.col1)
        ax1.set_title('Spectra')
        ax2.plot(R_D3s.wave,R_D3s.intrel)
        ax2.plot(G_D3s.wave,G_D3s.intrel)
        ax2.plot(B_D3s.wave,B_D3s.intrel)
        #ax2.plot(photo.wave,photo.intrel)
        ax2.plot(psas.wave,psas.intrel)
        #ax2.plot(scoto.wave,scoto.intrel)
        #ax2.plot(msas.wave,msas.intrel)
        ax2.set_title('Filter')
        plt.show()
    return Uj_Vj,Bj_Vj,Bj_Vj,Rj_Vj,SQM_Vj,m_Bn,m_Gn,m_Rn,i_sli,i_msi,i_ipi,Bn_Vj,Gn_Vj,Rn_Vj,m_VIIRS,m_photo,m_scoto,m_SQM,m_Halpha,m_Hbeta,m_Vj

def filtrog(other,g):
    a=np.where(other<g)
    for j in range(len(np.where(other<g)[0])):
                other[a[0][j],a[1][j]]=np.nan
    return other

def filtrom(other,m):
    a=np.where(other>m)
    for j in range(len(np.where(other>m)[0])):
                other[a[0][j],a[1][j]]=np.nan
    return other

path=os.getcwd()
lista=os.listdir(path)
lista2=" ".join(lista)
lista3=re.findall(r'[a-zA-Z_0-9]+\.dat',lista2)
lista3.sort()
listafiles=lista3






spec=np.recarray((len(listafiles),),dtype=[('num', int),('name', 'S100'),('Uj_Vj', float), ('Bj_Vj', float),('Rj_Vj', float), ('SQM_Vj', float), ('m_Bn', float), ('m_Gn', float),('m_Rn', float),('i_sli', float),('i_msi', float),('i_ipi', float),('Bn_Vj', float),('Gn_Vj', float),('Rn_Vj', float),('m_VIIRS', float),('m_photo', float),('m_scoto', float),('m_SQM', float),('m_Halpha', float),('m_Hbeta', float),('m_Vj', float)])
for X in range(len(listafiles)):
    spec[X].num=X
    spec[X].name=listafiles[X].split('.')[0]
    spec[X].Uj_Vj,spec[X].Bj_Vj,spec[X].Bj_Vj,spec[X].Rn_Vj,spec[X].SQM_Vj,spec[X].m_Bn,spec[X].m_Gn,spec[X].m_Rn,spec[X].i_sli,spec[X].i_msi,spec[X].i_ipi,spec[X].Bn_Vj,spec[X].Gn_Vj,spec[X].Rn_Vj,spec[X].m_VIIRS,spec[X].m_photo,spec[X].m_scoto,spec[X].m_SQM,spec[X].m_Halpha,spec[X].m_Hbeta,spec[X].m_Vj=calculamagnitudes(listafiles[X])

spec = spec[np.array(np.ones(len(spec))-np.isnan(spec['m_Bn'])).astype(bool)]



# Canales DMSP/OLS y VIIRS
# Espectros de lamparas
# Respuesta de Nikon D3s de LICA
# Colores esperados en ISS pictures
#spec= spec[1-np.isnan(spec['Uj_Vj'])]

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
BG=((10**((spec.m_Bn/-2.5)))/10**((spec.m_Gn/-2.5)))
GR=((10**((spec.m_Gn/-2.5)))/10**((spec.m_Rn/-2.5)))
BR=((10**((spec.m_Bn/-2.5)))+10**((spec.m_Rn/-2.5)))/(10**((spec.m_Gn/-2.5)))
BR1=((10**((spec.m_Bn/-2.5)))/10**((spec.m_VIIRS/-2.5)))
G_VIIRS=((10**((spec.m_Gn/-2.5)))/10**((spec.m_VIIRS/-2.5)))
VlG=((10**((spec.m_photo/-2.5)))/10**((spec.m_Gn/-2.5)))
SCG=((10**((spec.m_scoto/-2.5)))/10**((spec.m_Gn/-2.5)))
VlVIIRS=((10**((spec.m_photo/-2.5)))/10**((spec.m_VIIRS/-2.5)))
SQM_G=((10**((spec.m_SQM/-2.5)))/10**((spec.m_Gn/-2.5)))
SQM_Vj=spec.m_SQM-spec.m_Vj
Ha_G=((10**((spec.m_Halpha/-2.5)))/10**((spec.m_Gn/-2.5)))
Hb_G=((10**((spec.m_Hbeta/-2.5)))/10**((spec.m_Gn/-2.5)))

color1=[]
markers=mlines.Line2D.filled_markers
markers1=[]
for X in range(len(spec)):
    if BG[X]<0.05 and GR[X]<0.45:
        color1.append('r')
        markers1.append(markers[0])
    if BG[X]>0.05 and GR[X]<0.45:
        color1.append('k')
        markers1.append(markers[1])        
    if BG[X]<0.34 and GR[X]>0.45 and GR[X]<0.55:
        color1.append('y')
        markers1.append(markers[2])
    if BG[X]<0.34 and GR[X]>0.55:
        color1.append('g')
        markers1.append(markers[3])
    if BG[X]>0.34 and GR[X]>0.55:
        color1.append('b')
        markers1.append(markers[4])


patches = []
fig, ax = plt.subplots()
#ax.plot(BG,GR,c=color1,marker=markers1)
plt.scatter(BG,GR,s=30, c=color1, alpha=0.9)
#ax.scatter(BG,GR, c=c, cmap=cmap)
if 1:
 for X in range(len(listafiles)):
   try:
    plt.text(BG[X],GR[X]+0.01,str(spec[X].num),ha='right', va='bottom')
   except:
    print "fallo "+listafiles[X]
plt.xlabel("B/G")
plt.ylabel("G/R")
plt.title("ISS colors for typical lamps")
plt.xlim(0,0.9)
plt.ylim(0,2.1)
rect = mpatches.Rectangle([0., 0.55], 0.34, 5, ec="none")
patches.append(rect)
rect = mpatches.Rectangle([0.34, 0.55], 5, 5, ec="none")
patches.append(rect)
rect = mpatches.Rectangle([0, 0.45], 0.34,0.1 , ec="none")
patches.append(rect)
rect = mpatches.Rectangle([0.00, 0.], 0.05,0.45 , ec="none")
patches.append(rect)
rect = mpatches.Rectangle([0.05, 0.00], 0.29,0.45 , ec="none")
patches.append(rect)
rect = mpatches.Rectangle([0.34, 0.00], 5,0.55 , ec="none")
patches.append(rect)
colors = np.linspace(0, 1, len(patches))
collection = PatchCollection(patches, cmap=plt.cm.hsv, alpha=0.1)
collection.set_array(np.array(colors))
ax.add_collection(collection)
plt.savefig('colores_espectros2.png')
plt.show()







model_ransac_VlG_BG,line_y_ransac,line_X=ransacfit(BG,VlG)
zVlG_BG = np.polyfit(BG, VlG, 3)
p = np.poly1d(zVlG_BG)
xp = np.linspace(0, 1, 100)
plt.plot(line_X, line_y_ransac, '-b', label='RANSAC regressor')
plt.plot(xp, p(xp),'r-')
plt.scatter(BG,VlG,s=30, c=color1, alpha=0.9)
if 1:
 for X in range(len(listafiles)):
   try:
    plt.text(BG[X],VlG[X]+0.01,str(spec[X].num),ha='right', va='bottom')
   except:
    print "fallo"
plt.xlabel("B/G")
plt.ylabel("Vl/G")
plt.title("ISS colors for typical lamps")
plt.xlim(0,0.9)
plt.ylim(0,3)
plt.savefig('colores_espectros_VlG.png')
plt.show()

model_ransac_SCG_GR,line_y_ransac,line_X=ransacfit(GR,SCG)
zSCGGR = np.polyfit(GR, SCG, 5)
p = np.poly1d(zSCGGR)
xp = np.linspace(0, 2, 100)
plt.plot(line_X, line_y_ransac, '-b', label='RANSAC regressor')
plt.scatter(GR,SCG,s=30, c=color1, alpha=0.9)
plt.plot(xp, p(xp),'r-')
if 1:
 for X in range(len(listafiles)):
   try:
    plt.text(GR[X],SCG[X]+0.01,str(spec[X].num),ha='right', va='bottom')
   except:
    print "fallo"
plt.xlabel("G/R")
plt.ylabel("SC/G")
plt.xlim(0,1.6)
plt.ylim(0,1)
plt.title("ISS colors for typical lamps")
#plt.xlim(0,0.9)
#plt.ylim(0,2.1)
plt.savefig('colores_espectros_SCG_GR.png')
plt.show()

model_ransac_VlG_GR,line_y_ransac,line_X=ransacfit(GR,VlG)
zVlGGR = np.polyfit(GR, VlG, 3)
p = np.poly1d(zVlGGR)
xp = np.linspace(0, 2, 100)
plt.plot(line_X, line_y_ransac, '-b', label='RANSAC regressor')
plt.scatter(GR,VlG,s=30, c=color1, alpha=0.9)
plt.plot(xp, p(xp),'r-')
if 1:
 for X in range(len(listafiles)):
   try:
    plt.text(GR[X],VlG[X]+0.01,str(spec[X].num),ha='right', va='bottom')
   except:
    print "fallo"
plt.xlabel("G/R")
plt.ylabel("Vl/G")
plt.xlim(0,1.6)
plt.ylim(0.5,3)
plt.title("ISS colors for typical lamps")
#plt.xlim(0,0.9)
#plt.ylim(0,2.1)
plt.savefig('colores_espectros_SCG_GR.png')
plt.show()

model_ransac_SQM_G_GR,line_y_ransac,line_X=ransacfit(GR,SQM_G)
zSQM_GGR = np.polyfit(GR, SQM_G, 3)
p = np.poly1d(zSQM_GGR)
xp = np.linspace(0, 2, 100)
plt.plot(line_X, line_y_ransac, '-b', label='RANSAC regressor')
plt.scatter(GR,SQM_G,s=30, c=color1, alpha=0.9)
plt.plot(xp, p(xp),'r-')
if 1:
 for X in range(len(listafiles)):
   try:
    plt.text(GR[X],SQM_G[X]+0.01,str(spec[X].num),ha='right', va='bottom')
   except:
    print "fallo"
plt.xlabel("G/R")
plt.ylabel("SQM/G")
plt.xlim(0,1.6)
plt.ylim(0.5,1.3)
plt.title("ISS colors for typical lamps")
#plt.xlim(0,0.9)
#plt.ylim(0,2.1)
plt.savefig('colores_espectros_SQM_G_GR.png')
plt.show()

model_ransac_SQM_G_BG,line_y_ransac,line_X=ransacfit(BG,SQM_G)
zSQM_GBG = np.polyfit(BG, SQM_G, 3)
p = np.poly1d(zSQM_GBG)
xp = np.linspace(0, 2, 100)
plt.plot(line_X, line_y_ransac, '-b', label='RANSAC regressor')
plt.scatter(BG,SQM_G,s=30, c=color1, alpha=0.9)
plt.plot(xp, p(xp),'r-')
if 1:
 for X in range(len(listafiles)):
   try:
    plt.text(BG[X],SQM_G[X]+0.01,str(spec[X].num),ha='right', va='bottom')
   except:
    print "fallo"
plt.xlabel("B/G")
plt.ylabel("SQM/G")
plt.xlim(0,1.6)
plt.ylim(0.5,1.3)
plt.title("ISS colors for typical lamps")
#plt.xlim(0,0.9)
#plt.ylim(0,2.1)
plt.savefig('colores_espectros_SQM_G_BG.png')
plt.show()

model_ransac_SQM_G_BGj,line_y_ransac,line_X=ransacfit(spec.Bj_Vj,SQM_G)
zSQM_GBGj = np.polyfit(spec.Bj_Vj, SQM_G, 3)
p = np.poly1d(zSQM_GBGj)
xp = np.linspace(0, 2, 100)
plt.plot(line_X, line_y_ransac, '-b', label='RANSAC regressor')
plt.scatter(spec.Bj_Vj,SQM_G,s=30, c=color1, alpha=0.9)
plt.plot(xp, p(xp),'r-')
if 1:
 for X in range(len(listafiles)):
   try:
    plt.text(spec.Bj_Vj[X],SQM_G[X]+0.01,str(spec[X].num),ha='right', va='bottom')
   except:
    print "fallo"
plt.xlabel("Bj - Vj")
plt.ylabel("SQM/G")
plt.xlim(0,5)
plt.ylim(0.5,1.3)
plt.title("ISS colors for typical lamps")
#plt.xlim(0,0.9)
#plt.ylim(0,2.1)
plt.savefig('colores_espectros_SQM_G_BGj.png')
plt.show()

model_ransac_SQM_Vj_BGj,line_y_ransac,line_X=ransacfit(spec.Bj_Vj,SQM_Vj)
zSQM_VjBGj = np.polyfit(spec.Bj_Vj, SQM_Vj, 3)
p = np.poly1d(zSQM_VjBGj)
xp = np.linspace(0, 2, 100)
#plt.plot(line_X, line_y_ransac, '-b', label='RANSAC regressor')
plt.scatter(spec.Bj_Vj,SQM_Vj,s=30, c=color1, alpha=0.9)
#plt.plot(xp, p(xp),'r-')
if 1:
 for X in range(len(listafiles)):
   try:
    plt.text(spec.Bj_Vj[X],SQM_Vj[X]+0.01,str(spec[X].num),ha='right', va='bottom')
   except:
    print "fallo"
plt.xlabel("Bj - Vj")
plt.ylabel("SQM-Vj")
plt.xlim(0,8)
#plt.ylim(0.5,1.3)
plt.title("ISS colors for typical lamps")
#plt.xlim(0,0.9)
#plt.ylim(0,2.1)
plt.savefig('colores_espectros_SQM_Vj_BGj.png')
plt.show()

model_ransac_Ha_GR,line_y_ransac,line_X=ransacfit(GR,Ha_G)
zHa_GGR = np.polyfit(GR, Ha_G, 3)
p = np.poly1d(zHa_GGR)
xp = np.linspace(0, 2, 100)
plt.plot(line_X, line_y_ransac, '-b', label='RANSAC regressor')
plt.scatter(GR,Ha_G,s=30, c=color1, alpha=0.9)
plt.plot(xp, p(xp),'r-')
if 1:
 for X in range(len(listafiles)):
   try:
    plt.text(GR[X],Ha_G[X]+0.01,str(spec[X].num),ha='right', va='bottom')
   except:
    print "fallo"
plt.xlabel("G/R")
plt.ylabel("Ha/G")
plt.xlim(0,1.6)
plt.ylim(0.,2)
plt.title("ISS colors for typical lamps")
#plt.xlim(0,0.9)
#plt.ylim(0,2.1)
plt.savefig('colores_espectros_Ha_G_GR.png')
plt.show()

model_ransac_Hb_GR,line_y_ransac,line_X=ransacfit(GR,Hb_G)
zHb_GGR = np.polyfit(GR, Hb_G, 3)
p = np.poly1d(zHb_GGR)
xp = np.linspace(0, 2, 100)
plt.plot(line_X, line_y_ransac, '-b', label='RANSAC regressor')
plt.scatter(GR,Hb_G,s=30, c=color1, alpha=0.9)
plt.plot(xp, p(xp),'r-')
if 1:
 for X in range(len(listafiles)):
   try:
    plt.text(GR[X],Hb_G[X]+0.01,str(spec[X].num),ha='right', va='bottom')
   except:
    print "fallo"
plt.xlabel("G/R")
plt.ylabel("Hb/G")
plt.xlim(0,1.6)
plt.ylim(0.,1)
plt.title("ISS colors for typical lamps")
#plt.xlim(0,0.9)
#plt.ylim(0,2.1)
plt.savefig('colores_espectros_Hb_G_GR.png')
plt.show()


model_ransac_VIIRSBG,line_y_ransac,line_X=ransacfit(BG,G_VIIRS)
plt.plot(line_X, line_y_ransac, '-b', label='RANSAC regressor')
zgviirsBG = np.polyfit(BG, G_VIIRS, 2)
p = np.poly1d(zgviirsBG)
xp = np.linspace(0, 1, 100)
plt.scatter(BG,G_VIIRS,s=30, c=color1, alpha=0.9)
plt.plot(xp, p(xp),'r-')
if 1:
 for X in range(len(listafiles)):
   try:
    plt.text(BG[X],G_VIIRS[X]+0.01,str(spec[X].num),ha='right', va='bottom')
   except:
    print "fallo"
plt.xlabel("B/G")
plt.ylabel("G/VIIRS")
plt.title("VIIRS correction")
plt.xlim(0,0.9)
plt.ylim(0,2.4)
plt.savefig('colores_espectros_VIIRS_BG.png')
plt.show()

model_ransac_VIIRSGR,line_y_ransac,line_X=ransacfit(GR,G_VIIRS)
plt.plot(line_X, line_y_ransac, '-b', label='RANSAC regressor')
zgviirsGR = np.polyfit(GR, G_VIIRS, 2)
p = np.poly1d(zgviirsGR)
xp = np.linspace(0, 1, 100)
plt.scatter(GR,G_VIIRS,s=30, c=color1, alpha=0.9)
plt.plot(xp, p(xp),'r-')
if 1:
 for X in range(len(listafiles)):
   try:
    plt.text(GR[X],G_VIIRS[X]+0.01,str(spec[X].num),ha='right', va='bottom')
   except:
    print "fallo"
plt.xlabel("G/R")
plt.ylabel("G/VIIRS")
plt.title("VIIRS correction")
plt.xlim(0,1.8)
plt.ylim(0,2.4)
plt.savefig('colores_espectros_VIIRS_GR.png')
plt.show()


model_ransac_VlVIIRSGR,line_y_ransac,line_X=ransacfit(GR,VlVIIRS)
plt.plot(line_X, line_y_ransac, '-b', label='RANSAC regressor')
zgVlVIIRSGR = np.polyfit(GR, VlVIIRS, 2)
p = np.poly1d(zgVlVIIRSGR)
xp = np.linspace(0, 2, 100)
plt.scatter(GR,VlVIIRS,s=30, c=color1, alpha=0.9)
plt.plot(xp, p(xp),'r-')
if 1:
 for X in range(len(listafiles)):
   try:
    plt.text(GR[X],VlVIIRS[X]+0.01,str(spec[X].num),ha='right', va='bottom')
   except:
    print "fallo"
plt.xlabel("G/R")
plt.ylabel("Vl/VIIRS")
plt.title("VIIRS correction")
plt.xlim(0,1.8)
plt.ylim(0,2.4)
plt.savefig('colores_espectros_VlVIIRSGR.png')
plt.show()
        
    
model_ransac_sli,line_y_ransac,line_X=ransacfit(BG,spec.i_sli)
plt.plot(line_X, line_y_ransac, '-b', label='RANSAC regressor')
zslibg = np.polyfit(BG, spec.i_sli, 2)
p = np.poly1d(zslibg)
xp = np.linspace(0, 1, 100)
plt.scatter(BG,spec.i_sli,s=30, c=color1, alpha=0.9)
plt.plot(xp, p(xp),'r-')
if 1:
 for X in range(len(listafiles)):
   try:
    plt.text(BG[X],spec.i_sli[X]+0.01,str(spec[X].num),ha='right', va='bottom')
   except:
    print "fallo"
plt.xlabel("B/G")
plt.ylabel("SLI")
plt.title("SLI")
plt.xlim(0,0.9)
plt.ylim(0,1.2)
plt.savefig('colores_espectros_SLI_BG.png')
plt.show()

model_ransac_sli,line_y_ransac,line_X=ransacfit(GR,spec.i_sli)
plt.plot(line_X, line_y_ransac, '-b', label='RANSAC regressor')
zsligr = np.polyfit(GR, spec.i_sli, 3)
p = np.poly1d(zsligr)
xp = np.linspace(0, 2, 100)
plt.scatter(GR,spec.i_sli,s=30, c=color1, alpha=0.9)
plt.plot(xp, p(xp),'r-')
if 1:
 for X in range(len(listafiles)):
   try:
    plt.text(GR[X],spec.i_sli[X]+0.01,str(spec[X].num),ha='right', va='bottom')
   except:
    print "fallo"
plt.xlabel("G/R")
plt.ylabel("SLI")
plt.title("SLI")
plt.xlim(0,2.1)
plt.ylim(0,1.2)
plt.savefig('colores_espectros_SLI_GR.png')
plt.show()
###################################################3

    
model_ransac_ipi,line_y_ransac,line_X=ransacfit(BG,spec.i_ipi)
plt.plot(line_X, line_y_ransac, '-b', label='RANSAC regressor')
plt.scatter(BG,spec.i_ipi,s=30, c=color1, alpha=0.9)
if 1:
 for X in range(len(listafiles)):
   try:
    plt.text(BG[X],spec.i_ipi[X]+0.01,str(spec[X].num),ha='right', va='bottom')
   except:
    print "fallo"
plt.xlabel("B/G")
plt.ylabel("IPI")
plt.title("IPI")
plt.xlim(0,0.9)
plt.ylim(0,1.2)
plt.savefig('colores_espectros_IPI_BG.png')
plt.show()

plt.scatter(GR,spec.i_ipi,s=30, c=color1, alpha=0.9)
if 1:
 for X in range(len(listafiles)):
   try:
    plt.text(GR[X],spec.i_ipi[X]+0.01,str(spec[X].num),ha='right', va='bottom')
   except:
    print "fallo"
plt.xlabel("G/R")
plt.ylabel("IPI")
plt.title("IPI")
plt.xlim(0,2.1)
plt.ylim(0,1.2)
plt.savefig('colores_espectros_IPI_GR.png')
plt.show()

plt.scatter(BR,spec.i_ipi,s=30, c=color1, alpha=0.9)
if 1:
 for X in range(len(listafiles)):
   try:
    plt.text(BR[X],spec.i_ipi[X]+0.01,str(spec[X].num),ha='right', va='bottom')
   except:
    print "fallo"
plt.xlabel("(B+R)/G")
plt.ylabel("IPI")
plt.title("IPI")
#plt.xlim(0,1.2)
#plt.ylim(0,1.2)
plt.savefig('colores_espectros_IPI_BRG0.png')
plt.show()

plt.scatter(BR1,spec.i_ipi,s=30, c=color1, alpha=0.9)
if 0:
 for X in range(len(listafiles)):
   try:
    plt.text(BR1[X],spec.i_ipi[X]+0.01,str(spec[X].num),ha='right', va='bottom')
   except:
    print "fallo"
plt.xlabel("(B-G)/R")
plt.ylabel("IPI")
plt.title("IPI")
#plt.xlim(0,1.2)
#plt.ylim(0,1.2)
plt.savefig('colores_espectros_IPI_BRG1.png')
plt.show()   
###################################################3

    
model_ransac_msi,line_y_ransac,line_X=ransacfit(BG,spec.i_msi)
plt.plot(line_X, line_y_ransac, '-b', label='RANSAC regressor')
plt.scatter(BG,spec.i_msi,s=30, c=color1, alpha=0.9)
if 1:
 for X in range(len(listafiles)):
   try:
    plt.text(BG[X],spec.i_msi[X]+0.01,str(spec[X].num),ha='right', va='bottom')
   except:
    print "fallo"
plt.xlabel("B/G")
plt.ylabel("MSI")
plt.title("MSI")
plt.xlim(0,0.9)
plt.ylim(0,1.2)
plt.savefig('colores_espectros_MSI_BG.png')
plt.show()

model_ransac_msigr,line_y_ransac,line_X=ransacfit(GR,spec.i_msi)
plt.plot(line_X, line_y_ransac, '-b', label='RANSAC regressor')
zmsigr = np.polyfit(GR, spec.i_msi, 3)
p = np.poly1d(zsligr)
xp = np.linspace(0, 2, 100)
plt.scatter(GR,spec.i_msi,s=30, c=color1, alpha=0.9)
plt.plot(xp, p(xp),'r-')
if 1:
 for X in range(len(listafiles)):
   try:
    plt.text(GR[X],spec.i_msi[X]+0.01,str(spec[X].num),ha='right', va='bottom')
   except:
    print "fallo"
plt.xlabel("G/R")
plt.ylabel("MSI")
plt.title("MSI")
plt.xlim(0,2.1)
plt.ylim(0,1.2)
plt.savefig('colores_espectros_MSI_GR.png')
plt.show()   

asciitable.write(spec,'espectros.txt')


from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
#X, Y, Z = axes3d.get_test_data(0.05)
ax.scatter(BG, GR, spec.i_msi, c='r', marker='o')

ax.set_xlabel('B/G')
ax.set_ylabel('G/R')
ax.set_zlabel('MSI')
ax.set_xlim(0,1)
ax.set_ylim(0,2.1)
ax.set_zlim(0,1.2)

plt.show()

XYZ=np.vstack((BG,GR,spec.i_msi))


XYZ1=XYZ.transpose()
XYZ2=[]


plane=fitPlaneOptimize(XYZ1)

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

point  = np.mean(XYZ1,0)
normal = fitPLaneLTSQ(XYZ1)

# a plane is a*x+b*y+c*z+d=0
# [a,b,c] is the normal. Thus, we have to calculate
# d and we're set
d = -point.dot(normal)

# create x,y
xx, yy = np.meshgrid(range(3), range(3))

# calculate corresponding z
z = (-normal[0] * xx - normal[1] * yy - d) * 1. /normal[2]

# plot the surface
plt3d = plt.figure().gca(projection='3d')

#plt3d.plot_surface(xx, yy, z)
plt3d.plot_wireframe(xx, yy, z, rstride=1, cstride=1)
plt3d.scatter(BG, GR, spec.i_msi, c='r', marker='o')
plt3d.set_title('MSI fit by B/G and G/R')
plt3d.set_xlim(0,1)
plt3d.set_ylim(0,2.1)
plt3d.set_zlim(0,1.2)
plt3d.set_xlabel('B/G')
plt3d.set_ylabel('G/R')
plt3d.set_zlabel('MSI')
plt.show()


##################################
#
# Fiting plane IPI
#
####################################

XYZ=np.vstack((BG,GR,spec.i_ipi))


XYZ1=XYZ.transpose()
XYZ2=[]


plane=fitPlaneOptimize(XYZ1)

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

point  = np.mean(XYZ1,0)
normal = fitPLaneLTSQ(XYZ1)

# a plane is a*x+b*y+c*z+d=0
# [a,b,c] is the normal. Thus, we have to calculate
# d and we're set
d = -point.dot(normal)

# create x,y
xx, yy = np.meshgrid(range(3), range(3))

# calculate corresponding z
z = (-normal[0] * xx - normal[1] * yy - d) * 1. /normal[2]

# plot the surface
plt3d = plt.figure().gca(projection='3d')

#plt3d.plot_surface(xx, yy, z)
plt3d.plot_wireframe(xx, yy, z, rstride=1, cstride=1)
plt3d.scatter(BG, GR, spec.i_ipi, c='r', marker='o')
plt3d.set_title('IPI fit by B/G and G/R')
plt3d.set_xlim(0,1)
plt3d.set_ylim(0,2.1)
plt3d.set_zlim(0,1.2)
plt3d.set_xlabel('B/G')
plt3d.set_ylabel('G/R')
plt3d.set_zlabel('IPI')
plt.show()

plt.plot((-normal[0] * BG - normal[1] * GR - d) * 1. /normal[2],spec.i_ipi,'ko')
if 0:
 for X in range(len(listafiles)):
   try:
    plt.text((-normal[0] * BG - normal[1] * GR - d) * 1. /normal[2],spec.i_ipi[X]+0.01,str(spec[X].num),ha='right', va='bottom')
   except:
    print "fallo"
plt.xlabel("Combination(BGR)")
plt.ylabel("IPI")
plt.title("IPI")
#plt.xlim(0,1.2)
#plt.ylim(0,1.2)
plt.savefig('colores_espectros_IPI_BRG1.png')
plt.show() 


##################################
#
# Fiting plane SLI
#
####################################

XYZ=np.vstack((BG,GR,spec.i_sli))


XYZ1=XYZ.transpose()
XYZ2=[]


plane=fitPlaneOptimize(XYZ1)

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

point  = np.mean(XYZ1,0)
normal = fitPLaneLTSQ(XYZ1)

# a plane is a*x+b*y+c*z+d=0
# [a,b,c] is the normal. Thus, we have to calculate
# d and we're set
d = -point.dot(normal)

# create x,y
xx, yy = np.meshgrid(range(3), range(3))

# calculate corresponding z
z = (-normal[0] * xx - normal[1] * yy - d) * 1. /normal[2]

# plot the surface
plt3d = plt.figure().gca(projection='3d')

#plt3d.plot_surface(xx, yy, z)
plt3d.plot_wireframe(xx, yy, z, rstride=1, cstride=1)
plt3d.scatter(BG, GR, spec.i_sli, c='r', marker='o')
plt3d.set_title('SLI fit by B/G and G/R')
plt3d.set_xlim(0,1)
plt3d.set_ylim(0,2.1)
plt3d.set_zlim(0,1.2)
plt3d.set_xlabel('B/G')
plt3d.set_ylabel('G/R')
plt3d.set_zlabel('SLI')
plt.show()

XYZ=np.vstack((BG,GR,G_VIIRS))


XYZ1=XYZ.transpose()
XYZ2=[]


plane=fitPlaneOptimize(XYZ1)

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

point  = np.mean(XYZ1,0)
normal = fitPLaneLTSQ(XYZ1)

# a plane is a*x+b*y+c*z+d=0
# [a,b,c] is the normal. Thus, we have to calculate
# d and we're set
d = -point.dot(normal)

# create x,y
xx, yy = np.meshgrid(range(6), range(6))

# calculate corresponding z
z = (-normal[0] * xx - normal[1] * yy - d) * 1. /normal[2]

# plot the surface
plt3d = plt.figure().gca(projection='3d')

#plt3d.plot_surface(xx, yy, z)
plt3d.plot_wireframe(xx, yy, z, rstride=1, cstride=1)
plt3d.scatter(BG, GR,G_VIIRS , c='r', marker='o')
plt3d.set_title('VIIRS fit by B/G and G/R')
plt3d.set_xlim(0,1)
plt3d.set_ylim(0,2.1)
plt3d.set_zlim(0,2.6)
plt3d.set_xlabel('B/G')
plt3d.set_ylabel('G/R')
plt3d.set_zlabel('G/VIIRS')
plt.show()

print "RANSAC VlG"+str(model_ransac_sli.estimator_.predict(0))+' ' +str(model_ransac_sli.estimator_.coef_)
print "RANSAC MSI"+str(model_ransac_msi.estimator_.predict(0))+' ' +str(model_ransac_msi.estimator_.coef_)
print "RANSAC MSI GR"+str(model_ransac_msigr.estimator_.predict(0))+' ' +str(model_ransac_msigr.estimator_.coef_)
print "RANSAC IPI"+str(model_ransac_ipi.estimator_.predict(0))+' ' +str(model_ransac_ipi.estimator_.coef_)
print "RANSAC SLI"+str(model_ransac_sli.estimator_.predict(0))+' ' +str(model_ransac_sli.estimator_.coef_)
print "RANSAC VIIRS"+str(model_ransac_VIIRSGR.estimator_.predict(0))+' ' +str(model_ransac_VIIRSGR.estimator_.coef_)
print "SLI-BG"
print zslibg
print "SLI-GR"
print zsligr
print "MSI-GR"
print zmsigr
print "VlG-GR"
print zVlGGR
