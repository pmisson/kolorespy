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

#import matplotlib.style
#matplotlib.style.use('ggplot')
import matplotlib.pyplot as plt
import numpy as np
import random
import pandas as pd
from matplotlib import rcParams
rcParams.update({'figure.autolayout':True})
alpha=0.9
if 0:
    d = {'x': x, 'y': y, 'color': c}

    df = pd.DataFrame(d)

    mm = df.groupby('color')

    cat = ['R', 'G', 'B']
    char = ['*', '.', 'o']
def scatter2(cat,char,mm,size,alpha):
    for ct, ch in zip(cat, char):
       tab = mm.get_group(ct)
       plt.scatter(tab['x'], tab['y'], marker=ch,color=ct,s=size,alpha=alpha)

def scatter3(cat,char2,spec1,size,alpha):
    n=range(len(spec1))
    dict={'HPS':'*','INC':'o','HAL':'>','CMH':'<','CFL':'v','FL':'v','MV':'s','LPS':'d','LED':'p','MH':'1'}
    for ct, ch,X in zip(cat, char2,n):
       
       plt.scatter(ct,ch,marker=dict[spec1[X].tipo],color=spec1[X].color,s=size,alpha=alpha)
    for X in range(len(dict)):    
        plt.plot(0,0,dict[list(dict)[X]],c='b',label=list(dict)[X],ms=10)
        plt.plot(0,0,'.',c='w',ms=0)
    plt.legend(numpoints=1, borderaxespad=0.,fontsize=14,mode="expand",ncol=5,bbox_to_anchor=(0., 0., 1., .102))
    
    #plt.show()



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
alpha
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
NO2q=asciitable.read("NO2q.csv")
NOq=asciitable.read("NOq.csv")
NO3s=asciitable.read("NO3s.csv")

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

def jNO3(espectro,NO3s,NO2q,NOq,photo): #Calculate Index as is definide in http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0067798
    s = InterpolatedUnivariateSpline(espectro.col2, espectro.col1, k=3)
    NO3 =  InterpolatedUnivariateSpline(NO3s.wave*10, NO3s.flux, k=3)
    NO2 = InterpolatedUnivariateSpline(NO2q.wave*10, NO2q.flux,k=3)
    NO = InterpolatedUnivariateSpline(NOq.wave*10, NOq.flux,k=3)
    photo2 = InterpolatedUnivariateSpline(photo.wave, photo.intrel,k=3)
    xs = np.arange(4100, 6400, 10)
    s_1 = s(xs)
    ph_1=photo2(xs)
    S1n=s_1*ph_1
    NormS=np.sum(S1n)
    s_1x=s_1/NormS
    NO3_1 = NO3(xs)
    NO2_1= NO2(xs)
    NO_1=NO(xs)
    NOT=NO3_1*(NO2_1+NO_1)
    ANOT=np.sum(NOT)
    All_1=s_1x*NOT
    result=np.sum(All_1/ANOT)

    
    return result




def calculamagnitudes(filename,clouds):
    spectro=asciitable.read(filename)
    tipo=filename.split('_')[0][:3].upper()
    char2= ['*','o','>','<','V','V','s','d','D','p']
    tp =   ['HPS','INC','HAL','CMH','CFL','FL','MV','LPS','OTH','LED']
    specamply=asciitable.read(clouds)
    mask=spectro.col1<0
    spectro.col1[mask]=0
    cen=spectro.col2<7000
    cen2=spectro.col2<4000
    cen3=cen-cen2
    wave=np.sum(spectro.col2[cen3]*spectro.col1[cen3])/np.sum(spectro.col1[cen3])
    sA = InterpolatedUnivariateSpline(specamply.col1, specamply.col2, k=3)
    s_A=sA(spectro.col2)
    spectroA=spectro.copy()
    spectroA.col1=spectroA.col1*(s_A)
    m_Uj,I_Uj=mag(U_johnson,spectro,AB)
    m_Bj,I_Bj=mag(B_johnson,spectro,AB)
    m_BjA,I_BjA=mag(B_johnson,spectroA,AB)
    m_Vj,I_Vj=mag(V_johnson,spectro,AB)
    m_VjA,I_VjA=mag(V_johnson,spectroA,AB)
    m_Rj,I_Vj=mag(R_johnson,spectro,AB)
    m_Bn,I_Bn=mag(B_D3,spectro,AB)
    m_Gn,I_Gn=mag(G_D3,spectro,AB)
    m_Rn,I_Rn=mag(R_D3,spectro,AB)
    m_SQM,I_Rn=mag(SQM_LICA,spectro,AB)
    m_SQMA,I_RnA=mag(SQM_LICA,spectroA,AB)
    m_VIIRS,I_VIIRS=mag(VIIRS_2013,spectro,AB)
    m_DMSP,I_DMSP=mag(DMSP_1999,spectro,AB)
    m_photo,I_photo=mag(photo,spectro,AB)
    m_scoto,I_scoto=mag(scoto,spectro,AB)
    m_photoA,I_photoA=mag(photo,spectroA,AB)
    m_scotoA,I_scotoA=mag(scoto,spectroA,AB)
    m_Halpha,I_Halpha=mag(Halpha,spectro,AB)
    m_Hbeta,I_Hbeta=mag(Hbeta,spectro,AB)
    Uj_Vj=m_Uj-m_Vj
    Bj_Vj=m_Bj-m_Vj
    Bj_VjA=m_BjA-m_VjA
    Rj_Vj=m_Rj-m_Vj
    Bn_Vj=m_Bn-m_Vj
    Gn_Vj=m_Gn-m_Vj
    Rn_Vj=m_Rn-m_Vj
    SQM_Vj=m_SQM-m_Vj
    SQM_VjA=m_SQMA-m_VjA
    Bn_Gn=I_Bn/I_Gn
    Gn_Rn=I_Gn/I_Rn
    i_sli=IndiceIsoLux(scoto,spectro,sun_lambda,photo)
    i_msi=IndiceIsoLux(msas,spectro,sun_lambda,photo)
    i_ipi=IndiceIsoLux(psas,spectro,sun_lambda,photo)
    jNO=jNO3(spectro,NO3s,NO2q,NOq,photo)
    if 0:
        fig, (ax0,ax1,ax2) = plt.subplots(nrows=3)

        ax0.plot(WL_VEGA_AB, FLUX_VEGA_AB)
        ax0.plot(WL_VEGA_ST, FLUX_VEGA_ST)
        ax0.set_title(filename)
        ax1.plot(spectro.col2, spectro.col1)
        #ax1.set_title('Spectra')
        ax2.plot(R_D3s.wave,R_D3s.intrel)
        ax2.plot(G_D3s.wave,G_D3s.intrel)
        ax2.plot(B_D3s.wave,B_D3s.intrel)
        #ax2.plot(photo.wave,photo.intrel)
        ax2.plot(psas.wave,psas.intrel)
        #ax2.plot(scoto.wave,scoto.intrel)
        #ax2.plot(msas.wave,msas.intrel)
        ax2.set_title('Filter')
        plt.show()
    return Uj_Vj,Bj_Vj,Bj_VjA,Rj_Vj,SQM_Vj,SQM_VjA,m_Bn,m_Gn,m_Rn,i_sli,i_msi,i_ipi,Bn_Vj,Gn_Vj,Rn_Vj,m_VIIRS,m_photo,m_scoto,m_photoA,m_scotoA,m_SQM,m_SQMA,m_Halpha,m_Hbeta,m_Vj,m_VjA,jNO,wave,tipo

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

path=os.getcwd()
lista=os.listdir(path)
lista5=" ".join(lista)
lista4=re.findall(r'[a-zA-Z_0-9+.]+\.ds',lista5)
lista4.sort()
listafiles2=lista4

dict={'HPS':'*','INC':'o','HAL':'>','CMH':'<','CFL':'v','FL':'v','MV':'s','LPS':'d','LED':'p','MH':'1'}
clouds='sum_F0.05_AOD_0.1_Hc0.5km.ssc'

spec=np.recarray((len(listafiles),),dtype=[('num', int),('name', 'S100'),('Uj_Vj', float), ('Bj_Vj', float),('Bj_VjA', float),('Rj_Vj', float), ('SQM_Vj', float),('SQM_VjA', float), ('m_Bn', float), ('m_Gn', float),('m_Rn', float),('i_sli', float),('i_msi', float),('i_ipi', float),('Bn_Vj', float),('Gn_Vj', float),('Rn_Vj', float),('m_VIIRS', float),('m_photo', float),('m_scoto', float),('m_photoA', float),('m_scotoA', float),('m_SQM', float),('m_SQMA', float),('m_Halpha', float),('m_Hbeta', float),('m_Vj', float),('m_VjA', float),('jNO',float),('wave', float),('tipo', 'S3')])
for X in range(len(listafiles)):
    spec[X].num=X
    spec[X].name=listafiles[X].split('.')[0]
    spec[X].Uj_Vj,spec[X].Bj_Vj,spec[X].Bj_VjA,spec[X].Rn_Vj,spec[X].SQM_Vj,spec[X].SQM_VjA,spec[X].m_Bn,spec[X].m_Gn,spec[X].m_Rn,spec[X].i_sli,spec[X].i_msi,spec[X].i_ipi,spec[X].Bn_Vj,spec[X].Gn_Vj,spec[X].Rn_Vj,spec[X].m_VIIRS,spec[X].m_photo,spec[X].m_scoto,spec[X].m_photoA,spec[X].m_scotoA,spec[X].m_SQM,spec[X].m_SQMA,spec[X].m_Halpha,spec[X].m_Hbeta,spec[X].m_Vj,spec[X].m_VjA,spec[X].jNO,spec[X].wave,spec[X].tipo=calculamagnitudes(listafiles[X],clouds)

spec = spec[np.array(np.ones(len(spec))-np.isnan(spec['m_Bn'])).astype(bool)]

spec2=np.recarray((len(listafiles2),),dtype=[('num', int),('name', 'S100'),('Uj_Vj', float), ('Bj_Vj', float),('Bj_VjA', float),('Rj_Vj', float), ('SQM_Vj', float),('SQM_VjA', float), ('m_Bn', float), ('m_Gn', float),('m_Rn', float),('i_sli', float),('i_msi', float),('i_ipi', float),('Bn_Vj', float),('Gn_Vj', float),('Rn_Vj', float),('m_VIIRS', float),('m_photo', float),('m_scoto', float),('m_photoA', float),('m_scotoA', float),('m_SQM', float),('m_SQMA', float),('m_Halpha', float),('m_Hbeta', float),('m_Vj', float),('m_VjA', float),('jNO',float),('wave', float),('tipo', 'S3')])
for X in range(len(listafiles2)):
    spec2[X].num=X
    spec2[X].name=listafiles2[X].split('.')[0]
    spec2[X].Uj_Vj,spec2[X].Bj_Vj,spec2[X].Bj_VjA,spec2[X].Rn_Vj,spec2[X].SQM_Vj,spec2[X].SQM_VjA,spec2[X].m_Bn,spec2[X].m_Gn,spec2[X].m_Rn,spec2[X].i_sli,spec2[X].i_msi,spec2[X].i_ipi,spec2[X].Bn_Vj,spec2[X].Gn_Vj,spec2[X].Rn_Vj,spec2[X].m_VIIRS,spec2[X].m_photo,spec2[X].m_scoto,spec2[X].m_photoA,spec2[X].m_scotoA,spec2[X].m_SQM,spec2[X].m_SQMA,spec2[X].m_Halpha,spec2[X].m_Hbeta,spec2[X].m_Vj,spec2[X].m_VjA,spec2[X].jNO,spec2[X].wave,spec2[X].tipo=calculamagnitudes(listafiles2[X],clouds)

spec2 = spec2[np.array(np.ones(len(spec2))-np.isnan(spec2['m_Bn'])).astype(bool)]



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
SQMVl=((10**((spec.m_SQM/-2.5)))/10**((spec.m_photo/-2.5)))
SQMSC=((10**((spec.m_SQM/-2.5)))/10**((spec.m_scoto/-2.5)))
SQM_Vl=spec.m_SQM-spec.m_photo
SQM_SC=spec.m_SQM-spec.m_scoto
SQM_VjA=spec.m_SQMA-spec.m_VjA
SQM_VlA=spec.m_SQMA-spec.m_photoA
SQM_SCA=spec.m_SQMA-spec.m_scotoA
SQM_SQMA=spec.m_SQMA-spec.m_SQM
SQMASQM=((10**((spec.m_SQMA/-2.5)))/10**((spec.m_SQM/-2.5)))

BG2=((10**((spec2.m_Bn/-2.5)))/10**((spec2.m_Gn/-2.5)))
GR2=((10**((spec2.m_Gn/-2.5)))/10**((spec2.m_Rn/-2.5)))
BR2=((10**((spec2.m_Bn/-2.5)))+10**((spec2.m_Rn/-2.5)))/(10**((spec2.m_Gn/-2.5)))
BR12=((10**((spec2.m_Bn/-2.5)))/10**((spec2.m_VIIRS/-2.5)))
G_VIIRS2=((10**((spec2.m_Gn/-2.5)))/10**((spec2.m_VIIRS/-2.5)))
VlG2=((10**((spec2.m_photo/-2.5)))/10**((spec2.m_Gn/-2.5)))
SCG2=((10**((spec2.m_scoto/-2.5)))/10**((spec2.m_Gn/-2.5)))
VlVIIRS2=((10**((spec2.m_photo/-2.5)))/10**((spec2.m_VIIRS/-2.5)))
SQM_G2=((10**((spec2.m_SQM/-2.5)))/10**((spec2.m_Gn/-2.5)))
SQM_Vj2=spec2.m_SQM-spec2.m_Vj
SQMVl2=((10**((spec2.m_SQM/-2.5)))/10**((spec2.m_photo/-2.5)))
SQMSC2=((10**((spec2.m_SQM/-2.5)))/10**((spec2.m_scoto/-2.5)))

SQM_Vl2=spec2.m_SQM-spec2.m_photo
SQM_SC2=spec2.m_SQM-spec2.m_scoto
SQM_VjA2=spec2.m_SQMA-spec2.m_VjA
SQM_VlA2=spec2.m_SQMA-spec2.m_photoA
SQM_SCA2=spec2.m_SQMA-spec2.m_scotoA
SQM_SQMA2=spec2.m_SQMA-spec2.m_SQM


color1=[]
markers=mlines.Line2D.filled_markers
markers=['o','v','^','<','>','1','*','h']
markers1=[]

SQM_SC_r=[]
SQM_Vl_r=[]
SQM_Vj_r=[]
SQM_SQMA_r=[]
SQMASQM_r=[]
B_V_r=[]
wave_r=[]
SQM_SCA_r=[]
SQM_VlA_r=[]
SQM_VjA_r=[]
SQM_SQMA_r=[]
B_VA_r=[]
SQM_SC_y=[]
SQM_Vl_y=[]
SQM_Vj_y=[]
SQM_SQMA_y=[]
SQMASQM_y=[]
B_V_y=[]
wave_y=[]
SQM_SCA_y=[]
SQM_VlA_y=[]
SQM_VjA_y=[]
SQM_SQMA_y=[]
B_VA_y=[]
SQM_SC_g=[]
SQM_Vl_g=[]
SQM_Vj_g=[]
SQM_SQMA_g=[]
SQMASQM_g=[]
B_V_g=[]
wave_g=[]
SQM_SCA_g=[]
SQM_VlA_g=[]
SQM_VjA_g=[]
SQM_SQMA_g=[]
B_VA_g=[]
SQM_SC_b=[]
SQM_Vl_b=[]
SQM_Vj_b=[]
SQM_SQMA_b=[]
SQMASQM_b=[]
B_V_b=[]
wave_b=[]
SQM_SCA_b=[]
SQM_VlA_b=[]
SQM_VjA_b=[]
SQM_SQMA_b=[]
B_VA_b=[]
SQM_SC_m=[]
SQM_Vl_m=[]
SQM_Vj_m=[]
SQM_SQMA_m=[]
SQMASQM_m=[]
B_V_m=[]
wave_m=[]
SQM_SCA_m=[]
SQM_VlA_m=[]
SQM_VjA_m=[]
SQM_SQMA_m=[]
B_VA_m=[]
SQM_SC_w=[]
SQM_Vl_w=[]
SQM_Vj_w=[]
SQM_SQMA_w=[]
SQMASQM_w=[]
B_V_w=[]
wave_w=[]
SQM_SCA_w=[]
SQM_VlA_w=[]
SQM_VjA_w=[]
SQM_SQMA_w=[]
B_VA_w=[]

for X in range(len(spec)):
    if BG[X]<0.05 and GR[X]<0.45:
        color1.append('r')
        markers1.append(markers[0])
        SQM_SC_r.append(SQM_SC[X])
        SQM_Vl_r.append(SQM_Vl[X])
        SQM_Vj_r.append(SQM_Vj[X])
        B_V_r.append(spec.Bj_Vj[X])
        SQM_SCA_r.append(SQM_SCA[X])
        SQM_VlA_r.append(SQM_VlA[X])
        SQM_VjA_r.append(SQM_VjA[X])
        SQM_SQMA_r.append(SQM_SQMA[X])
        B_VA_r.append(spec.Bj_VjA[X])
        wave_r.append(spec.wave[X])
        SQMASQM_r.append(SQMASQM[X])
    elif BG[X]>0.05 and GR[X]<0.36 and BG[X]<0.25:
        color1.append('k')
        markers1.append(markers[1])
        SQM_SC_y.append(SQM_SC[X])
        SQM_Vl_y.append(SQM_Vl[X])
        SQM_Vj_y.append(SQM_Vj[X])
        B_V_y.append(spec.Bj_Vj[X])
        SQM_SCA_y.append(SQM_SCA[X])
        SQM_VlA_y.append(SQM_VlA[X])
        SQM_VjA_y.append(SQM_VjA[X])
        SQM_SQMA_y.append(SQM_SQMA[X])
        B_VA_y.append(spec.Bj_VjA[X])
        wave_y.append(spec.wave[X])
        SQMASQM_y.append(SQMASQM[X])
    elif BG[X]>0.25 and GR[X]<0.56:
        color1.append('c')
        markers1.append(markers[6])
        SQM_SC_w.append(SQM_SC[X])
        SQM_Vl_w.append(SQM_Vl[X])
        SQM_Vj_w.append(SQM_Vj[X])
        B_V_w.append(spec.Bj_Vj[X])
        SQM_SCA_w.append(SQM_SCA[X])
        SQM_VlA_w.append(SQM_VlA[X])
        SQM_VjA_w.append(SQM_VjA[X])
        SQM_SQMA_w.append(SQM_SQMA[X])
        B_VA_w.append(spec.Bj_VjA[X])
        wave_w.append(spec.wave[X])
        SQMASQM_w.append(SQMASQM[X])                  
    elif BG[X]<0.25 and GR[X]>0.36 and GR[X]<0.55:
        color1.append('y')
        markers1.append(markers[2])
        SQM_SC_y.append(SQM_SC[X])
        SQM_Vl_y.append(SQM_Vl[X])
        SQM_Vj_y.append(SQM_Vj[X])
        B_V_y.append(spec.Bj_Vj[X])
        SQM_SCA_y.append(SQM_SCA[X])
        SQM_VlA_y.append(SQM_VlA[X])
        SQM_VjA_y.append(SQM_VjA[X])
        SQM_SQMA_y.append(SQM_SQMA[X])
        B_VA_y.append(spec.Bj_VjA[X])    
        wave_y.append(spec.wave[X])
        SQMASQM_y.append(SQMASQM[X])
    elif BG[X]>0.46 and GR[X]>0 and GR[X]<0.55:
        color1.append('w')
        markers1.append(markers[6])
        SQM_SC_w.append(SQM_SC[X])
        SQM_Vl_w.append(SQM_Vl[X])
        SQM_Vj_w.append(SQM_Vj[X])
        B_V_w.append(spec.Bj_Vj[X])
        SQM_SQMA_w.append(SQM_SQMA[X])
        SQM_SCA_w.append(SQM_SCA[X])
        SQM_VlA_w.append(SQM_VlA[X])
        SQM_VjA_w.append(SQM_VjA[X])
        B_VA_w.append(spec.Bj_VjA[X])
        wave_w.append(spec.wave[X])
        SQMASQM_w.append(SQMASQM[X])  
    elif BG[X]<0.34 and GR[X]>0.55:
        color1.append('g')
        markers1.append(markers[3])
        SQM_SC_g.append(SQM_SC[X])
        SQM_Vl_g.append(SQM_Vl[X])
        SQM_Vj_g.append(SQM_Vj[X])
        B_V_g.append(spec.Bj_Vj[X])
        SQM_SCA_g.append(SQM_SCA[X])
        SQM_VlA_g.append(SQM_VlA[X])
        SQM_VjA_g.append(SQM_VjA[X])
        B_VA_g.append(spec.Bj_VjA[X])
        SQM_SQMA_g.append(SQM_SQMA[X])
        wave_g.append(spec.wave[X])
        SQMASQM_g.append(SQMASQM[X])    
    elif BG[X]>0.34 and GR[X]>0.55:
        color1.append('b')
        markers1.append(markers[4])
        SQM_SC_b.append(SQM_SC[X])
        SQM_Vl_b.append(SQM_Vl[X])
        SQM_Vj_b.append(SQM_Vj[X])
        B_V_b.append(spec.Bj_Vj[X])
        SQM_SCA_b.append(SQM_SCA[X])
        SQM_VlA_b.append(SQM_VlA[X])
        SQM_VjA_b.append(SQM_VjA[X])
        B_VA_b.append(spec.Bj_VjA[X])
        SQM_SQMA_b.append(SQM_SQMA[X])
        wave_b.append(spec.wave[X])
        SQMASQM_b.append(SQMASQM[X]) 
    elif BG[X]<0.34 and GR[X]>0.2 and GR[X]<0.46:
        color1.append('m')
        markers1.append(markers[5])
        SQM_SC_m.append(SQM_SC[X])
        SQM_Vl_m.append(SQM_Vl[X])
        SQM_Vj_m.append(SQM_Vj[X])
        B_V_m.append(spec.Bj_Vj[X])
        SQM_SCA_m.append(SQM_SCA[X])
        SQM_VlA_m.append(SQM_VlA[X])
        SQM_VjA_m.append(SQM_VjA[X])
        SQM_SQMA_m.append(SQM_SQMA[X])
        B_VA_m.append(spec.Bj_VjA[X])
        wave_m.append(spec.wave[X])
        SQMASQM_m.append(SQMASQM[X])
    else:
       color1.append('m') 
       markers1.append(markers[5])
       SQM_SC_m.append(SQM_SC[X])
       SQM_Vl_m.append(SQM_Vl[X])
       SQM_Vj_m.append(SQM_Vj[X])
       B_V_m.append(spec.Bj_Vj[X])
       SQM_SCA_m.append(SQM_SCA[X])
       SQM_VlA_m.append(SQM_VlA[X])
       SQM_VjA_m.append(SQM_VjA[X])
       SQM_SQMA_m.append(SQM_SQMA[X])
       B_VA_m.append(spec.Bj_VjA[X])
       wave_m.append(spec.wave[X])
       SQMASQM_m.append(SQMASQM[X])

spec1=mlab.rec_append_fields(spec,'color',np.array(color1))
#spec.color=np.array(color1)
SQMASQM_t=np.array([np.mean(SQMASQM_r),np.mean(SQMASQM_y),np.mean(SQMASQM_g),np.mean(SQMASQM_b),np.mean(SQMASQM_w)])#,np.mean(SQMASQM_m)])
SQM_SC_t=np.array([np.mean(SQM_SC_r),np.mean(SQM_SC_y),np.mean(SQM_SC_g),np.mean(SQM_SC_b),np.mean(SQM_SC_w)])#,np.mean(SQM_SC_m)])
SQM_Vl_t=np.array([np.mean(SQM_Vl_r),np.mean(SQM_Vl_y),np.mean(SQM_Vl_g),np.mean(SQM_Vl_b),np.mean(SQM_Vl_w)])#,np.mean(SQM_Vl_m)])
SQM_Vj_t=np.array([np.mean(SQM_Vj_r),np.mean(SQM_Vj_y),np.mean(SQM_Vj_g),np.mean(SQM_Vj_b),np.mean(SQM_Vj_w)])#,np.mean(SQM_Vj_m)])
B_V_t=np.array([np.mean(B_V_r),np.mean(B_V_y),np.mean(B_V_g),np.mean(B_V_b),np.mean(B_V_w)])#,np.mean(B_V_m)])
wave_t=np.array([np.mean(wave_r),np.mean(wave_y),np.mean(wave_g),np.mean(wave_b),np.mean(wave_w)])#,np.mean(wave_m)])
SQM_SCA_t=np.array([np.mean(SQM_SCA_r),np.mean(SQM_SCA_y),np.mean(SQM_SCA_g),np.mean(SQM_SCA_b),np.mean(SQM_SCA_w)])#,np.mean(SQM_SCA_m)])
SQM_VlA_t=np.array([np.mean(SQM_VlA_r),np.mean(SQM_VlA_y),np.mean(SQM_VlA_g),np.mean(SQM_VlA_b),np.mean(SQM_VlA_w)])#,np.mean(SQM_VlA_m)])
SQM_VjA_t=np.array([np.mean(SQM_VjA_r),np.mean(SQM_VjA_y),np.mean(SQM_VjA_g),np.mean(SQM_VjA_b),np.mean(SQM_VjA_w)])#,np.mean(SQM_VjA_m)])
SQM_SQMA_t=np.array([np.mean(SQM_SQMA_r),np.mean(SQM_SQMA_y),np.mean(SQM_SQMA_g),np.mean(SQM_SQMA_b),np.mean(SQM_SQMA_w)])#,np.mean(SQM_SQMA_m)])
B_VA_t=np.array([np.mean(B_VA_r),np.mean(B_VA_y),np.mean(B_VA_g),np.mean(B_VA_b),np.mean(B_VA_w)])#,np.mean(B_VA_m)])
SQM_SC_tmax=np.array([max(SQM_SC_r),max(SQM_SC_y),max(SQM_SC_g),max(SQM_SC_b),max(SQM_SC_w)])#,max(SQM_SC_m)])
SQM_Vl_tmax=np.array([max(SQM_Vl_r),max(SQM_Vl_y),max(SQM_Vl_g),max(SQM_Vl_b),max(SQM_Vl_w)])#,max(SQM_Vl_m)])
SQM_Vj_tmax=np.array([max(SQM_Vj_r),max(SQM_Vj_y),max(SQM_Vj_g),max(SQM_Vj_b),max(SQM_Vj_w)])#,max(SQM_Vj_m)])
SQM_SQMA_tmax=np.array([max(SQM_SQMA_r),max(SQM_SQMA_y),max(SQM_SQMA_g),max(SQM_SQMA_b),max(SQM_SQMA_w)])#,max(SQM_SQMA_m)])
SQMASQM_tmax=np.array([max(SQMASQM_r),max(SQMASQM_y),max(SQMASQM_g),max(SQMASQM_b),max(SQMASQM_w)])#,max(SQMASQM_m)])
B_V_tmax=np.array([max(B_V_r),max(B_V_y),max(B_V_g),max(B_V_b),max(B_V_w)])#,max(B_V_m)])
SQM_SCA_tmax=np.array([max(SQM_SCA_r),max(SQM_SCA_y),max(SQM_SCA_g),max(SQM_SCA_b),max(SQM_SCA_w)])#,max(SQM_SCA_m)])
SQM_VlA_tmax=np.array([max(SQM_VlA_r),max(SQM_VlA_y),max(SQM_VlA_g),max(SQM_VlA_b),max(SQM_VlA_w)])#,max(SQM_VlA_m)])
SQM_VjA_tmax=np.array([max(SQM_VjA_r),max(SQM_VjA_y),max(SQM_VjA_g),max(SQM_VjA_b),max(SQM_VjA_w)])#,max(SQM_VjA_m)])
B_VA_tmax=np.array([max(B_VA_r),max(B_VA_y),max(B_VA_g),max(B_VA_b),max(B_VA_w)])#,max(B_VA_m)])
SQM_SC_tmin=np.array([min(SQM_SC_r),min(SQM_SC_y),min(SQM_SC_g),min(SQM_SC_b),min(SQM_SC_w)])#,min(SQM_SC_m)])
SQM_Vl_tmin=np.array([min(SQM_Vl_r),min(SQM_Vl_y),min(SQM_Vl_g),min(SQM_Vl_b),min(SQM_Vl_w)])#,min(SQM_Vl_m)])
SQM_Vj_tmin=np.array([min(SQM_Vj_r),min(SQM_Vj_y),min(SQM_Vj_g),min(SQM_Vj_b),min(SQM_Vj_w)])#,min(SQM_Vj_m)])
B_V_tmin=np.array([min(B_V_r),min(B_V_y),min(B_V_g),min(B_V_b),min(B_V_w)])#,min(B_V_m)])
SQM_SCA_tmin=np.array([min(SQM_SCA_r),min(SQM_SCA_y),min(SQM_SCA_g),min(SQM_SCA_b),min(SQM_SCA_w)])#,min(SQM_SCA_m)])
SQM_VlA_tmin=np.array([min(SQM_VlA_r),min(SQM_VlA_y),min(SQM_VlA_g),min(SQM_VlA_b),min(SQM_VlA_w)])#,min(SQM_VlA_m)])
SQM_VjA_tmin=np.array([min(SQM_VjA_r),min(SQM_VjA_y),min(SQM_VjA_g),min(SQM_VjA_b),min(SQM_VjA_w)])#,min(SQM_VjA_m)])
SQM_SQMA_tmin=np.array([min(SQM_SQMA_r),min(SQM_SQMA_y),min(SQM_SQMA_g),min(SQM_SQMA_b),min(SQM_SQMA_w)])#,min(SQM_SQMA_m)])
SQMASQM_tmin=np.array([min(SQMASQM_r),min(SQMASQM_y),min(SQMASQM_g),min(SQMASQM_b),min(SQMASQM_w)])#,min(SQM_SQMA_m)])
B_VA_tmin=np.array([min(B_VA_r),min(B_V_y),min(B_VA_g),min(B_VA_b),min(B_VA_w)])#,min(B_VA_m)])
wave_tmin=np.array([min(wave_r),min(wave_y),min(wave_g),min(wave_b),min(wave_w)])#,min(wave_m)])
wave_tmax=np.array([max(wave_r),max(wave_y),max(wave_g),max(wave_b),max(wave_w)])#,min(wave_m)])
SQM_SC_tmax=SQM_SC_tmax-SQM_SC_t
SQM_Vl_tmax=SQM_Vl_tmax-SQM_Vl_t
SQM_Vj_tmax=SQM_Vj_tmax-SQM_Vj_t
SQM_SQMA_tmax=SQM_SQMA_tmax-SQM_SQMA_t
B_V_tmax=B_V_tmax-B_V_t
SQM_SCA_tmax=SQM_SCA_tmax-SQM_SCA_t
SQM_VlA_tmax=SQM_VlA_tmax-SQM_VlA_t
SQM_VjA_tmax=SQM_VjA_tmax-SQM_VjA_t
B_VA_tmax=B_VA_tmax-B_VA_t
SQM_SC_tmin=SQM_SC_t-SQM_SC_tmin
SQM_Vl_tmin=SQM_Vl_t-SQM_Vl_tmin
SQM_Vj_tmin=SQM_Vj_t-SQM_Vj_tmin
SQM_SQMA_tmin=SQM_SQMA_t-SQM_SQMA_tmin
B_V_tmin=B_V_t-B_V_tmin
SQM_SCA_tmin=SQM_SCA_t-SQM_SCA_tmin
SQM_VlA_tmin=SQM_VlA_t-SQM_VlA_tmin
SQM_VjA_tmin=SQM_VjA_t-SQM_VjA_tmin
B_VA_tmin=B_VA_t-B_VA_tmin

patches = []
fig, ax = plt.subplots()
#ax.plot(BG,GR,40,marker=markers1)
#scatter3(BG,GR,spec1, 40,marker=markers1, alpha)

d = {'x': BG, 'y': GR, 'color': color1}
df = pd.DataFrame(d)
mm = df.groupby('color')
cat = ['r','c','b','y','k','g']
char = ['*','o','>','<','v','s']
char2= ['*','o','>','>','V','V','s','d','D','p']
tp =   ['HPS','INC','HAL','CMH','CFL','FL','MV','LPS','OTH','LED']
#scatter2(cat,char,mm,40,alpha)
#ax.scatter(BG,GR, c=c, cmap=cmap)
scatter3(BG,GR,spec1,40,alpha)
if 0:
 for X in range(len(listafiles)):
   try:
    plt.text(BG[X],GR[X]+0.01,str(spec[X].num),ha='right', va='bottom')
   except:
    print "fallo "+listafiles[X]
plt.xlabel("B/G")
plt.ylabel("G/R")
#plt.title("ISS colors for typical lamps")
plt.xlim(0,1.1)
plt.ylim(0,1.1)
rect = mpatches.Rectangle([0., 0.56], 0.34, 5, ec="none")#MV
patches.append(rect)
rect = mpatches.Rectangle([0.34, 0.56], 5, 5, ec="none")#MH
patches.append(rect)
rect = mpatches.Rectangle([0, 0.36], 0.25,0.2 , ec="none")#White sodioum
patches.append(rect)
rect = mpatches.Rectangle([0.00, 0.], 0.05,0.36 , ec="none")#LPS
patches.append(rect)
rect = mpatches.Rectangle([0.05, 0.00], 0.20,0.36 , ec="none")#HPS
patches.append(rect)
rect = mpatches.Rectangle([0.25, 0.00], 0.2,0.56 , ec="none")#pure White
patches.append(rect)
rect = mpatches.Rectangle([0.45, 0.00], 5,0.56 , ec="none")#No
patches.append(rect)
colors = np.linspace(0, 1, len(patches))
collection = PatchCollection(patches, cmap=plt.cm.hsv, alpha=0.1)
collection.set_array(np.array(colors))
ax.add_collection(collection)
plt.savefig('colores_espectros2.png')
plt.show()



asciitable.write(spec,'espectros.txt')
asciitable.write(spec2,'espectros2.txt')




model_ransac_VlG_BG,line_y_ransac,line_X=ransacfit(BG,VlG)
zVlG_BG = np.polyfit(BG, VlG, 3)
p = np.poly1d(zVlG_BG)
xp = np.linspace(0, 1, 100)
plt.plot(line_X, line_y_ransac, '-b')#, label='RANSAC regressor')
plt.plot(xp, p(xp),'r-')
d = {'x': BG, 'y': VlG, 'color': color1}
df = pd.DataFrame(d)
mm = df.groupby('color')
cat = ['r','c','b','y','k','g']
char = ['*','o','>','<','v','s']
scatter2(cat,char,mm,40,alpha)
#scatter3(BG,VlG,spec1, 40, alpha)
if 0:
 for X in range(len(listafiles)):
   try:
    plt.text(BG[X],VlG[X]+0.01,str(spec[X].num),ha='right', va='bottom')
   except:
    print "fallo"
plt.xlabel("B/G")
plt.ylabel("Vl/G")
#plt.title("ISS colors for typical lamps")
plt.xlim(0,0.9)
plt.ylim(0,3)
plt.savefig('colores_espectros_VlG.png')
plt.show()

model_ransac_SCG_GR,line_y_ransac,line_X=ransacfit(GR,SCG)
zSCGGR = np.polyfit(GR, SCG, 5)
p = np.poly1d(zSCGGR)
xp = np.linspace(0, 2, 100)
plt.plot(line_X, line_y_ransac, '-b')#, label='RANSAC regressor')
d = {'x': GR, 'y': SCG, 'color': color1}
df = pd.DataFrame(d)
mm = df.groupby('color')
cat = ['r','c','b','y','k','g']
char = ['*','o','>','<','v','s']
scatter2(cat,char,mm,40,alpha)
#scatter3(GR,SCG,spec1, 40, alpha)
plt.plot(xp, p(xp),'r-')
if 0:
 for X in range(len(listafiles)):
   try:
    plt.text(GR[X],SCG[X]+0.01,str(spec[X].num),ha='right', va='bottom')
   except:
    print "fallo"
plt.xlabel("G/R")
plt.ylabel("SC/G")
plt.xlim(0,1.2)
plt.ylim(0,1)
#plt.title("ISS colors for typical lamps")
#plt.xlim(0,0.9)
#plt.ylim(0,2.1)
plt.savefig('colores_espectros_SCG_GR.png')
plt.show()

model_ransac_VlG_GR,line_y_ransac,line_X=ransacfit(GR,VlG)
zVlGGR = np.polyfit(GR, VlG, 3)
p = np.poly1d(zVlGGR)
xp = np.linspace(0, 2, 100)
plt.plot(line_X, line_y_ransac, '-b')#, label='RANSAC regressor')
d = {'x': GR, 'y': VlG, 'color': color1}
df = pd.DataFrame(d)
mm = df.groupby('color')
cat = ['r','c','b','y','k','g']
char = ['*','o','>','<','v','s']
scatter2(cat,char,mm,40,alpha)
#scatter3(GR,VlG,spec1, 40, alpha)
plt.plot(xp, p(xp),'r-')
if 0:
 for X in range(len(listafiles)):
   try:
    plt.text(GR[X],VlG[X]+0.01,str(spec[X].num),ha='right', va='bottom')
   except:
    print "fallo"
plt.xlabel("G/R")
plt.ylabel("Vl/G")
plt.xlim(0,1.2)
plt.ylim(0.5,3)
#plt.title("ISS colors for typical lamps")
#plt.xlim(0,0.9)
#plt.ylim(0,2.1)
plt.savefig('colores_espectros_SCG_GR.png')
plt.show()

model_ransac_SQM_G_GR,line_y_ransac,line_X=ransacfit(GR,SQM_G)
zSQM_GGR = np.polyfit(GR, SQM_G, 3)
p = np.poly1d(zSQM_GGR)
xp = np.linspace(0, 2, 100)
plt.plot(line_X, line_y_ransac, '-b')#, label='RANSAC regressor')
d = {'x': GR, 'y': SQM_G, 'color': color1}
df = pd.DataFrame(d)
mm = df.groupby('color')
cat = ['r','c','b','y','k','g']
char = ['*','o','>','<','v','s']
scatter2(cat,char,mm,40,alpha)
#scatter3(GR,SQM_G,spec1, 40, alpha)
plt.plot(xp, p(xp),'r-')
if 0:
 for X in range(len(listafiles)):
   try:
    plt.text(GR[X],SQM_G[X]+0.01,str(spec[X].num),ha='right', va='bottom')
   except:
    print "fallo"
plt.xlabel("G/R")
plt.ylabel("SQM/G")
plt.xlim(0,1.2)
plt.ylim(0.5,1.3)
#plt.title("ISS colors for typical lamps")
#plt.xlim(0,0.9)
#plt.ylim(0,2.1)
plt.savefig('colores_espectros_SQM_G_GR.png')
plt.show()

model_ransac_SQM_G_BG,line_y_ransac,line_X=ransacfit(BG,SQM_G)
zSQM_GBG = np.polyfit(BG, SQM_G, 3)
p = np.poly1d(zSQM_GBG)
xp = np.linspace(0, 2, 100)
plt.plot(line_X, line_y_ransac, '-b')#, label='RANSAC regressor')
d = {'x': BG, 'y': SQM_G, 'color': color1}
df = pd.DataFrame(d)
mm = df.groupby('color')
cat = ['r','c','b','y','k','g']
char = ['*','o','>','<','v','s']
scatter2(cat,char,mm,40,alpha)
#scatter3(BG,SQM_G,spec1, 40, alpha)
plt.plot(xp, p(xp),'r-')
if 0:
 for X in range(len(listafiles)):
   try:
    plt.text(BG[X],SQM_G[X]+0.01,str(spec[X].num),ha='right', va='bottom')
   except:
    print "fallo"
plt.xlabel("B/G")
plt.ylabel("SQM/G")
plt.xlim(0,1.2)
plt.ylim(0.5,1.3)
#plt.title("ISS colors for typical lamps")
#plt.xlim(0,0.9)
#plt.ylim(0,2.1)
plt.savefig('colores_espectros_SQM_G_BG.png')
plt.show()

#model_ransac_SQM_G_BGj,line_y_ransac,line_X=ransacfit(spec.Bj_Vj,SQM_G)
#zSQM_GBGj = np.polyfit(spec.Bj_Vj, SQM_G, 3)
#p = np.poly1d(zSQM_GBGj)
#xp = np.linspace(0, 2, 100)
plt.plot(line_X, line_y_ransac, '-b')#, label='RANSAC regressor')
d = {'x': spec.Bj_Vj, 'y': SQM_G, 'color': color1}
df = pd.DataFrame(d)
mm = df.groupby('color')
cat = ['r','c','b','y','k','g']
char = ['*','o','>','<','v','s']
scatter2(cat,char,mm,40,alpha)
#scatter3(spec.Bj_Vj,SQM_G,spec1, 40, alpha)
plt.plot(xp, p(xp),'r-')
if 0:
 for X in range(len(listafiles)):
   try:
    plt.text(spec.Bj_Vj[X],SQM_G[X]+0.01,str(spec[X].num),ha='right', va='bottom')
   except:
    print "fallo"
plt.xlabel("Bj - Vj")
plt.ylabel("SQM/G")
plt.xlim(0,5)
plt.ylim(0.5,1.3)
#plt.title("ISS colors for typical lamps")
#plt.xlim(0,0.9)
#plt.ylim(0,2.1)
plt.savefig('colores_espectros_SQM_G_BGj.png')
plt.show()

fig, ax1 = plt.subplots()

#model_ransac_SQM_Vj_BGj,line_y_ransac,line_X=ransacfit(spec.Bj_Vj,SQM_Vj)
#zSQM_VjBGj = np.polyfit(spec.Bj_Vj, SQM_Vj, 3)
#p = np.poly1d(zSQM_VjBGj)
#xp = np.linspace(0, 2, 100)
#plt.plot(line_X, line_y_ransac, '-b')#, label='RANSAC regressor')
d = {'x': spec.Bj_Vj, 'y': SQM_Vj, 'color': color1}
df = pd.DataFrame(d)
mm = df.groupby('color')
cat = ['r','c','b','y','k','g']
char = ['*','o','>','<','v','s']
scatter2(cat,char,mm,40,alpha)
#scatter3(spec.Bj_Vj,SQM_Vj,spec1, 40, alpha)
ax1.errorbar(B_V_t, SQM_Vj_t,ecolor='k',capthick=2, yerr=[SQM_Vj_tmin, SQM_Vj_tmax], xerr=[B_V_tmin, B_V_tmax], fmt='.')
#scatter3(B_V_t, SQM_Vj_t,s=500, c=['r','y','g','b','c'], alpha)

d = {'x': B_V_t, 'y':  SQM_Vj_t, 'color': ['r','y','g','b','c']}
df = pd.DataFrame(d)
mm = df.groupby('color')
cat = ['r','y','g','b','c']
char = ['*','v','s','>','o']
scatter2(cat,char,mm,500,0.5)

#plt.plot(xp, p(xp),'r-')
if 0:
 for X in range(len(listafiles)):
   try:
    ax1.text(spec.Bj_Vj[X],SQM_Vj[X]+0.01,str(spec[X].num),ha='right', va='bottom')
   except:
    print "fallo"
ax1.set_xlabel("Bj - Vj")
ax1.set_ylabel("SQM - Vj")
ax1.set_xlim(0,7)
ax1.set_ylim(0.2,1)
#ax1.set_title("SQM - V vs B-V Johnson color")
ax2 = ax1.twinx()


ax2.set_ylim(0.2-0.35,1-0.35)
ax2.set_ylabel("SQM* - Vj")
#plt.xlim(0,0.9)
#plt.ylim(0,2.1)
fig.savefig('colores_espectros_SQM_Vj_BGj.png')
plt.show()

########################################
#With clouds SQM -V vs B-V
########################################

#model_ransac_SQM_VjA_BGAj,line_y_ransac,line_X=ransacfit(spec.Bj_VjA,SQM_VjA)
#zSQM_VjABGjA = np.polyfit(spec.Bj_VjA, SQM_VjA, 3)
#p = np.poly1d(zSQM_VjABGjA)
#xp = np.linspace(0, 2, 100)
#plt.plot(line_X, line_y_ransac, '-b')#, label='RANSAC regressor')
d = {'x': spec.Bj_VjA, 'y': SQM_VjA, 'color': color1}
df = pd.DataFrame(d)
mm = df.groupby('color')
cat = ['r','c','b','y','k','g']
char = ['*','o','>','<','v','s']
scatter2(cat,char,mm,40,alpha)
#scatter3(spec.Bj_VjA,SQM_VjA,spec1, 40, alpha)
plt.errorbar(B_VA_t, SQM_VjA_t,ecolor='k',capthick=2, yerr=[SQM_VjA_tmin, SQM_VjA_tmax], xerr=[B_VA_tmin, B_VA_tmax], fmt='.')
#scatter3(B_VA_t, SQM_VjA_t,s=500, c=['r','y','g','b','c'], alpha)

d = {'x': B_VA_t, 'y':  SQM_VjA_t, 'color': ['r','y','g','b','c']}
df = pd.DataFrame(d)
mm = df.groupby('color')
cat = ['r','y','g','b','c']
char = ['*','v','s','>','o']
scatter2(cat,char,mm,500,0.5)

#plt.plot(xp, p(xp),'r-')
if 0:
 for X in range(len(listafiles)):
   try:
    plt.text(spec.Bj_VjA[X],SQM_VjA[X]+0.01,str(spec[X].num),ha='right', va='bottom')
   except:
    print "fallo"
plt.xlabel("Bj - Vj")
plt.ylabel("SQM - Vj")
plt.xlim(0,7)
plt.ylim(0.2,1)
#plt.title("SQM - V vs B-V Johnson color Cloud amplified")
#plt.xlim(0,0.9)
#plt.ylim(0,2.1)
plt.savefig('colores_espectros_SQM_Vj_BGjA_'+clouds[:-4]+'.png')
plt.show()


########################################

########################################
#model_ransac_SQMVl_BGj,line_y_ransac,line_X=ransacfit(spec.Bj_Vj,SQMVl)
#zSQMVlBGj = np.polyfit(spec.Bj_Vj, SQMVl, 3)
#p = np.poly1d(zSQMVlBGj)
#xp = np.linspace(0, 2, 100)
#plt.plot(line_X, line_y_ransac, '-b')#, label='RANSAC regressor')
d = {'x': spec.Bj_Vj, 'y': SQMVl, 'color': color1}
df = pd.DataFrame(d)
mm = df.groupby('color')
cat = ['r','c','b','y','k','g']
char = ['*','o','>','<','v','s']
scatter2(cat,char,mm,40,alpha)
#scatter3(spec.Bj_Vj,SQMVl,spec1, 40, alpha)
#plt.plot(xp, p(xp),'r-')
if 1:
 for X in range(len(listafiles)):
   try:
    plt.text(spec.Bj_Vj[X],SQMVl[X]+0.01,str(spec[X].tipo),ha='right', va='bottom')
   except:
    print "fallo"
plt.xlabel("Bj - Vj")
plt.ylabel("SQM/Vl")
plt.xlim(0,7)
#plt.ylim(0.5,1.3)
#plt.title("SQM/Photopic vs B-V Johnson")
#plt.xlim(0,0.9)
#plt.ylim(0,2.1)
plt.savefig('colores_espectros_SQMVl_BGj.png')
plt.show()


###########################################

###########################################

#model_ransac_SQMSC_BGj,line_y_ransac,line_X=ransacfit(spec.Bj_Vj,SQMSC)
#zSQMSCBGj = np.polyfit(spec.Bj_Vj, SQMSC, 3)
#p = np.poly1d(zSQMSCBGj)
#xp = np.linspace(0, 2, 100)
#plt.plot(line_X, line_y_ransac, '-b')#, label='RANSAC regressor')
#scatter3(spec.Bj_Vj,SQMSC,spec1, 40, alpha)
d = {'x': spec.Bj_Vj, 'y': SQMSC, 'color': color1}
df = pd.DataFrame(d)
mm = df.groupby('color')
cat = ['r','c','b','y','k','g']
char = ['*','o','>','<','v','s']
scatter2(cat,char,mm,40,alpha)
#plt.plot(xp, p(xp),'r-')
if 0:
 for X in range(len(listafiles)):
   try:
    plt.text(spec.Bj_Vj[X],SQMSC[X]+0.01,str(spec[X].num),ha='right', va='bottom')
   except:
    print "fallo"
plt.xlabel("Bj - Vj")
plt.ylabel("SQM/SC")
plt.xlim(0,7)
#plt.ylim(0.5,1.3)
#plt.title("SQM/Scotopic vs B-V Johnson color")
#plt.xlim(0,0.9)
#plt.ylim(0,2.1)
plt.savefig('colores_espectros_SQMSC_BGj.png')
plt.show()


###########################################

###########################################

########################################

########################################
#model_ransac_SQM_Vl_BGj,line_y_ransac,line_X=ransacfit(spec.Bj_Vj,SQM_Vl)
#zSQM_VlBGj = np.polyfit(spec.Bj_Vj, SQM_Vl, 3)
#p = np.poly1d(zSQM_VlBGj)
#xp = np.linspace(0, 2, 100)
#plt.plot(line_X, line_y_ransac, '-b')#, label='RANSAC regressor')
fig, ax1 = plt.subplots()
d = {'x': spec.Bj_Vj, 'y': SQM_Vl, 'color': color1}
df = pd.DataFrame(d)
mm = df.groupby('color')
cat = ['r','c','b','y','k','g']
char = ['*','o','>','<','v','s']
scatter2(cat,char,mm,40,alpha)
#scatter3(spec.Bj_Vj,SQM_Vl,spec1, 40, alpha)
ax1.errorbar(B_V_t, SQM_Vl_t,ecolor='k',capthick=2, yerr=[SQM_Vl_tmin, SQM_Vl_tmax], xerr=[B_V_tmin, B_V_tmax], fmt='.')
#scatter3(B_V_t, SQM_Vl_t,s=500, c=['r','y','g','b','c'], alpha)


d = {'x': B_V_t, 'y':  SQM_Vl_t, 'color': ['r','y','g','b','c']}
df = pd.DataFrame(d)
mm = df.groupby('color')
cat = ['r','y','g','b','c']
char = ['*','v','s','>','o']
scatter2(cat,char,mm,500,0.5)

ax2 = ax1.twinx()


#plt.plot(xp, p(xp),'r-')
if 0:
 for X in range(len(listafiles)):
   try:
    plt.text(spec.Bj_Vj[X],SQM_Vl[X]+0.01,str(spec[X].num),ha='right', va='bottom')
   except:
    print "fallo"
ax1.set_xlabel("Bj - Vj")
ax1.set_ylabel("SQM - V($\lambda$) Photopic vision")
ax1.set_xlim(0,7)
ax1.set_ylim(0.2,1.1)
#ax1.set_title("SQM-Photopic vs B-V Johnson color")
ax2.set_ylim(0.2-0.35,1-0.35)
ax2.set_ylabel("SQM* - V($\lambda$) Photopic vision")
#plt.xlim(0,0.9)
#plt.ylim(0,2.1)
plt.savefig('colores_espectros_SQM_Vl_BVj.png')
plt.show()
########################################
#Clouds
########################################
#model_ransac_SQM_Vl_BGj,line_y_ransac,line_X=ransacfit(spec.Bj_VjA,SQM_VlA)
#zSQM_VlBGjA = np.polyfit(spec.Bj_VjA, SQM_VlA, 3)
#p = np.poly1d(zSQM_VlBGjA)
#xp = np.linspace(0, 2, 100)
#plt.plot(line_X, line_y_ransac, '-b')#, label='RANSAC regressor')
d = {'x': spec.Bj_VjA, 'y': SQM_VlA, 'color': color1}
df = pd.DataFrame(d)
mm = df.groupby('color')
cat = ['r','c','b','y','k','g']
char = ['*','o','>','<','v','s']
scatter2(cat,char,mm,40,alpha)
#scatter3(spec.Bj_VjA,SQM_VlA,spec1, 40, alpha)
plt.errorbar(B_VA_t, SQM_VlA_t,ecolor='k',capthick=2, yerr=[SQM_VlA_tmin, SQM_VlA_tmax], xerr=[B_VA_tmin, B_VA_tmax], fmt='.')
#scatter3(B_VA_t, SQM_VlA_t,s=500, c=['r','y','g','b','c'], alpha)

d = {'x': B_VA_t, 'y':  SQM_VlA_t, 'color': ['r','y','g','b','c']}
df = pd.DataFrame(d)
mm = df.groupby('color')
cat = ['r','y','g','b','c']
char = ['*','v','s','>','o']
scatter2(cat,char,mm,500,0.5)

#plt.plot(xp, p(xp),'r-')
if 0:
 for X in range(len(listafiles)):
   try:
    plt.text(spec.Bj_VjA[X],SQM_VlA[X]+0.01,str(spec[X].num),ha='right', va='bottom')
   except:
    print "fallo"
plt.xlabel("Bj - Vj")
plt.ylabel("SQM - V($\lambda$) Photopic vision")
plt.xlim(0,7)
plt.ylim(0.2,1.1)
#plt.title("SQM-Photopic vs B-V Johnson cloud amplified")
#plt.xlim(0,0.9)
#plt.ylim(0,2.1)
plt.savefig('colores_espectros_SQM_Vl_BVjA.png')
plt.show()

###########################################

###########################################
if 0:
    model_ransac_SQM_SC_BGj,line_y_ransac,line_X=ransacfit(spec.Bj_Vj,SQM_SC)
    zSQM_SCBGj = np.polyfit(spec.Bj_Vj, SQM_SC, 3)
    p = np.poly1d(zSQM_SCBGj)
    xp = np.linspace(0, 2, 100)
    #plt.plot(line_X, line_y_ransac, '-b')#, label='RANSAC regressor')
#scatter3(spec.Bj_Vj,SQM_SC,spec1, 40, alpha)
d = {'x': spec.Bj_Vj, 'y': SQM_SC, 'color': color1}
df = pd.DataFrame(d)
mm = df.groupby('color')
cat = ['r','c','b','y','k','g']
char = ['*','o','>','<','v','s']
scatter2(cat,char,mm,40,alpha)
plt.errorbar(B_V_t, SQM_SC_t,ecolor='k',capthick=2, yerr=[SQM_SC_tmin, SQM_SC_tmax], xerr=[B_V_tmin, B_V_tmax], fmt='.')
#scatter3(B_V_t, SQM_SC_t,s=500, c=['r','y','g','b','c'], alpha)

d = {'x': B_V_t, 'y':  SQM_SC_t, 'color': ['r','y','g','b','c']}
df = pd.DataFrame(d)
mm = df.groupby('color')
cat = ['r','y','g','b','c']
char = ['*','v','s','>','o']
scatter2(cat,char,mm,500,0.5)

#plt.plot(xp, p(xp),'r-')
if 0:
 for X in range(len(listafiles)):
   try:
    plt.text(spec.Bj_Vj[X],SQM_SC[X]+0.01,str(spec[X].num),ha='right', va='bottom')
   except:
    print "fallo"
plt.xlabel("Bj - Vj")
plt.ylabel("SQM - V'($\lambda$) Soctopic vision")
plt.xlim(0,7)
plt.ylim(-2,0.4)
#plt.title("SQM-Scotopic vs B-V Johnson color")
#plt.xlim(0,0.9)
#plt.ylim(0,2.1)
plt.savefig('colores_espectros_SQM_SC_BGj.png')
plt.show()

###########################################
# Clouds
###########################################
if 0:
    model_ransac_SQM_SC_BGj,line_y_ransac,line_X=ransacfit(spec.Bj_VjA,SQM_SCA)
    zSQM_SCBGjA = np.polyfit(spec.Bj_VjA, SQM_SCA, 3)
    p = np.poly1d(zSQM_SCBGjA)
    xp = np.linspace(0, 2, 100)
    #plt.plot(line_X, line_y_ransac, '-b')#, label='RANSAC regressor')
#scatter3(spec.Bj_VjA,SQM_SCA,spec1, 40, alpha)
d = {'x': spec.Bj_VjA, 'y': SQM_SCA, 'color': color1}
df = pd.DataFrame(d)
mm = df.groupby('color')
cat = ['r','c','b','y','k','g']
char = ['*','o','>','<','v','s']
scatter2(cat,char,mm,40,alpha)
plt.errorbar(B_VA_t, SQM_SCA_t,ecolor='k',capthick=2, yerr=[SQM_SCA_tmin, SQM_SCA_tmax], xerr=[B_VA_tmin, B_VA_tmax], fmt='.')
#scatter3(B_VA_t, SQM_SCA_t,s=500, c=['r','y','g','b','c'], alpha)


d = {'x': B_VA_t, 'y':  SQM_SCA_t, 'color': ['r','y','g','b','c']}
df = pd.DataFrame(d)
mm = df.groupby('color')
cat = ['r','y','g','b','c']
char = ['*','v','s','>','o']
scatter2(cat,char,mm,500,0.5)

#plt.plot(xp, p(xp),'r-')
if 0:
 for X in range(len(listafiles)):
   try:
    plt.text(spec.Bj_Vj[X],SQM_SC[X]+0.01,str(spec[X].num),ha='right', va='bottom')
   except:
    print "fallo"
plt.xlabel("Bj - Vj")
plt.ylabel("SQM - V'($\lambda$) Soctopic vision")
plt.xlim(0,7)
plt.ylim(-2,0.4)
#plt.title("SQM-Scotopic vs B-V Johnson cloud amplified")
#plt.xlim(0,0.9)
#plt.ylim(0,2.1)
plt.savefig('colores_espectros_SQM_SC_BGjA.png')
plt.show()
###########################################
# Clouds
###########################################
if 0:
    model_ransac_SQM_SQMA_BGj,line_y_ransac,line_X=ransacfit(spec.Bj_VjA,SQM_SQMA)
    zSQM_SQMABGj = np.polyfit(spec.Bj_Vj, SQM_SQMA, 3)
    p = np.poly1d(zSQM_SQMABGj)
    xp = np.linspace(0, 2, 100)
    #plt.plot(line_X, line_y_ransac, '-b')#, label='RANSAC regressor')
#scatter3(spec.Bj_Vj,SQM_SQMA,spec1, 40, alpha)
d = {'x': spec.Bj_Vj, 'y': SQM_SQMA, 'color': color1}
df = pd.DataFrame(d)
mm = df.groupby('color')
cat = ['r','c','b','y','k','g']
char = ['*','o','>','<','v','s']
scatter2(cat,char,mm,40,alpha)
plt.errorbar(B_V_t, SQM_SQMA_t,ecolor='k',capthick=2, yerr=[SQM_SQMA_tmin, SQM_SQMA_tmax], xerr=[B_V_tmin, B_V_tmax], fmt='.')
#scatter3(B_V_t, SQM_SQMA_t,s=500, c=['r','y','g','b','c'], alpha)

d = {'x': B_V_t, 'y':  SQM_SQMA_t, 'color': ['r','y','g','b','c']}
df = pd.DataFrame(d)
mm = df.groupby('color')
cat = ['r','y','g','b','c']
char = ['*','v','s','>','o']
scatter2(cat,char,mm,500,0.5)

#plt.plot(xp, p(xp),'r-')
if 0:
 for X in range(len(listafiles)):
   try:
    plt.text(spec.Bj_Vj[X],SQM_SQMA[X]+0.01,str(spec[X].num),ha='right', va='bottom')
   except:
    print "fallo"
plt.xlabel("Bj - Vj")
plt.ylabel("SQM' - SQM")
plt.xlim(0,7)
#plt.ylim(-2,0.4)
#plt.title("SQM - SQM & clouds* vs B-V Johnson")
#plt.xlim(0,0.9)
#plt.ylim(0,2.1)
plt.savefig('colores_espectros_SQM_SQMA_BGjA'+clouds[:-4]+'.png')
plt.show()

###########################################
# Clouds wave
###########################################
fig, ax1 = plt.subplots()

if 0:
    model_ransac_SQM_SQMA_BGj,line_y_ransac,line_X=ransacfit(spec.Bj_VjA,SQM_SQMA)
    zSQM_SQMABGj = np.polyfit(spec.Bj_Vj, SQM_SQMA, 3)
    p = np.poly1d(zSQM_SQMABGj)
    xp = np.linspace(0, 2, 100)
    #plt.plot(line_X, line_y_ransac, '-b')#, label='RANSAC regressor')
d = {'x': spec.wave/10, 'y': SQMASQM, 'color': color1}
df = pd.DataFrame(d)
mm = df.groupby('color')
cat = ['r','c','b','y','k','g']
char = ['*','o','>','<','v','s']
scatter2(cat,char,mm,40,alpha)
#scatter3(spec.wave,SQM_SQMA,spec1, 40, alpha)
plt.errorbar(wave_t/10, SQMASQM_t,ecolor='k',capthick=2, yerr=[SQMASQM_t-SQMASQM_tmin, SQMASQM_tmax-SQMASQM_t], xerr=[wave_t/10 -wave_tmin/10, wave_tmax/10-wave_t/10], fmt='.')

d = {'x': wave_t/10, 'y':  SQMASQM_t, 'color': ['r','y','g','b','c']}
df = pd.DataFrame(d)
mm = df.groupby('color')
cat = ['r','y','g','b','c']
char = ['*','v','s','>','o']
scatter2(cat,char,mm,500,0.5)

#scatter3(wave_t/10, SQMASQM_t,s=500, c=['r','y','g','b','c'], alpha)
#plt.plot(xp, p(xp),'r-')
if 0:
 for X in range(len(listafiles)):
   try:
    plt.text(spec.Bj_Vj[X],SQM_SQMA[X]+0.01,str(spec[X].num),ha='right', va='bottom')
   except:
    print "fallo"

ax1.set_xlabel("Effective wavelength (nm)")
ax1.set_ylabel("SQM'/SQM")
ax2 = ax1.twinx()

ax2.set_ylabel("Amplification")
ax1.set_xlim(500,650)
ax2.set_xlim(500,650)
ax1.set_ylim(15,30)
ax2.set_ylim(15,30)
#plt.ylim(-2,0.4)
##plt.title("SQM/SQM & clouds* vs Avg wavelenght")
specamply=asciitable.read(clouds)
if 0:
    specamply.col1*SQM_LICA
    amp = InterpolatedUnivariateSpline(specamply.col1, specamply.col2, k=3)
    SQM_f = InterpolatedUnivariateSpline(SQM_LICA.col1, SQM_LICA.col2, k=3)
    ampf=amp(spectro.col2)
    SQMff=SQM_f(spectro.col2)
    SQMampf=SQMff*ampf
    AB_r =  InterpolatedUnivariateSpline(AB.wave, AB.flux, k=3)
    ABf=AB_r(spectro.col2)
ax2.plot(specamply.col1/10,specamply.col2)
#plt.xlim(0,0.9)
#plt.ylim(0,2.1)
plt.savefig('colores_espectros_SQM_SQMA_wave'+clouds[:-4]+'.png')
plt.show()


###################################################
####################################################
#######################################################3
#####################################################

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()


#model_ransac_SQM_Vj_BGj,line_y_ransac,line_X=ransacfit(spec.Bj_Vj,SQM_Vj)
#zSQM_VjBGj = np.polyfit(spec.Bj_Vj, SQM_Vj, 3)
#p = np.poly1d(zSQM_VjBGj)
#xp = np.linspace(0, 2, 100)
#plt.plot(line_X, line_y_ransac, '-b')#, label='RANSAC regressor')
plt.scatter(spec2.Bj_Vj,SQM_Vj2,marker='H', s=100, alpha=alpha)
if 0:
    d = {'x': spec2.Bj_Vj, 'y':SQM_Vj2, 'color': color1}
    df = pd.DataFrame(d)
    mm = df.groupby('color')
    cat = ['r','c','b','y','k','g']
    char = ['*','o','>','<','v','s']
    scatter2(cat,char,mm,40,alpha)
#ax1.errorbar(B_V_t, SQM_Vj_t,ecolor='k',capthick=2, yerr=[SQM_Vj_tmin, SQM_Vj_tmax], xerr=[B_V_tmin, B_V_tmax], fmt='.')
#scatter3(B_V_t, SQM_Vj_t,s=500, c=['r','y','g','b','c'], alpha)
#plt.plot(xp, p(xp),'r-')
if 1:
 for X in range(len(listafiles2)):
   try:
    ax1.text(spec2.Bj_Vj[X],SQM_Vj2[X]+0.01,str(spec2[X].num),ha='right', va='bottom')
   except:
    print "fallo2"
ax1.set_xlabel("Bj - Vj")
ax1.set_ylabel("SQM - Vj")
ax1.set_xlim(0,2)
ax1.set_ylim(0,0.6)
#ax1.set_title("Natural sources")

ax2.set_ylim(0-0.35,0.6-0.35)
ax2.set_ylabel("SQM* - Vj")
#plt.xlim(0,0.9)
#plt.ylim(0,2.1)
fig.savefig('colores_espectros_SQM_Vj_BGj2.png')
plt.show()
'''
########################################
#With clouds SQM -V vs B-V
########################################

#model_ransac_SQM_VjA_BGAj,line_y_ransac,line_X=ransacfit(spec.Bj_VjA,SQM_VjA)
#zSQM_VjABGjA = np.polyfit(spec.Bj_VjA, SQM_VjA, 3)
#p = np.poly1d(zSQM_VjABGjA)
#xp = np.linspace(0, 2, 100)
#plt.plot(line_X, line_y_ransac, '-b')#, label='RANSAC regressor')
#scatter3(spec.Bj_VjA,SQM_VjA,spec1, 40, alpha)

d = {'x': spec.Bj_VjA, 'y': SQM_VjA, 'color': color1}
df = pd.DataFrame(d)
mm = df.groupby('color')
cat = ['r','c','b','y','k','g']
char = ['*','o','>','<','v','s']
scatter2(cat,char,mm,40,alpha)
plt.errorbar(B_VA_t, SQM_VjA_t,ecolor='k',capthick=2, yerr=[SQM_VjA_tmin, SQM_VjA_tmax], xerr=[B_VA_tmin, B_VA_tmax], fmt='.')
scatter3(B_VA_t, SQM_VjA_t,s=500, c=['r','y','g','b','c'], alpha)
#plt.plot(xp, p(xp),'r-')
if 0:
 for X in range(len(listafiles)):
   try:
    plt.text(spec.Bj_VjA[X],SQM_VjA[X]+0.01,str(spec[X].num),ha='right', va='bottom')
   except:
    print "fallo"
plt.xlabel("Bj - Vj")
plt.ylabel("SQM - Vj")
plt.xlim(0,7)
plt.ylim(0.2,1)
#plt.title("SQM - V vs B-V Johnson color Cloud amplified")
#plt.xlim(0,0.9)
#plt.ylim(0,2.1)
plt.savefig('colores_espectros_SQM_Vj_BGjA.png')
plt.show()


########################################

########################################
#model_ransac_SQMVl_BGj,line_y_ransac,line_X=ransacfit(spec.Bj_Vj,SQMVl)
#zSQMVlBGj = np.polyfit(spec.Bj_Vj, SQMVl, 3)
#p = np.poly1d(zSQMVlBGj)
#xp = np.linspace(0, 2, 100)
#plt.plot(line_X, line_y_ransac, '-b')#, label='RANSAC regressor')
#scatter3(spec.Bj_Vj,SQMVl,spec1, 40, alpha)
d = {'x': spec.Bj_Vj, 'y': SQMVl, 'color': color1}
df = pd.DataFrame(d)
mm = df.groupby('color')
cat = ['r','c','b','y','k','g']
char = ['*','o','>','<','v','s']
scatter2(cat,char,mm,40,alpha)
#plt.plot(xp, p(xp),'r-')
if 0:
 for X in range(len(listafiles)):
   try:
    plt.text(spec.Bj_Vj[X],SQMVl[X]+0.01,str(spec[X].num),ha='right', va='bottom')
   except:
    print "fallo"
plt.xlabel("Bj - Vj")
plt.ylabel("SQM/Vl")
plt.xlim(0,7)
#plt.ylim(0.5,1.3)
#plt.title("SQM/Photopic vs B-V Johnson")
#plt.xlim(0,0.9)
#plt.ylim(0,2.1)
plt.savefig('colores_espectros_SQMVl_BGj.png')
plt.show()


###########################################

###########################################

#model_ransac_SQMSC_BGj,line_y_ransac,line_X=ransacfit(spec.Bj_Vj,SQMSC)
#zSQMSCBGj = np.polyfit(spec.Bj_Vj, SQMSC, 3)
#p = np.poly1d(zSQMSCBGj)
#xp = np.linspace(0, 2, 100)
#plt.plot(line_X, line_y_ransac, '-b')#, label='RANSAC regressor')
#scatter3(spec.Bj_Vj,SQMSC,spec1, 40, alpha)

d = {'x': spec.Bj_Vj, 'y': SQMSC, 'color': color1}
df = pd.DataFrame(d)
mm = df.groupby('color')
cat = ['r','c','b','y','k','g']
char = ['*','o','>','<','v','s']
scatter2(cat,char,mm,40,alpha)
#plt.plot(xp, p(xp),'r-')
if 0:
 for X in range(len(listafiles)):
   try:
    plt.text(spec.Bj_Vj[X],SQMSC[X]+0.01,str(spec[X].num),ha='right', va='bottom')
   except:
    print "fallo"
plt.xlabel("Bj - Vj")
plt.ylabel("SQM/SC")
plt.xlim(0,7)
#plt.ylim(0.5,1.3)
#plt.title("SQM/Scotopic vs B-V Johnson color")
#plt.xlim(0,0.9)
#plt.ylim(0,2.1)
plt.savefig('colores_espectros_SQMSC_BGj.png')
plt.show()


###########################################

###########################################

########################################

########################################
#model_ransac_SQM_Vl_BGj,line_y_ransac,line_X=ransacfit(spec.Bj_Vj,SQM_Vl)
#zSQM_VlBGj = np.polyfit(spec.Bj_Vj, SQM_Vl, 3)
#p = np.poly1d(zSQM_VlBGj)
#xp = np.linspace(0, 2, 100)
#plt.plot(line_X, line_y_ransac, '-b')#, label='RANSAC regressor')
fig, ax1 = plt.subplots()
#scatter3(spec.Bj_Vj,SQM_Vl,spec1, 40, alpha)

d = {'x': spec.Bj_Vj, 'y': SQM_Vl, 'color': color1}
df = pd.DataFrame(d)
mm = df.groupby('color')
cat = ['r','c','b','y','k','g']
char = ['*','o','>','<','v','s']
scatter2(cat,char,mm,40,alpha)

ax1.errorbar(B_V_t, SQM_Vl_t,ecolor='k',capthick=2, yerr=[SQM_Vl_tmin, SQM_Vl_tmax], xerr=[B_V_tmin, B_V_tmax], fmt='.')
scatter3(B_V_t, SQM_Vl_t,s=500, c=['r','y','g','b','c'], alpha)
ax2 = ax1.twinx()


#plt.plot(xp, p(xp),'r-')
if 0:
 for X in range(len(listafiles)):
   try:
    plt.text(spec.Bj_Vj[X],SQM_Vl[X]+0.01,str(spec[X].num),ha='right', va='bottom')
   except:
    print "fallo"
ax1.set_xlabel("Bj - Vj")
ax1.set_ylabel("SQM - V($\lambda$) Photopic vision")
ax1.set_xlim(0,7)
ax1.set_ylim(0.2,1)
#ax1.set_title("SQM-Photopic vs B-V Johnson color")
ax2.set_ylim(0.2-0.35,1-0.35)
ax2.set_ylabel("SQM* - V($\lambda$) Photopic vision")
#plt.xlim(0,0.9)
#plt.ylim(0,2.1)
plt.savefig('colores_espectros_SQM_Vl_BVj.png')
plt.show()
########################################
#Clouds
########################################
#model_ransac_SQM_Vl_BGj,line_y_ransac,line_X=ransacfit(spec.Bj_VjA,SQM_VlA)
#zSQM_VlBGjA = np.polyfit(spec.Bj_VjA, SQM_VlA, 3)
#p = np.poly1d(zSQM_VlBGjA)
#xp = np.linspace(0, 2, 100)
#plt.plot(line_X, line_y_ransac, '-b')#, label='RANSAC regressor')

scatter3(spec.Bj_VjA,SQM_VlA,spec1, 40, alpha)
plt.errorbar(B_VA_t, SQM_VlA_t,ecolor='k',capthick=2, yerr=[SQM_VlA_tmin, SQM_VlA_tmax], xerr=[B_VA_tmin, B_VA_tmax], fmt='.')
scatter3(B_VA_t, SQM_VlA_t,s=500, c=['r','y','g','b','c'], alpha)
#plt.plot(xp, p(xp),'r-')
if 0:
 for X in range(len(listafiles)):
   try:
    plt.text(spec.Bj_VjA[X],SQM_VlA[X]+0.01,str(spec[X].num),ha='right', va='bottom')
   except:
    print "fallo"
plt.xlabel("Bj - Vj")
plt.ylabel("SQM - V($\lambda$) Photopic vision")
plt.xlim(0,7)
plt.ylim(0.2,1)
#plt.title("SQM-Photopic vs B-V Johnson cloud amplified")
#plt.xlim(0,0.9)
#plt.ylim(0,2.1)
plt.savefig('colores_espectros_SQM_Vl_BVjA.png')
plt.show()

###########################################

###########################################
if 0:
    model_ransac_SQM_SC_BGj,line_y_ransac,line_X=ransacfit(spec.Bj_Vj,SQM_SC)
    zSQM_SCBGj = np.polyfit(spec.Bj_Vj, SQM_SC, 3)
    p = np.poly1d(zSQM_SCBGj)
    xp = np.linspace(0, 2, 100)
    #plt.plot(line_X, line_y_ransac, '-b')#, label='RANSAC regressor')
scatter3(spec.Bj_Vj,SQM_SC,spec1, 40, alpha)
plt.errorbar(B_V_t, SQM_SC_t,ecolor='k',capthick=2, yerr=[SQM_SC_tmin, SQM_SC_tmax], xerr=[B_V_tmin, B_V_tmax], fmt='.')
scatter3(B_V_t, SQM_SC_t,s=500, c=['r','y','g','b','c'], alpha)
#plt.plot(xp, p(xp),'r-')
if 0:
 for X in range(len(listafiles)):
   try:
    plt.text(spec.Bj_Vj[X],SQM_SC[X]+0.01,str(spec[X].num),ha='right', va='bottom')
   except:
    print "fallo"
plt.xlabel("Bj - Vj")
plt.ylabel("SQM - V'($\lambda$) Soctopic vision")
plt.xlim(0,7)
plt.ylim(-2,0.4)
#plt.title("SQM-Scotopic vs B-V Johnson color")
#plt.xlim(0,0.9)
#plt.ylim(0,2.1)
plt.savefig('colores_espectros_SQM_SC_BGj.png')
plt.show()

###########################################
# Clouds
###########################################
if 0:
    model_ransac_SQM_SC_BGj,line_y_ransac,line_X=ransacfit(spec.Bj_VjA,SQM_SCA)
    zSQM_SCBGjA = np.polyfit(spec.Bj_VjA, SQM_SCA, 3)
    p = np.poly1d(zSQM_SCBGjA)
    xp = np.linspace(0, 2, 100)
    #plt.plot(line_X, line_y_ransac, '-b')#, label='RANSAC regressor')
scatter3(spec.Bj_VjA,SQM_SCA,spec1, 40, alpha)
plt.errorbar(B_VA_t, SQM_SCA_t,ecolor='k',capthick=2, yerr=[SQM_SCA_tmin, SQM_SCA_tmax], xerr=[B_VA_tmin, B_VA_tmax], fmt='.')
scatter3(B_VA_t, SQM_SCA_t,s=500, c=['r','y','g','b','c'], alpha)
#plt.plot(xp, p(xp),'r-')
if 0:
 for X in range(len(listafiles)):
   try:
    plt.text(spec.Bj_Vj[X],SQM_SC[X]+0.01,str(spec[X].num),ha='right', va='bottom')
   except:
    print "fallo"
plt.xlabel("Bj - Vj")
plt.ylabel("SQM - V'($\lambda$) Soctopic vision")
plt.xlim(0,7)
plt.ylim(-2,0.4)
#plt.title("SQM-Scotopic vs B-V Johnson cloud amplified")
#plt.xlim(0,0.9)
#plt.ylim(0,2.1)
plt.savefig('colores_espectros_SQM_SC_BGjA.png')
plt.show()
###########################################
# Clouds
###########################################
if 0:
    model_ransac_SQM_SQMA_BGj,line_y_ransac,line_X=ransacfit(spec.Bj_VjA,SQM_SQMA)
    zSQM_SQMABGj = np.polyfit(spec.Bj_Vj, SQM_SQMA, 3)
    p = np.poly1d(zSQM_SQMABGj)
    xp = np.linspace(0, 2, 100)
    #plt.plot(line_X, line_y_ransac, '-b')#, label='RANSAC regressor')
scatter3(spec.Bj_Vj,SQM_SQMA,spec1, 40, alpha)
plt.errorbar(B_V_t, SQM_SQMA_t,ecolor='k',capthick=2, yerr=[SQM_SQMA_tmin, SQM_SQMA_tmax], xerr=[B_V_tmin, B_V_tmax], fmt='.')
scatter3(B_V_t, SQM_SQMA_t,s=500, c=['r','y','g','b','c'], alpha)
#plt.plot(xp, p(xp),'r-')
if 0:
 for X in range(len(listafiles)):
   try:
    plt.text(spec.Bj_Vj[X],SQM_SQMA[X]+0.01,str(spec[X].num),ha='right', va='bottom')
   except:
    print "fallo"
plt.xlabel("Bj - Vj")
plt.ylabel("SQM' - SQM")
plt.xlim(0,7)
#plt.ylim(-2,0.4)
#plt.title("SQM - SQM & clouds vs B-V Johnson")
#plt.xlim(0,0.9)
#plt.ylim(0,2.1)
plt.savefig('colores_espectros_SQM_SQMA_BGjA2.png')
plt.show()


############################################################
###################################################3
####################################################
'''
if 1:
    ###########################################

    ###########################################
    '''
    if 0:
        model_ransac_Ha_GR,line_y_ransac,line_X=ransacfit(GR,Ha_G)
        zHa_GGR = np.polyfit(GR, Ha_G, 3)
        p = np.poly1d(zHa_GGR)
        xp = np.linspace(0, 2, 100)
    plt.plot(line_X, line_y_ransac, '-b')#, label='RANSAC regressor')
    scatter3(GR,Ha_G,spec1, 40, alpha)
    plt.plot(xp, p(xp),'r-')
    if 0:
     for X in range(len(listafiles)):
       try:
        plt.text(GR[X],Ha_G[X]+0.01,str(spec[X].num),ha='right', va='bottom')
       except:
        print "fallo"
    plt.xlabel("G/R")
    plt.ylabel("Ha/G")
    plt.xlim(0,1.2)
    plt.ylim(0.,2)
    #plt.title("ISS colors for typical lamps")
    #plt.xlim(0,0.9)
    #plt.ylim(0,2.1)
    plt.savefig('colores_espectros_Ha_G_GR.png')
    plt.show()
    if 0:
        model_ransac_Hb_GR,line_y_ransac,line_X=ransacfit(GR,Hb_G)
        zHb_GGR = np.polyfit(GR, Hb_G, 3)
        p = np.poly1d(zHb_GGR)
        xp = np.linspace(0, 2, 100)
        plt.plot(line_X, line_y_ransac, '-b')#, label='RANSAC regressor')
    scatter3(GR,Hb_G,spec1, 40, alpha)
    plt.plot(xp, p(xp),'r-')
    if 0:
     for X in range(len(listafiles)):
       try:
        plt.text(GR[X],Hb_G[X]+0.01,str(spec[X].num),ha='right', va='bottom')
       except:
        print "fallo"
    plt.xlabel("G/R")
    plt.ylabel("Hb/G")
    plt.xlim(0,1.2)
    plt.ylim(0.,1)
    #plt.title("ISS colors for typical lamps")
    #plt.xlim(0,0.9)
    #plt.ylim(0,2.1)
    plt.savefig('colores_espectros_Hb_G_GR.png')
    plt.show()
    '''
    if 0:
        model_ransac_VIIRSBG,line_y_ransac,line_X=ransacfit(BG,G_VIIRS)
        plt.plot(line_X, line_y_ransac, '-b')#, label='RANSAC regressor')
        zgviirsBG = np.polyfit(BG, G_VIIRS, 2)
        p = np.poly1d(zgviirsBG)
        xp = np.linspace(0, 1, 100)
    scatter3(BG,G_VIIRS,spec1, 40, alpha)
    plt.plot(xp, p(xp),'r-')
    if 0:
     for X in range(len(listafiles)):
       try:
        plt.text(BG[X],G_VIIRS[X]+0.01,str(spec[X].num),ha='right', va='bottom')
       except:
        print "fallo"
    plt.xlabel("B/G")
    plt.ylabel("G/VIIRS")
    #plt.title("VIIRS correction")
    plt.xlim(0,0.9)
    plt.ylim(0,2.4)
    plt.savefig('colores_espectros_VIIRS_BG.png')
    plt.show()

    model_ransac_VIIRSGR,line_y_ransac,line_X=ransacfit(GR,G_VIIRS)
    plt.plot(line_X, line_y_ransac, '-b')#, label='RANSAC regressor')
    zgviirsGR = np.polyfit(GR, G_VIIRS, 2)
    p = np.poly1d(zgviirsGR)
    xp = np.linspace(0, 1, 100)
    scatter3(GR,G_VIIRS,spec1, 40, alpha)
    plt.plot(xp, p(xp),'r-')
    if 0:
     for X in range(len(listafiles)):
       try:
        plt.text(GR[X],G_VIIRS[X]+0.01,str(spec[X].num),ha='right', va='bottom')
       except:
        print "fallo"
    plt.xlabel("G/R")
    plt.ylabel("G/VIIRS")
    #plt.title("VIIRS correction")
    plt.xlim(0,1.8)
    plt.ylim(0,2.4)
    plt.savefig('colores_espectros_VIIRS_GR.png')
    plt.show()


    model_ransac_VlVIIRSGR,line_y_ransac,line_X=ransacfit(GR,VlVIIRS)
    plt.plot(line_X, line_y_ransac, '-b')#, label='RANSAC regressor')
    zgVlVIIRSGR = np.polyfit(GR, VlVIIRS, 2)
    p = np.poly1d(zgVlVIIRSGR)
    xp = np.linspace(0, 2, 100)
    scatter3(GR,VlVIIRS,spec1, 40, alpha)
    plt.plot(xp, p(xp),'r-')
    if 0:
     for X in range(len(listafiles)):
       try:
        plt.text(GR[X],VlVIIRS[X]+0.01,str(spec[X].num),ha='right', va='bottom')
       except:
        print "fallo"
    plt.xlabel("G/R")
    plt.ylabel("Vl/VIIRS")
    #plt.title("VIIRS correction")
    plt.xlim(0,1.8)
    plt.ylim(0,2.4)
    plt.savefig('colores_espectros_VlVIIRSGR.png')
    plt.show()
            
        
    model_ransac_sli,line_y_ransac,line_X=ransacfit(BG,spec.i_sli)
    plt.plot(line_X, line_y_ransac, '-b')#, label='RANSAC regressor')
    zslibg = np.polyfit(BG, spec.i_sli, 2)
    p = np.poly1d(zslibg)
    xp = np.linspace(0, 1, 100)
    scatter3(BG,spec.i_sli,spec1, 40, alpha)
    plt.plot(xp, p(xp),'r-')
    if 0:
     for X in range(len(listafiles)):
       try:
        plt.text(BG[X],spec.i_sli[X]+0.01,str(spec[X].num),ha='right', va='bottom')
       except:
        print "fallo"
    plt.xlabel("B/G")
    plt.ylabel("SLI")
    #plt.title("SLI")
    plt.xlim(0,0.9)
    plt.ylim(0,1.4)
    plt.savefig('colores_espectros_SLI_BG.png')
    plt.show()

    model_ransac_sli,line_y_ransac,line_X=ransacfit(GR,spec.i_sli)
    plt.plot(line_X, line_y_ransac, '-b')#, label='RANSAC regressor')
    zsligr = np.polyfit(GR, spec.i_sli, 3)
    p = np.poly1d(zsligr)
    xp = np.linspace(0, 2, 100)
    scatter3(GR,spec.i_sli,spec1, 40, alpha)
    plt.plot(xp, p(xp),'r-')
    if 0:
     for X in range(len(listafiles)):
       try:
        plt.text(GR[X],spec.i_sli[X]+0.01,str(spec[X].num),ha='right', va='bottom')
       except:
        print "fallo"
    plt.xlabel("G/R")
    plt.ylabel("SLI")
    #plt.title("SLI")
    plt.xlim(0,1.2)
    plt.ylim(0,1.4)
    plt.savefig('colores_espectros_SLI_GR.png')
    plt.show()
    ###################################################3

        
    model_ransac_ipi,line_y_ransac,line_X=ransacfit(BG,spec.i_ipi)
    plt.plot(line_X, line_y_ransac, '-b')#, label='RANSAC regressor')
    scatter3(BG,spec.i_ipi,spec1, 40, alpha)
    if 0:
     for X in range(len(listafiles)):
       try:
        plt.text(BG[X],spec.i_ipi[X]+0.01,str(spec[X].num),ha='right', va='bottom')
       except:
        print "fallo"
    plt.xlabel("B/G")
    plt.ylabel("IPI")
    #plt.title("IPI")
    plt.xlim(0,0.9)
    plt.ylim(0,1.4)
    plt.savefig('colores_espectros_IPI_BG.png')
    plt.show()

    scatter3(GR,spec.i_ipi,spec1, 40, alpha)
    if 0:
     for X in range(len(listafiles)):
       try:
        plt.text(GR[X],spec.i_ipi[X]+0.01,str(spec[X].num),ha='right', va='bottom')
       except:
        print "fallo"
    plt.xlabel("G/R")
    plt.ylabel("IPI")
    #plt.title("IPI")
    plt.xlim(0,1.2)
    plt.ylim(0,1.4)
    plt.savefig('colores_espectros_IPI_GR.png')
    plt.show()

    scatter3(BR,spec.i_ipi,spec1, 40, alpha)
    if 0:
     for X in range(len(listafiles)):
       try:
        plt.text(BR[X],spec.i_ipi[X]+0.01,str(spec[X].num),ha='right', va='bottom')
       except:
        print "fallo"
    plt.xlabel("(B+R)/G")
    plt.ylabel("IPI")
    #plt.title("IPI")
    #plt.xlim(0,1.2)
    #plt.ylim(0,1.4)
    plt.savefig('colores_espectros_IPI_BRG0.png')
    plt.show()

    scatter3(BR1,spec.i_ipi,spec1, 40, alpha)
    if 0:
     for X in range(len(listafiles)):
       try:
        plt.text(BR1[X],spec.i_ipi[X]+0.01,str(spec[X].num),ha='right', va='bottom')
       except:
        print "fallo"
    plt.xlabel("(B-G)/R")
    plt.ylabel("IPI")
    #plt.title("IPI")
    #plt.xlim(0,1.2)
    #plt.ylim(0,1.4)
    plt.savefig('colores_espectros_IPI_BRG1.png')
    plt.show()   
    ###################################################3

        
    model_ransac_msi,line_y_ransac,line_X=ransacfit(BG,spec.i_msi)
    plt.plot(line_X, line_y_ransac, '-b')#, label='RANSAC regressor')
    scatter3(BG,spec.i_msi,spec1, 40, alpha)
    if 0:
     for X in range(len(listafiles)):
       try:
        plt.text(BG[X],spec.i_msi[X]+0.01,str(spec[X].num),ha='right', va='bottom')
       except:
        print "fallo"
    plt.xlabel("B/G")
    plt.ylabel("MSI")
    #plt.title("MSI")
    plt.xlim(0,0.9)
    plt.ylim(0,1.4)
    plt.savefig('colores_espectros_MSI_BG.png')
    plt.show()

    model_ransac_msigr,line_y_ransac,line_X=ransacfit(GR,spec.i_msi)
    plt.plot(line_X, line_y_ransac, '-b')#, label='RANSAC regressor')
    zmsigr = np.polyfit(GR, spec.i_msi, 3)
    p = np.poly1d(zsligr)
    xp = np.linspace(0, 2, 100)
    scatter3(GR,spec.i_msi,spec1, 40, alpha)
    plt.plot(xp, p(xp),'r-')
    if 0:
     for X in range(len(listafiles)):
       try:
        plt.text(GR[X],spec.i_msi[X]+0.01,str(spec[X].num),ha='right', va='bottom')
       except:
        print "fallo"
    plt.xlabel("G/R")
    plt.ylabel("MSI")
    #plt.title("MSI")
    plt.xlim(0,1.2)
    plt.ylim(0,1.4)
    plt.savefig('colores_espectros_MSI_GR.png')
    plt.show()

    model_ransac_jNO3gr,line_y_ransac,line_X=ransacfit(GR,spec.jNO)
    plt.plot(line_X, line_y_ransac, '-b')#, label='RANSAC regressor')
    zmsigr = np.polyfit(GR, spec.jNO, 3)
    p = np.poly1d(zsligr)
    xp = np.linspace(0, 2, 100)
    #scatter3(GR,spec.jNO,spec1, 40, alpha)
    scatter3(GR,spec.jNO,spec1,40,alpha)
    #plt.plot(xp, p(xp),'r-')
    if 0:
     for X in range(len(listafiles)):
       try:
        plt.text(GR[X],spec.jNO[X]+0.0001,str(spec[X].tipo),ha='right', va='bottom')
       except:
        print "fallo"
    plt.xlabel("G/R")
    plt.ylabel("jNO3/lumen")
    #plt.title("MSI")
    plt.xlim(0,1.2)
    #plt.ylim(0,1.4)
    plt.savefig('colores_espectros_jNO3_GR.png')
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
    #plt.title("IPI")
    #plt.xlim(0,1.2)
    #plt.ylim(0,1.4)
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
    

