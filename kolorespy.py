#!/usr/bin/env python

# IPython log file

import asciitable,sys
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad
from scipy.interpolate import UnivariateSpline
from scipy.interpolate import InterpolatedUnivariateSpline
def normaliza(X):
    max=np.max(X.intrel)
    X.intrel=X.intrel/max
    return X

waveref=np.arange(250,15000,1)
NPMAX=10000   
WL_VEGA_AB=[]
FLUX_VEGA_AB=[]
for X in range(NPMAX):
    WL_VEGA1=1000.0+(float(X)/(NPMAX-1)*(30000.0-1000.0))
    WL_VEGA_AB.append(WL_VEGA1)
    FLUX_VEGA1=0.1088/(WL_VEGA1*WL_VEGA1)
    FLUX_VEGA_AB.append(FLUX_VEGA1)

WL_VEGA_ST=[]
FLUX_VEGA_ST=[]
for X in range(NPMAX):
    WL_VEGA1=1000.0+(float(X)/(NPMAX-1)*(30000.0-1000.0))
    WL_VEGA_ST.append(WL_VEGA1)
    FLUX_VEGA1=3.63E-9
    FLUX_VEGA_ST.append(FLUX_VEGA1)

#FLUX_VEGA_AB=np.array(FLUX_VEGA_AB,dtype=[('flux', float)])
#WL_VEGA_AB=np.array(WL_VEGA_AB,dtype=[('wave', float)])

#AB=zip(WL_VEGA_AB,FLUX_VEGA_AB)
AB=np.recarray((len(FLUX_VEGA_AB),),dtype=[('wave', float), ('flux', float)])
AB.wave=np.array(WL_VEGA_AB)
AB.flux=np.array(FLUX_VEGA_AB) #erg s-1 cm-2 Amstrong-1
sun_lambda=asciitable.read('SUN_CIE.txt')

#filename='outfile.dat'
filename=sys.argv[1]
spectro=asciitable.read(filename)
B_D3=asciitable.read("B_D3.csv")
G_D3=asciitable.read("G_D3.csv")
R_D3=asciitable.read("R_D3.csv")
B_D3s=asciitable.read("B_D3s.csv")
G_D3s=asciitable.read("G_D3s.csv")
R_D3s=asciitable.read("R_D3s.csv")
U_johnson=asciitable.read("U_johnson.csv")
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

def normaliza(X):
    max=np.max(X.intrel)
    X.intrel=X.intrel/max
    return X

fig, (ax0,ax1,ax2) = plt.subplots(nrows=3)

ax0.plot(WL_VEGA_AB, FLUX_VEGA_AB)
ax0.plot(WL_VEGA_ST, FLUX_VEGA_ST)
ax0.set_title('Reference AB')
ax1.plot(spectro.col2, spectro.col1)
ax1.set_title('Spectra')
ax2.plot(photo.wave,photo.intrel)
ax2.plot(psas.wave,psas.intrel)
ax2.plot(scoto.wave,scoto.intrel)
ax2.plot(msas.wave,msas.intrel)
ax2.set_title('Filter')
plt.show()
#AB_ref = np.array(WL_VEGA_AB, dtype=[('wave', float)])
def mag(filtroX,espectro,referencia):
    if 0:
     if filtroX.wave[-1]>espectro.col2[-1]:
        dif=filtroX.wave[-1]-espectro.col2[-1]
        n=espectro.size
        new=n+int(dif)+2
        espectro=np.resize(espectro,new)
        espectro[n:-1]['col2']=np.linspace(espectro['col2'][n-1],filtroX.wave[-1],num=espectro[n:-1]['col2'].size)

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
    I=Up[0]/Dn[0]
    return m,I

def IndiceIsoLux(filtroX,espectro,referencia,photo):
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
Uj_Vj=m_Uj-m_Vj
Bj_Vj=m_Bj-m_Vj
Rj_Vj=m_Rj-m_Vj
Bn_Vj=m_Bn-m_Vj
Gn_Vj=m_Gn-m_Vj
Rn_Vj=m_Rn-m_Vj
SQM_Vj=m_SQM-m_Vj
Bn_Gn=I_Bn/I_Gn
Gn_Rn=I_Gn/I_Rn

c=299792458.0
frec=540*10**12
lamnda=c/frec
Inten= 1.0/683#[W.str]
waveref=np.arange(250,15000,1)
intesity=np.arange(250.0,15000.0,1)*0
ncd=list(waveref).index(int(lamnda*1E10))
intesity[ncd]=1.0/683

intensityAB=Wtoergs(X)

cd=np.recarray((len(FLUX_VEGA_AB),),dtype=[('wave', float), ('flux', float)])
cd.wave=np.array(waveref)
cd.flux=np.array(intesity)

def srtoarcsec2(X):
    X=X* 4.254517e+10
    return X

def srtoarcsec2(X):
    X=X* 2.3504431e-11
    return X

def ergstoW(X):
      X=X* 1e-07
      return X
def Wtoergs(X):
      X=X/1e-07
      return X

def cm_2tom_2(X):
      X=X*10000
      return X

def m_2tocm_2(X):
      X=X/10000
      return X

def Angstrontonm(X):
      X=X/10
      return X

def nmtoAngstron(X):
      X=X*10
      return X


def ABtoFlux(AB,filtro):
    w=np.sum(filtro.intrel)/float(filtro.size)
    Flux=10^((AB-5*np.log10(w)-2.401)/-2.5)
    return Flux
def FluxtoAB(Flux,filtro):
    w=np.sum(filtro.intrel)/float(filtro.size)
    AB = -2.5*np.log10(Flux)-5*np.log10(w)+2.401 #[W cm -2 arsec-2 amstrong-1]
    return AB

#def FluxtoCandela(Flux,filtro,photo):
    

i_sli=IndiceIsoLux(scoto,spectro,sun_lambda,photo)
i_msi=IndiceIsoLux(msas,spectro,sun_lambda,photo)
i_ipi=IndiceIsoLux(psas,spectro,sun_lambda,photo)

