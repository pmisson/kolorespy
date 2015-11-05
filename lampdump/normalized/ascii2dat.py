#!/usr/bin/env python

# IPython log file
import asciitable,sys
import numpy as np
#filename='LED_cree.ascii'
filename=sys.argv[1]
datos=asciitable.read(filename)
datos.wavelength=datos.wavelength*10
if datos.wavelength[-1]<10396:
    size1=datos.size
    nextsize=np.int((10396.-datos.wavelength[-1])/50.0)
    datoswavelength=np.zeros(size1+nextsize)
    datosflux=np.zeros(size1+nextsize)
    datoswavelength[0:size1]=datos.wavelength
    datosflux[0:size1]=datos.intensity
    newflux=np.zeros(nextsize,)
    newwave=np.linspace(datos.wavelength[-1]+50,10396.,nextsize)
    datoswavelength[size1:]=newwave
    datosflux[size1:]=newflux
    asciitable.write({'intensity': datosflux, 'wavelength': datoswavelength}, filename.split('.')[0]+'.dat',Writer=asciitable.NoHeader, names=['intensity', 'wavelength'])
    asciitable.write({ 'wavelength': datoswavelength,'intensity': datosflux}, filename.split('.')[0]+'.txt',Writer=asciitable.NoHeader, names=[ 'wavelength','intensity'])
else:    
    asciitable.write({'intensity': datos.intensity, 'wavelength': datos.wavelength}, filename.split('.')[0]+'.dat',Writer=asciitable.NoHeader, names=['intensity', 'wavelength'])
    asciitable.write({ 'wavelength': datos.wavelength,'intensity': datos.intensity}, filename.split('.')[0]+'.txt',Writer=asciitable.NoHeader, names=[ 'wavelength','intensity'])

