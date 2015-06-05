#!/usr/bin/env python

# IPython log file
import asciitable,sys
filename=sys.argv[1]
datos=asciitable.read(filename)
datos.wavelength=datos.wavelength*10
asciitable.write({'relativeIntensity': datos.relativeIntensity, 'wavelength': datos.wavelength}, filename.split('.')[0]+'.dat',Writer=asciitable.NoHeader, names=['relativeIntensity', 'wavelength'])
asciitable.write({ 'wavelength': datos.wavelength,'relativeIntensity': datos.relativeIntensity}, filename.split('.')[0]+'.txt',Writer=asciitable.NoHeader, names=[ 'wavelength','relativeIntensity'])

