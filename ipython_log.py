# IPython log file

import asciitable
filename='LED_cree.dat'
datos=asciitable.read(filename)
datos.wavelength=datos.wavelength*10
asciitable.write({'relativeIntensity': datos.relativeIntensity, 'wavelength': datos.wavelength}, filename.split('.')[0]+'.dat',Writer=asciitable.NoHeader, names=['relativeIntensity', 'wavelength'])

