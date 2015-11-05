#!/usr/bin/env python

#By A. Sanchez de Miguel v1.2
import asciitable
import os,sys
from scipy.interpolate import UnivariateSpline
from scipy.interpolate import InterpolatedUnivariateSpline
namefile=sys.argv[1]

def gen(name1,STWV,disper):

	fichero="genSlit_B.com"

	f=open(fichero,"w")

	print >> f,"genimage<<teescribo" 
	print >> f,"n" 		#Initialize frame to existing file (y/n) [n] ?
	print >> f,"1"	#NSCAN  (1,...,2700)?
	print >> f,"2048"	#NSCAN  (1,...,2700)?
	print >> f,"5"		#Option?
	print >> f,"1,2048"	#Scan region (0,0=EXIT)...?
	print >> f,"2"	#2) Table X,Y: X=wavelength, Y=flux, linear interpolation
	print >> f,"outfile.dat"	#Spectrum file name?
	print >> f,STWV	#STWV     [0.0000000] ?
	print >> f,disper	#DISP     [0.0000000] ?
	print >> f,"1"		#Column number where wavelength is located  (1,...,100) [1] ?#WARNING: STWV
	print >> f,"2"		#Column number where flux is located......  (1,...,100) [2] ?
	print >> f,"1"      #Factor to be applied to the wavelength values [1.0] ?
	print >> f,"0"      #Radial velocity (km/sec) [0.0] ?
	print >> f,"1"      #Factor to be applied to the flux values [1.0] ?
	print >> f,"0"		#Option?
	print >> f,"n"		#Change head information (y/n) [n] ?
	print >> f,name1	#Output file name?
	print >> f,"teescribo"
	
	f.close()
	comando2="chmod +x genSlit_B.com"
	comando1="./genSlit_B.com"

	os.system(comando2)
	os.system(comando1)

def writefits(name1):

	fichero="write_B.com"

	f=open(fichero,"w")

	print >> f,"writefits<<teescribo" 
	print >> f,name1	#Spectrum file name?
	print >> f,name1[:-2]+'.fits'	#Output file name?
	print >> f,"n"      #Are you including the FITS header from another file (y/n) [n] ?
	print >> f,"teescribo"



	f.close()
	comando2="chmod +x write_B.com"
	comando1="./write_B.com"

	os.system(comando2)
	os.system(comando1)
	
espc12=asciitable.read(namefile,header_start=17,data_start=18,data_end=2066,guess=False,delimiter='\t')
curvaresp=asciitable.read('CorrectionJAZ.cor')
x=espc12.W*10
STWV=x[0]
disper=(x[-1]-x[0])/2048
y=espc12.S#-espc12.D
curvaresp.col2[curvaresp.col1<300]=0
y=y*curvaresp.col2
asciitable.write({'x': x, 'y': y}, 'outfile.dat',Writer=asciitable.NoHeader )
asciitable.write({'x': x, 'y': y}, namefile[:-4]+'.dat',Writer=asciitable.NoHeader,names=['y','x'] )
name1=namefile[:-4]+'.u'
gen(name1,STWV,disper)
writefits(name1)

exit()
