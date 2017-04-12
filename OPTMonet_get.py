#! /usr/bin/env python

# OPTMonet_get.py
#
# - J.M.C. Court, 2017
#
# Call as 'python OPTMonet_get.py'
#
# Produces a lightcurve from a number of given Monet files.
#
# Requires 3 inputs:
#
#  Path to the directory containing all Monet fits files
#  Path to the file containing list of desired Monet fits files
#  Radius of extraction region in arcseconds
#
# File containing list of desired files should be formatted as such:
#  One filename per line, with extension
#  Full path from the directory given as input 1
#

from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
import pylab as pl
import time as ti
import calendar as ca

indire=raw_input('Full path to directory of files  : ')
print('')
infile=raw_input('File containing infile names     : ')
print('')

try:
   exreg=float(raw_input('Radius of extraction region (sec): '))/3600.0
except:
   print('Must give number!')
   exit()

try:
   obsf=open(infile)
except:
   print('Cannot find',infile+'!')
   print('Aborting!')
   exit()

obss=[]

for line in obsf:
   obss.append(indire+'/'+line[:-1])

rates={}
times={}

for obs in obss:
   cts=0
   try:
      f=fits.open(obs)[0]
   except:
      print('Warning!  '+obs+' does not appear to be Monet fits file!')
      continue
   w=WCS(obs)
   objra=float(f.header['OBJRA'])
   obdec=float(f.header['OBJDEC'])
   filtr=f.header['FILTER']
   if filtr not in rates.keys():
      rates[filtr]=[]
      times[filtr]=[]
   xdelt=abs(float(f.header['CDELT1']))
   ydelt=abs(float(f.header['CDELT2']))
   timedate=f.header['DATE-OBS'].split('.')
   timel=ti.strptime(timedate[0], "%Y-%m-%dT%H:%M:%S")
   timeins=ca.timegm(timel)+float('0.'+timedate[1])
   table=f.data
   objcoords=w.all_world2pix([(objra,obdec)],1)[0]
   objx=objcoords[1]
   objy=objcoords[0]

   #---Diagnostic Block 1---

   #mn=np.mean(table)
   #ptable=f.data
   #ptable[ptable>mn+100]=mn+100
   #ptable[ptable<mn]=mn
   #pl.figure()
   #pl.imshow(ptable)
   #pl.ylim(0,1000)
   #pl.axhline(objx)
   #pl.axvline(objy)
   #pl.show(block=False)

   mini=objx-(exreg/xdelt)
   minj=objy-(exreg/ydelt)
   maxi=objx+(exreg/xdelt)+1
   maxj=objy+(exreg/ydelt)+1
   table=table[int(mini):int(maxi),int(minj):int(maxj)]
   for i in np.arange(len(table),dtype=int):
      for j in np.arange(len(table[0]),dtype=int):
         if np.sqrt(((i+mini-objx)*xdelt)**2+((j+minj-objy)*ydelt)**2)<exreg:
            cts+=table[i,j]

   #---Diagnostic Block 2---

   #      else:
   #         table[i,j]=mn

   times[filtr].append(timeins)
   rates[filtr].append(cts)

   #---Diagnostic Block 3---

   #pl.figure()
   #pl.imshow(table)
   #pl.show(block=True)

for filtr in times.keys():
   times[filtr]=(np.array(times[filtr])+3506716799)/(24*3600)

pl.figure()
for filtr in times.keys():
   pl.plot(times[filtr],rates[filtr],label=filtr)
pl.legend()
pl.xlabel('MJD')
pl.ylabel('Counts')
pl.show(block=True)

