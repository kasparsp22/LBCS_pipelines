#Note: Load your fits files into AIPS (with the same name as the fits files on your disk) before running this code.
#if the LOTSS sources are named L1234+5678, make sure that the corresponding LBCS sources (on the disk as well as on AIPS)are named as L1234+5678_l (any character after the underscore is okay).

#This code accepts UV datasets loaded in AIPS as inputs and gives complex uv visiblities for each dataset at a particular baseline through the subroutine 'getbaseline'. It then calculates the signal to noise ratio by dividing the mean of the visibilities by the noise (standard deviation). It finally generates a log file containing the name of the sources and the Signal to Noise ratio for that particular baseline. 

from math import *
from AIPS import AIPS, AIPSDisk
from AIPSTask import AIPSTask, AIPSList, AIPSMessageLog
from AIPSData import AIPSUVData, AIPSImage, AIPSCat
from Wizardry.AIPSData import AIPSUVData as WizAIPSUVData
from scipy import ndimage; from scipy.ndimage import measurements
import matplotlib; from matplotlib import pyplot as plt
import pyfits; from pyfits import getdata,getheader
import re,sys,pickle,numpy as np,os,glob,time,warnings; from numpy import fft
import cmath

aipsno = 3
baseline_DE = ['ST001','DE605']
baseline_UK = ['ST001','UK608']
incl = 'FITS'

#subroutine to give the visibilites of the uv datasets loaded into aips
#if this one doesnt work, use the getbaseline subroutine given in snr.py

def getbaseline (aipsno,baseline,inna,incl,indisk=1,inseq=1,minmatch=5,ipol=0):
    AIPS.userno = aipsno
    print 'AIPS no is', aipsno
    h = AIPSUVData (inna,incl,indisk,inseq).header
    print 'got AIPSUVData'
    f0, deltaf = h['crval'][2], h['cdelt'][2]
    data = WizAIPSUVData (inna, incl, indisk, inseq)
    print "****** FITS file ->"
    print data
    antnum = np.array([],dtype='int')
    foundit = False
    for i in range (len(data.antennas)):
        if data.antennas[i][:minmatch] in baseline:
            antnum = np.append(antnum,i)
#            print data.antennas[i],i
    antnum = np.sort(antnum)
#    print antnum
    atable = data.table('AN',1)
    newantnum = [atable[antnum[0]]['nosta'],atable[antnum[1]]['nosta']]
#    print newantnum
    for v in data:
#	print 'for v in data'
#	print v.baseline, newantnum
        
        if v.baseline[0]!=newantnum[0] or v.baseline[1]!=newantnum[1]:
            continue
#	print 'here now'
        vu = v.visibility
	#print 'vu is', vu
        chans = np.arange(vu.shape[0]*vu.shape[1])
        thisuvw = np.outer(1.+chans*deltaf/f0,v.uvw)
        vv = vu.reshape(vu.shape[0]*vu.shape[1],vu.shape[2],vu.shape[3])
        vthis = vv[:,ipol,0] + 1j*vv[:,ipol,1]
#	print 'reached here!'
        try:

#	    print 'try'
            d = np.vstack((d,vthis))
            uvw = np.dstack((uvw,thisuvw))
	    
        except:
#	    print 'except'
            d = np.copy (vthis)
            uvw = np.copy(thisuvw)
	    
    uvw = np.rollaxis(uvw,2,0)
    return d,uvw

print 'The UV datasets are:'
print '\n'

#the line below assumes that the names of the datasets loaded into AIPS are the same as the names of their corresponding fits files on the disk

for i in os.listdir('/data/data-disk1/data/LBCS/scripts_k/LOTSS_sources'): #put location of the fits files on the disk
#	print i.split('.')[0]


	vis1, uvw1 = getbaseline(aipsno,baseline_UK,i.split('.')[0],incl)
	print i.split('.')[0]
	r, phi = cmath.polar(np.mean(vis1)) 
	print 'abs value of mean is', r
	print 'mean/std dev', i.split('.')[0], 'is', cmath.polar(np.mean(vis1)/np.std(vis1))[0]
	r_sc, phi_sc = cmath.polar(np.std(vis1-np.mean(vis1))) #just doing np.std(vis1) would give the same thing... 
	print r_sc
	print 'SNR is', r/r_sc #signal to noise
	print '###'

	file_1 = open('log_lbcs_UK_1.log','a') #log file to dump the sources names and their SNRs
	file_1.write(str(i.split('.')[0]+ '    ' + str(r/r_sc)+ '\n'))
	file_1.close()
	

print 'Done for both datasets'
