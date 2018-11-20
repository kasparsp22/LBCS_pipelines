from math import *
from AIPS import AIPS, AIPSDisk
from AIPSTask import AIPSTask, AIPSList, AIPSMessageLog
from AIPSData import AIPSUVData, AIPSImage, AIPSCat
from Wizardry.AIPSData import AIPSUVData as WizAIPSUVData
from scipy import ndimage; from scipy.ndimage import measurements
import matplotlib; from matplotlib import pyplot as plt
import pyfits; from pyfits import getdata,getheader
import re,sys,pickle,numpy as np,os,glob,time,warnings; from numpy import fft

def pmulti (aipsname, indisk, aipsclass='FITS'):
    multi = AIPSTask('MULTI')
    multi.indata = AIPSUVData (aipsname, aipsclass, indisk, 1)
    multi.outdata = AIPSUVData (aipsname, 'FITS', indisk, 2)
    multi.go()

def ptacop (name1, class1, name2, class2, ttype, innum, outnum, disk):
    tacop = AIPSTask('TACOP')
    tacop.indata = AIPSUVData(name1, class1, disk, 1)
    tacop.outdata = AIPSUVData(name2, class2, disk, 1)
    tacop.inext = ttype
    tacop.invers = innum
    tacop.outvers = outnum
    tacop.go()
    
def pwtmod (aipsname, refant, antennas, supweight=50.0, indisk=1):
    uvdata = AIPSUVData (aipsname, 'FITS', indisk, 1)
    wtmod = AIPSTask ('WTMOD')
    wtmod.indata = uvdata
    wtmod.outdata = uvdata
    wtmod.antwt[1:] = [1]*len(antennas)
    wtmod.antwt[antennas.index(refant)+1] = supweight
#    wtmod.inp()
    wtmod.go()

# dparm(8): binary with 1=rates, 2=delays, 4=phase

def pfring (aipsname,refant,antennas,source,indisk=1,delaywin=600,ratewin=20,\
            solint=1,snr=2,dparm4=0,logdir='./',weightit=1,zero=0,aipsclass='FITS',\
            aipsseq=1,suppress_rate=0):
    uvdata = AIPSUVData (aipsname,aipsclass,indisk,aipsseq)
    fring = AIPSTask ('FRING')
    fring.refant = refant
    fring.indata = uvdata
    fring.calsour[1:] = [source]
    fring.antennas[1:] = antennas
    fring.solint = solint
    fring.aparm[1:] = [0,0,0,0,0,2,snr,0,0]
    fring.dparm[1:] = [0,delaywin,ratewin,dparm4,0,0,0,zero,suppress_rate]
    fring.weightit = weightit
    stdout = sys.stdout; sys.stdout = open(logdir+aipsname+'.log','a')
#    fring.inp()
    fring.go()
    sys.stdout.close(); sys.stdout = stdout

def pclcal (aipsname,indisk,inver,aipsclass='FITS',logdir='./',\
            snver=-1,gainver=0,gainuse=0, aipsseq=1):
    uvdata = AIPSUVData (aipsname,aipsclass,indisk,aipsseq)
    clcal = AIPSTask ('clcal')
    clcal.indata = uvdata
    clcal.inver = inver
    clcal.snver = inver if snver==-1 else snver
    clcal.gainver = gainver
    clcal.gainuse = gainuse
    stdout = sys.stdout; sys.stdout = open(logdir+aipsname+'.log','a')
    clcal.go()
    sys.stdout.close(); sys.stdout = stdout

def psplit (aipsname,source,indisk,logdir='./'):
    uvdata = AIPSUVData (aipsname,'FITS',indisk,1)
    split = AIPSTask ('split')
    split.indata = uvdata
    split.outclass = 'SPLIT'
    split.docalib = 1
    stdout = sys.stdout; sys.stdout = open(logdir+aipsname+'.log','a')
    split.go()
    sys.stdout.close(); sys.stdout = stdout
    uvdata = AIPSUVData(source,'SPLIT',indisk,1)
    uvdata.rename(aipsname,'SPLIT',1)

def pload (filename,aipsname,indisk,outcl,logdir='./',doindxr=True):
    fitld = AIPSTask ('FITLD')
    fitld.datain = str(filename)
    fitld.outna = aipsname
    fitld.outcl = outcl
    fitld.dokeep = 1
    fitld.outdisk = indisk
#    fitld.inp()
    stdout = sys.stdout; sys.stdout = open(logdir+aipsname+'.log','a')
    fitld.go()
    if doindxr:
        uvdata = AIPSUVData (aipsname,'FITS',1,1)
        indxr = AIPSTask ('INDXR')
        indxr.cparm[1:] = [0,0,0.1,0,0,0,0,0,0,0]
        indxr.indata = uvdata
        indxr.go()
    sys.stdout.close(); sys.stdout = stdout

def stars (aipsname, incl, indisk, intext='./starsfile',logdir='./'):
    stars = AIPSTask('stars')
    stars.inname = aipsname
    stars.inclass = incl
    stars.indisk = indisk
    try:
        stars.stvers = 0    # does not exist in some AIPS versions
    except:
        pass
    stars.intext = './starsfile'
    stdout = sys.stdout; sys.stdout = open(logdir+aipsname+'.log','a')
#    stars.inp()
    try:
        stars.go()
    except:
        pass
    sys.stdout.close(); sys.stdout = stdout

def greys (aipsname, incl, indisk, pmin, pmax, stfac, stvers, logdir='./'):
    greys = AIPSTask('greys')
    greys.inname = aipsname
    greys.inclass = incl
    greys.indisk = indisk
    greys.pixrange[1:] = [float(pmin),float(pmax)]
    greys.dotv = -1
    greys.stfac = stfac
    try:
        greys.stvers = stvers  # does not exist in some aips versions
    except:
        pass
    stdout = sys.stdout; sys.stdout = open(logdir+aipsname+'.log','a')
    try:
        greys.go()
    except:
        pass
    sys.stdout.close(); sys.stdout = stdout

def lwpla (aipsname,incl,indisk,outfile,logdir='./'):
    lwpla = AIPSTask('lwpla')
    lwpla.inname = aipsname
    lwpla.inclass = incl
    lwpla.indisk = indisk
    lwpla.outfile = outfile
    stdout = sys.stdout; sys.stdout = open(logdir+aipsname+'.log','a')
    try:
        lwpla.go()
    except:
        pass
    sys.stdout.close(); sys.stdout = stdout

def pimagr (aipsname,aipsclass,docalib,source='',imsize=256,cellsize=0.1,\
            nchav=-1,disk=1,gainuse=0,robust=0,niter=100,outname='',\
            antennas = [], stokes='I'):
    imagr = AIPSTask('imagr')
    imagr.inname = aipsname
    imagr.inclass = aipsclass
    imagr.docalib = docalib
    imagr.indisk = disk
    imagr.stokes = stokes
    imagr.imsize[1:] = [imsize,imsize]
    imagr.cellsize[1:] = [cellsize,cellsize]
    if len(antennas):
        imagr.antennas[1:] = antennas
#        imagr.baseline[1:] = antennas
    if source=='':
        try:
            su = AIPSUVData(aipsname,aipsclass,disk,1).table('SU',1)
            imagr.source[1] = su[0]['source']
        except:
            print 'Unable to extract source from source table'
    if nchav==-1:
        try:
            h = AIPSUVData(aipsname,aipsclass,disk,1).header
            imagr.nchav = h['naxis'][h['ctype'].index('FREQ')]
        except:
            imagr.nchav=0
    imagr.gainuse = gainuse
    imagr.robust = robust
    imagr.niter = niter
    imagr.outname = outname
#    imagr.inp()
    imagr.go()
