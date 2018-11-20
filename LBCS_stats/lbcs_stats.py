#
#   Processes a directory full of .pic files into a LBCS calibrator list
#   See Jackson et al. 2016 for full details
#   Normally the routines getstats() and plotstats() will need to be run
#
import numpy as np, os, sys, glob, matplotlib, scipy; from scipy import stats
from matplotlib import pyplot as plt, path as mpath
from correlate import *
CORE2EFF = 266000.0
tels = ['DE601','DE602','DE603','DE604','DE605','FR606','SE607','UK608','DE609','PL610','PL611','PL612','IE613']
telc = ['r','y','m','g','k','b','c','#008000','#444444','#99EEAA','#CCBB66','#0000E0','#55AA22']
nstats = 8
# 2-d numpy array of vertices, then bbPath = mpath.Path(array), then bbPath.contains_point((x,y))
ppath = np.array([[-0.01,-0.01],[-0.01,20],[50,45],[100,8],[50,-0.01]])
spath = np.array([[100,8],[50,45],[85,60],[150,16]])
tpath = np.array([[-0.01,20],[-0.01,130],[85,130],[85,60]])
upath = np.array([[50,-0.01],[150,16],[350,25],[350,-0.01]])
pPath,sPath,tPath,uPath = mpath.Path(ppath),mpath.Path(spath),mpath.Path(tpath),mpath.Path(upath)
tlabel = ['Delay scatter/ns','Delay scatter L-R/ns','Delay L-R/ns',\
           'Rate scatter/mHz','Rate L-R/mHz','SNR','SNR-4min','Phase scatter L-R/deg']
tlim = [[0,350],[0,350],[0,350],[0,100],[0,100],[0,100],[0,1000],[0,130]]
titlestring = '       DELAY  DSCAT_LR  D_DELAY_4   RATE   RATE_4  SNR  SNR_4 PHSCAT_LR  UV  N\n'

def hms2decimal (c, sep):
    cs = np.asarray(c.split(sep),dtype='f')
    return cs[0]*15.0 + cs[1]/4.0 + cs[2]/240.0

def dms2decimal (c, sep):
    cs = c.split(sep)
    issouth = -1.0 if '-' in cs[0] else 1.0
    cs = np.asarray(cs,dtype='f')
    return issouth*(abs(cs[0])+cs[1]/60.0+cs[2]/3600.0)

def lbcsjvas (jvascat='/home/njj/catalogues/jvas.txt'):
    lbcs = np.loadtxt('lbcs_stats.sum',dtype='S')
    jvas = np.loadtxt(jvascat, dtype='S')
    for i in lbcs:
        new = np.array([hms2decimal(i[1],':'),dms2decimal(i[2],':')])
        try:
            clbcs = np.vstack((clbcs,new))
        except:
            clbcs = np.copy(new)
    for i in jvas:    #2-4,5-7
        new = np.array([hms2decimal(i[2]+':'+i[3]+':'+i[4],':'),\
                        dms2decimal(i[5]+':'+i[6]+':'+i[7],':')])
        try:
            cjvas = np.vstack((cjvas,new))
        except:
            cjvas = np.copy(new)
    a=correlate(clbcs,0,1,cjvas,0,1,10./3600.)
    f=open('lbcs_in_jvas','w')
    for i in a:
        f.write('%s\n'%(lbcs[i[0],0]))
    f.close()

def transfer (picdir):
    os.system('grep P lbcs_stats.sum >transfer_temp')
    dat = np.loadtxt('transfer_temp',dtype='S')
    for i in dat:
        new = np.array([hms2decimal(i[1],':'),dms2decimal(i[2],':')])
        try:
            coords = np.vstack((coords,new))
        except:
            coords = np.copy(new)

    a = correlate (coords,0,1,coords,0,1,2.0)
    a = a[a[:,0]<a[:,1]]
    # Select correlations between observations at the same time

    for i in a:
        i0,i1 = int(i[0]), int(i[1])
        condtime = dat[i0,3]==dat[i1,3] and dat[i0,4]==dat[i1,4]
        condP = dat[i0,5].count('P')>5 and dat[i1,5].count('P')>5
        if condtime:
            try:
                atime = np.vstack((atime,i))
            except:
                atime = np.copy(i)
            if condP:
                try:
                    aP = np.vstack((aP,i))
                except:
                    aP = np.copy(i)

    fo = open('transfer_sources','w')
    for i in aP:
        i0,i1 = int(i[0]),int(i[1])
        pcount = dat[i0,5].count('P')+dat[i1,5].count('P')
        scount = dat[i0,5].count('S')+dat[i1,5].count('S')
        xcount = dat[i0,5].count('X')+dat[i1,5].count('X')
        prop = float(pcount)/float(pcount+scount+xcount)
        fo.write('%s %s %.3f %.3f\n'%(dat[i0,0],dat[i1,0],i[2],prop))
    fo.close()


def coherence (picdir,solint=0.1,dtype='phase',dmax=500,outdir='',elevlim=0.0):
    coh = [[] for i in range(len(tels))]
    picdir = picdir+'/' if not picdir[-1]=='/' else picdir
    a = np.sort(glob.glob(picdir+'*'))
    dat = np.loadtxt('lbcs_stats.sum',dtype='S')
    fo = open('coherence.log','w')
    for i in dat:
        try:
            pidx = np.argwhere(a==picdir+i[0]+'.pic')[0][0]
        except:
            continue
        if 'S' in i[5] or 'X' in i[5] or 'D' in i[5]:  # only very good sources
            continue
        pic = np.load(a[pidx])
        y = []; py = []
        fo.write(i[0]+' ')
        for j in range(len(tels)):
#            if elev(pic,tels[j])<elevlim or elev(pic,'CS001')<elevlim:
#                fo.write('   nan    nan ')
#                continue
            if tels[j] in pic['antennas']:
                tidx=pic['antennas'].index(tels[j])
            else:
                y.append(np.array([]));py.append(np.array([]))
                fo.write('%6f %6f '%(np.nan,np.nan))
                continue
            pl,pr = pic[dtype][tidx,:,0],pic[dtype][tidx,:,1]
            pl,pr = pl[~np.isnan(pl)],pr[~np.isnan(pr)]
            if dtype=='phase':
                pl,pr = np.unwrap(pl),np.unwrap(pr)
            el = np.polyfit(np.arange(len(pl)),pl,2)
            er = np.polyfit(np.arange(len(pr)),pr,2)
            fl = np.poly1d(el)(np.arange(len(pl)))
            fr = np.poly1d(er)(np.arange(len(pr)))
            y.append (fl if np.std(pl-fl) < np.std(pr-fr) else fr)
            py.append (pl if np.std(pl-fl) < np.std(pr-fr) else pr)
            if dtype == 'phase':    # in s (for phase change 1 rad)
                val = solint*60./np.mean(np.abs(np.gradient(y[-1])))
            else:
                val = np.mean(y[-1])
            coh[j].append(val)
            vstd = 1000.*np.std(py[-1]-y[-1])
            fo.write('%6.2f '%(val) if val<1000. else '%6.1f '%(val))
            fo.write('%6.2f '%(vstd) if vstd<1000. else '%6.1f '%(vstd))
        fo.write('\n')
        if len(outdir):
            maxrange = 0.0
            for j in range(len(tels)):
                if len(py[j]):
                    maxrange = max(maxrange,1.1*(py[j].max()-py[j].min()))
            for j in range(len(tels)):
                if len(py[j]):
                    plt.subplot(521+j,xticks=[])
                    scat = 1000.*np.std(py[j]-y[j])
                    cohv = coh[j][-1]
                    scat = 'nan' if (np.isnan(scat) or np.isinf(scat)) else str(int(scat))
                    cohv = 'nan' if (np.isnan(cohv) or np.isinf(cohv)) else str(int(cohv))
                    label=tels[j]+' '+cohv+' '+scat
                    yhalf = 0.5*(py[j].max()+py[j].min())
                    plt.plot(py[j],'bo')
                    plt.plot(y[j],'r-',label=label)
                    plt.ylim (yhalf-0.5*maxrange, yhalf+0.5*maxrange)
                    plt.legend(handlelength=0.0001,loc=0,fontsize=7)
            plt.suptitle(i[0])
            outdir = outdir+'/' if outdir[-1]!='/' else outdir
            if not os.path.isdir(outdir):
                os.system('mkdir '+outdir)
            plt.savefig(outdir+'C'+i[0]+'.png')
            plt.clf()
    fo.close()

def cohstats (stdlim=300,doprob=False,outfile='lbcs_coh.png',dmax=500):
    a = np.asarray(np.loadtxt('coherence.log',dtype='S')[:,1:],dtype='float')
    ks = np.zeros((len(tels),len(tels)))
    sys.stdout.write('      ')
    for tel in tels:
        sys.stdout.write(tel+' ')
    sys.stdout.write('\n')
    lengths = np.loadtxt('IBlengths')[1:,1:]
    
    boxes = [7,8,2,5,1,9,3,6,4]
    for i in range(len(tels)):
        samp1,std1 = a[:,2*i],a[:,2*i+1]
        std1 = std1[(~np.isnan(samp1))&(~np.isinf(samp1))]
        samp1 = samp1[(~np.isnan(samp1))&(~np.isinf(samp1))]
        samp1 = samp1[std1<stdlim]
        plt.subplot(330+boxes[i],yticks=[])
#        plt.subplot(331+i,yticks=[])
        plt.hist(samp1,range=[0,dmax],bins=25,label=tels[i])
        print tels[i],np.median(samp1)
        plt.legend(handlelength=0.0)
        for j in range(len(tels)):
            samp2,std2 = a[:,2*j],a[:,2*j+1]
            std2 = std2[(~np.isnan(samp2))&(~np.isinf(samp2))]
            samp2 = samp2[(~np.isnan(samp2))&(~np.isinf(samp2))]
            samp2 = samp2[std2<stdlim]
            ks[i,j] = (-np.log10(stats.ks_2samp(samp1,samp2)[1])) if doprob \
                      else stats.ks_2samp(samp1,samp2)[0]
    plt.savefig(outfile,bbox_inches='tight'); plt.clf()
    px = []; py = []
    for i in range(9):
        sys.stdout.write(tels[i]+' ')
        for j in range(9):
            sys.stdout.write('%5.1f '%ks[i,j] if doprob else '%5.2f '%ks[i,j])
            if j>i:
                px.append(lengths[i,j])
                py.append(ks[i,j])
        sys.stdout.write('\n')
    plt.subplot(111)
    plt.plot(px,py,'bo')
    print stats.spearmanr(px,py)
    plt.show()            
    plt.clf()

def complotstats (xcol,ycol,doannot=False):
    f = open('lbcs_stats.log')
    s = [np.array([])]*len(tels)
    for line in f:
        if len(line)<5 or 'RATE' in line:
            continue
        if line[0] == 'L':
            source, ra, dec, day, time = line.split()
        if line[:5] in tels:
            new = np.asarray(line.split()[1:],dtype='f')
            idx = tels.index(line[:5])
            try:
                s[idx] = np.vstack((s[idx],new))
            except:
                s[idx] = np.copy(new)

    for i in range(len(tels)):
        plt.plot(s[i][:,xcol],s[i][:,ycol],color=telc[i],marker='.',\
                 linestyle='',ms=0.7)
    plt.xlabel(tlabel[xcol],fontsize=16)
    plt.ylabel(tlabel[ycol],fontsize=16)
    plt.xlim(tlim[xcol]); plt.ylim(tlim[ycol])
    if doannot:
        plt.text(50.0,5.0,'Good calibrators',fontsize=16)
        plt.text(150.0,100.0,'Non-calibrators',fontsize=16)
        plt.text(10,20,'P',fontsize=16)
        plt.text(85,40,'S',fontsize=16)
        plt.text(280,100,'X',fontsize=16)
        plt.text(30,80,'D',fontsize=16)
        plt.text(210,8,'D',fontsize=16)
        plt.plot(ppath[:,0],ppath[:,1],'k-')
        plt.plot(spath[:,0],spath[:,1],'k-')
        plt.plot(tpath[:,0],tpath[:,1],'k-')
        plt.plot(upath[:,0],upath[:,1],'k-')
    plt.show()

def lstats (picdir='./picfiles/',logfile='lbcs_stats.log',\
           sumfile='lbcs_stats.temp',access='w'):
    piclist = np.sort(glob.glob(picdir+'*.pic'))
    rstats = np.ones((nstats,len(tels),len(piclist)))*np.nan
    fo = open(logfile,access)
    fs = open(sumfile,access)
    tlist = np.array([''])
    for i in range(len(piclist)):
        if not i%1000:
            print 'Processing line',i
        a = np.load(piclist[i])
        wav = CORE2EFF/np.sqrt((a['uvw'][0]**2).sum())
        t = a['time_start']*24
        th = int(t); t-=th; t*=60.
        tm = int(t); t-=tm; tstr = '%02d:%02d:%02d'%(th,tm,int(t*60.))
        idstr = a['date_obs']+'  '+tstr
        if idstr not in tlist:
            tlist = np.append (tlist,idstr)
        nstr = piclist[i].split('L')[-1].split('.')[0]
        zstr = 'L'+nstr+'  '+a['point_RA']+'  '+a['point_dec']+'  '+idstr
        zgood = ['-']*len(tels)
        fo.write('\n\n'+zstr+'\n')
        fs.write(zstr+'  ')
        fo.write(titlestring)
        for j in range(len(a['delay_4'])):   # loop over telescopes
            jt = tels.index(a['antennas'][j])
            bproj = np.hypot(a['uvw'][j][0],a['uvw'][j][1])/(1000.*wav)
            if np.isnan(bproj):   # telescope not working
                for k in range(nstats):
                    rstats[k,jt,i] = np.nan
                continue
            dcond = ~np.isnan(a['delay'][j,:,0]) & ~np.isnan(a['delay'][j,:,1])
            delay0 = a['delay'][j,:,0][dcond]
            delay1 = a['delay'][j,:,1][dcond]
            phcond = ~np.isnan(a['phase'][j,:,0]) & ~np.isnan(a['phase'][j,:,1])
            phase0 = np.unwrap(a['phase'][j,:,0][phcond])
            phase1 = np.unwrap(a['phase'][j,:,1][phcond])
            rstats[0,jt,i] = np.std(delay0)
            rstats[1,jt,i] = np.median(abs(np.gradient(delay0-delay1))) if len(delay0)>9 else np.nan
            rstats[2,jt,i] = np.abs(a['delay_4'][j,0]-a['delay_4'][j,1])
            rstats[3,jt,i] = 1.E12*np.nanstd(a['rate'][j,:,0])
            rstats[4,jt,i] = 1.E12*(np.abs(a['rate_4'][j,0]-a['rate_4'][j,1]))
            rstats[5,jt,i] = np.nanmin(a['snr'][j,:,0])
            rstats[6,jt,i] = a['snr_4'][j,0]
            rstats[7,jt,i] = np.nan
            if len(phase0)>9:
                diffphgdt = abs(np.gradient(phase0-phase1))*180./np.pi
                rstats[7,jt,i] = np.median(diffphgdt)
# 2-d numpy array of vertices, then bbPath = mpath.Path(array), then bbPath.contains_point((x,y))
            if pPath.contains_point((rstats[1,jt,i],rstats[7,jt,i])):
                zgood[jt]='P'
            elif sPath.contains_point((rstats[1,jt,i],rstats[7,jt,i])):
                zgood[jt]='S'
            elif (tPath.contains_point((rstats[1,jt,i],rstats[7,jt,i])) or
                  uPath.contains_point((rstats[1,jt,i],rstats[7,jt,i]))):
                zgood[jt]='D'
            else:
                zgood[jt]='X'
            fo.write(a['antennas'][j]+'  ')
            for k in range(nstats):  # loop over statistics
                fo.write('%7.1f '%(rstats[k,jt,i]))
            fo.write(' %4.0f %d\n'%(bproj,len(delay0)))
        for k in zgood:
            fs.write('%c'%k)
        sra = np.asarray(a['point_RA'].split(':'),dtype='f')
        fs.write(' %10.6f' % (15.*sra[0]+sra[1]/4.+sra[2]/240.))
        sdec = np.asarray(a['point_dec'].split(':'),dtype='f')
        fs.write(' %10.6f\n' % (sdec[0]+sdec[1]/60.+sdec[2]/3600.))
    fo.close()
    fs.close()

#---------------------------------------------
# Combine the pic files and extract the statistics
def getstats (picdir='/data/data-disk1/data/LBCS/picfiles'):
    lstats(picdir)
    print picdir 
    os.system('sort -k 4,5 lbcs_stats.temp >lbcs_stats.temp1')
    # Quality control: add quality flag if fewer than 8 (4) sources detected
    # by at least 2 of DE601, DE605, DE609
    a = np.loadtxt('lbcs_stats.temp1',dtype='S')
    fo = open('lbcs_stats.temp2','w')
    icou=0
    while True:
        try:
            b=a[(a[:,3]==a[icou,3])&(a[:,4]==a[icou,4])]
        except:
            break
        qthis = np.zeros(len(b))
        for i in range(len(b)):
            if (b[i,5][0]+b[i,5][4]+b[i,5][8]).count('P')>1:
                qthis[i]=1
        for ib in b:
            fo.write(ib[0]+'  '+ib[1]+'  '+ib[2]+'  '+ib[3]+'  '+ib[4]+'  '+ib[5]+\
                     '  '+'%02d'%int(100.*qthis.sum()/len(b))+'  '+ib[6]+'  '+ib[7]+'\n')
        icou += len(b)
    fo.close()
    os.system('sort -k 2 lbcs_stats.temp2 >lbcs_stats.temp3')   # sort to RA
    os.system('cat lbcs_stats.head lbcs_stats.temp3 >lbcs_stats.sum')  # add header
    os.system('rm lbcs_stats.temp1;rm lbcs_stats.temp2;rm lbcs_stats.temp3')

# Get coordinates from the sum file (columns 2 and 3)
def sum2coords (a):
    for obs in a:
        new = np.array([hms2decimal(obs[1],':'),dms2decimal(obs[2],':')])
        try:
            coords = np.vstack((coords,new))
        except:
            coords = np.copy(new)
    return coords

# Look for duplicate observations and produce a list of two-character 
# representations from pairs of observations - e.g. 'PP' means detected
# on a particular baseline in both
def reproducibility ():   
    a = np.loadtxt('lbcs_stats.sum',dtype='S')
    coords = sum2coords(a)
    acorr = correlate(coords,0,1,coords,0,1,0.003)
    acorr = acorr[acorr[:,0]<acorr[:,1]]
    pp=ps=px=sx=xx=0; rr=[]
    for i in acorr:
        r1,r2 = a[int(i[0])][5], a[int(i[1])][5]
        for j in range(len(r1)):
            print r1[j]+r2[j]
            rr.append(r1[j]+r2[j])

def plotit():
    plotdensity(0,'P',False,'lbcs_den_0_P.png','DE601/P')
    plotdensity(0,'PS',False,'lbcs_den_0_PS.png','DE601/PS') 
    plotdensity(7,'P',False,'lbcs_den_0_P.png','UK608/P') 
    plotdensity(7,'PS',False,'lbcs_den_7_PS.png','UK608/PS') 
   
    #plotdensity(0,'X',False,'lbcs_den_0_X.png','DE601/X') 
    #plotdensity(0,'PSX',False,'lbcs_den_0_PSX.png','DE601/PSX') 
    #plotdensity(0,'S',False,'lbcs_den_0_S.png','DE601/S') 
    #plotdensity(7,'X',False,'lbcs_den_7_X.png','UK608/X') 
    #plotdensity(7,'PSX',False,'lbcs_den_7_PSX.png','UK608/PSX)
    #plotdensity(7,'S',False,'lbcs_den_7_S.png','UK608/S') 

def plotdensity (station=0,reqtype='PSX-',dobar=True,\
                 outfile='lbcs_den.png',text='',vmin=0,vmax=2):
    radius = 3.0
    ain = np.loadtxt('lbcs_stats.sum',dtype='S')
    for i in ain:
        if i[5][station] in reqtype:
            try:
                a = np.vstack((a,i))
            except:
                a = np.copy(i)
          
    coords = sum2coords(a)
    y,x = np.meshgrid(np.arange(0.,90,2),np.arange(0.,360,2))
    cgrid = np.dstack((x,y)).reshape(45*180,2)
    acorr = correlate(cgrid,0,1,coords,0,1,radius)
    ngrid = 0.0*cgrid[:,0]
    for i in range(len(cgrid)):
        ngrid[i] = len(np.argwhere(acorr[:,0]==i))
    cplot = ngrid.reshape(180,45).T/(np.pi*radius**2)
    cplot = cplot[:,::-1]
    cplot = cplot[1::,:]
    plt.imshow(np.flip(cplot, axis=0),extent=[24.0,0.0,0,90],cmap=matplotlib.cm.gray_r, vmin=vmin,vmax=vmax,aspect=0.0692)
    plt.contour(cplot,[0.5, 1.0],extent=[24.0,0.0,0,90]) ## Added contour of 0.5 calibrators for square degree
    plt.text(23.8,5,text,{'color':'black','fontproperties':'bold',\
	     'bbox':dict(boxstyle="square", facecolor='white')})
    eq = ga2eq(10); plt.plot(eq[:,0],eq[:,1],'r-')
    eq = ga2eq(-10); plt.plot(eq[:,0],eq[:,1],'r-')
    plt.xlabel ('Right ascension/hr')
    plt.ylabel ('Declination/deg')
    plt.ylim(0.0,90.0)
    plt.grid()
    if dobar:
        plt.colorbar(orientation='horizontal')
    plt.savefig (outfile,bbox_inches='tight')
    plt.clf()

def ga2eq (ga):
    b = np.deg2rad(ga)
    pole_ra,pole_dec,posangle = 3.366033,0.473478,0.574772
    eq = np.zeros((360,2))
    l = np.arange(0.0,2.*np.pi,2.*np.pi/360.0)
    eq[:,0] = np.arctan2((np.cos(b)*np.cos(l-posangle)), (np.sin(b)*np.cos(pole_dec) - np.cos(b)*np.sin(pole_dec)*np.sin(l-posangle))) + pole_ra
    np.putmask(eq[:,0],eq[:,0]>2*np.pi,eq[:,0]-2.*np.pi)
    eq[:,1] = np.arcsin(np.cos(b)*np.cos(pole_dec)*np.sin(l-posangle)+np.sin(b)*np.sin(pole_dec))
    eq = eq[np.argsort(eq[:,0])]
    eq[:,0] *= 24./(2.*np.pi)
    eq[:,1] *= 180./np.pi
    return eq   
   
# Make plots of detectability against various things
# Columns: lbcs_list(N): name ra dec vlss msss c6 wenss isvlba
#          lbcs_south: name ra dec vlss nvss
def plotdetect (dfile = 'lbcs_stats.sum'):
    a = np.loadtxt(dfile ,dtype='S')
    coords = sum2coords(a)
    listN = np.loadtxt('/data/data-disk1/data/LBCS/scheduling/lbcs_list',dtype='S')
    coordsN = np.asarray(listN[:,1:3],dtype='float')
    acorrN = correlate(coords,0,1,coordsN,0,1,0.003)
#    listS = np.loadtxt('../scheduling/lbcs_south',dtype='S')
#    coordsS = np.asarray(listS[:,1:3],dtype='float')
#    acorrS = correlate(coords,0,1,coordsS,0,1,0.003)
    ndet,wflux,spind,ncorr = np.array([]),np.array([]),np.array([]), 0
    det = np.zeros((9,len(acorrN)),dtype='S')
    for cN in acorrN:
        for j in range(9):
            det[j,ncorr] = a[int(cN[0])][5][j]
        ncorr+=1
        ndet = np.append(ndet,a[int(cN[0])][6])
        lineN = listN[int(cN[1])]
        fv,fm,fc,fw = np.asarray(lineN[3:7],dtype='f')
        sp_vc = np.log10(fv/fc)/np.log10(70/150.)
        sp_vw = np.log10(fv/fw)/np.log10(70/325.)
        sp_cw = np.log10(fc/fw)/np.log10(150/325.)
        spind = np.append(spind,np.nanmean([sp_vc,sp_vw,sp_cw]))
        wflux = np.append(wflux,fw)
    plt.semilogx(wflux,spind,'bx')
    plt.semilogx(wflux[det[0]=='P'],spind[det[0]=='P'],'rx')
    plt.semilogx(wflux[det[7]=='P'],spind[det[7]=='P'],'gx')
    aplot0 = np.zeros((10,10))*np.nan
    aplot7 = np.zeros((10,10))*np.nan
    n = len(wflux)
    for i in range(10):
        x = -1.0+0.2*i
        for j in range(10):
            nall = n0 = n7 = 0.0
            y = -1.0+0.1*j
            for k in range(n):
                if x<np.log10(wflux[k])<x+0.2 and y<spind[k]<y+0.1:
                    nall+=1.0
                    if det[0,k]=='P':
                        n0+=1.0
                    if det[7,k]=='P':
                        n7+=1.0
            if nall and n0 and (n0/nall)*np.sqrt(1./n0+1./nall)<0.1:
                aplot0[i,j] = n0/nall
            if nall and n7 and (n7/nall)*np.sqrt(1./n7+1./nall)<0.1:
                aplot7[i,j] = n7/nall
            if i==0 and j==1:
                print 'n0,n7,nall',n0,n7,nall

    matplotlib.rcParams.update({'font.size': 14})
    plt.subplot(211,xticks=[])                
    plt.imshow(np.flip(aplot0, axis=0),extent=[-1.0,1.0,-1.0,0.0],vmax=0.6,cmap=matplotlib.cm.gray_r)
    plt.colorbar()
    plt.ylabel('Low-f spectral index')
    plt.plot([-1.,0.0],[-0.5,-1],'k-')
    plt.text(0.5,-0.2,'DE601')
    plt.subplot(212)
    plt.colorbar()
    plt.ylabel('Low-f spectral index')
    plt.plot([-1.,0.0],[-0.5,-1],'k-')
    plt.xlabel('log 325-MHz flux density/Jy')
    plt.imshow(np.flip(aplot7, axis=0),extent=[-1.0,1.0,-1.0,0.0],vmax=0.6,cmap=matplotlib.cm.gray_r)
    plt.text(0.5,-0.2,'UK608')
#    plt.show()
    #plt.savefig('lbcs_detfig.png',bbox_inches='tight') # Kaspars
    plt.savefig('lbcs_detfig.png',dpi=300) # Kaspars

def plotstats(doshow=False):
    #  next part of code is public-release along with lbcs_stats.sum
    a = np.loadtxt('lbcs_stats.sum',dtype='S')
    ra,dec,status,q = np.zeros(len(a)),np.zeros(len(a)),np.zeros(len(a)),a[:,6]
    m = []
    colour = []
    for i in range(len(a)):
        t1,t2,t3 = a[i,1].split(':')
        ra[i] = float(t1)*15. + float(t2)/4. + float(t3)/240.
        t1,t2,t3 = a[i,2].split(':')
        dec[i] = float(t1)+float(t2)/60.+float(t3)/3600.
        if '-' in a[i,2]:
            dec[i]*=-1.
        # exclude PL stations as present in minority of obs
        Np,Nx,Ns,No = a[i,5][:9].count('P'), a[i,5][:9].count('X'),\
                      a[i,5][:9].count('S'), a[i,5][:9].count('-')
        Nt = Np+Ns+Nx
        if Np>=Nt-1 and Np>=4:
            m.append('o')
            colour.append ('g')
        elif Np+Ns>=Nt-1 and Np+Ns>=4:
            m.append('o')
            colour.append ( 'y')
        elif Np+Ns>=Nt-2 and Np+Ns>=3:
            m.append('o')
            colour.append ('orange')
        else:
            m.append('x')
            if int(q[i]) > 20:
                colour.append ('b')
            elif int(q[i]) > 10:
                colour.append ('m')
            else:
                colour.append ('r')
    
    m=np.asarray(m)
    colour=np.asarray(colour)
    fig = plt.figure(figsize=(12.02,3.0))
    ra/=15.0   # hours
    plt.scatter(ra[m=='o'],dec[m=='o'],color=colour[m=='o'],marker='o')
    plt.scatter(ra[m=='x'],dec[m=='x'],color=colour[m=='x'],marker='x')
    plt.grid()
    plt.xlabel('Right ascension/hours')
    plt.ylabel('Declination/degrees')
    plt.xlim(-0.1,24.1);plt.ylim(0.0,90.0)
    plt.plot([19.98,23.38],[40.4,58.8],'k*',ms=12)
    plt.xlim(plt.xlim()[::-1])
    plt.show() if doshow else plt.savefig('lbcs_stats.png',bbox_inches='tight')
    plt.clf()
    ra_g, dec_g = ra[colour=='g'], dec[colour=='g']
    plt.scatter(ra_g,dec_g,color='g',marker='o')
    plt.xlim(-0.1,24.1);plt.ylim(0.0,90.0)
    plt.xlim(plt.xlim()[::-1])
    plt.grid()
    plt.show() if doshow else plt.savefig('lbcs_green.png')
    plt.clf()
    plt.scatter(ra[m=='o'],dec[m=='o'],color=colour[m=='o'],marker='o')
    plt.scatter(ra[m=='x'],dec[m=='x'],color=colour[m=='x'],marker='x')
    plt.xlim(100./15,120./15);plt.ylim(43.58,56.42)
    plt.grid()
    plt.xlabel('Right ascension/hours')
    plt.ylabel('Declination/degrees')
    plt.xlim(plt.xlim()[::-1])
    plt.show() if doshow else plt.savefig('lbcs_stats1.png',bbox_inches='tight')
    plt.clf()
    
def elev (a, tel):
    Y,M,D = np.asarray(a['date_obs'].split('-'),dtype='int')
    if M<3:
        Y -= 1
        M += 12
    A = int(Y/100)
    B = 2-A+int(A/4)
    C = int(365.25*float(Y))
    E = int(30.6001*(M+1.0))
    JD = B+C+D+E+1720994.5
    T = (JD-2451545.0)/36525.0
    UT = 24.0*a['time_start']
    GST = 6.697374558+(2400.051336*T)+(0.000025862*T**2)+(UT*1.0027379093)
    while GST<0.0:
        GST+=24.0
    while GST>24.0:
        GST-=24.0
    ibpos = np.loadtxt('IBpos',dtype='S')
    lat = float(ibpos[:,1][np.argwhere(ibpos[:,0]==tel)[0][0]])
    longit = float(ibpos[:,2][np.argwhere(ibpos[:,0]==tel)[0][0]])
    lst = 15.0*GST+longit
    rastr = a['point_RA'].split(':')
    ra = 15.*float(rastr[0])+float(rastr[1])/4.0+float(rastr[2])/240.0
    decstr = a['point_dec'].split(':')
    dec = float(decstr[0])+float(decstr[1])/60.0+float(decstr[2])/3600.0
    ha = np.deg2rad(lst-ra)
    lat = np.deg2rad(lat)
    dec = np.deg2rad(dec)
    return np.rad2deg(np.arcsin(np.sin(lat)*np.sin(dec)+np.cos(lat)*np.cos(dec)*np.cos(ha)))
