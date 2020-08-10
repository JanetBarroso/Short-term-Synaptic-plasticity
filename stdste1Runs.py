from stdste1 import *



# ------------------------------------------------------------
# Alpha function 
if 0: 
    tau_alpha=1.0
    t = sc.arange(0,10,0.001)
    y = alphaFunction(x=t,tau=tau_alpha)
    f= plt.figure()
    plt.ioff()
    ax=f.add_subplot(111)
    ax.plot(t, y) 
    ax.plot([0,t.max()], [0,0],'k:') 
    ax.set_xlabel('ms')
    ax.text(8,0.25,r'$\alpha(t)$')
    for n in range(1,4):
        if n>1:
            str0=r"$%d \tau_{\alpha}$"%(n)
        else:
            str0=r"$\tau_{\alpha}$"%()
            
        ax.annotate(str0,
                xy=(n*tau_alpha,0), xycoords='data',
                xytext=(n*tau_alpha, 0.2), textcoords='data',
                arrowprops=dict(arrowstyle="->",
                                connectionstyle="arc3"))
    ax.set_xlim([-0.1,t.max()])
    ax.set_ylim([-0.1,1.1])
    plt.ion(); plt.draw()
    if 1: 
        f.savefig('./figuresDynamicalSystems/alphaFunctionTau1.png')



def stimTrains(pars):
    #print 'Obtaining presynaptic spike times and stimulus train'
    sampleTimes=sc.arange(pars['timeMin'],pars['timeMax'],pars['timeStep'])
    pars['stimPeriod']=1000/pars['stimHz']
    preAPs=sc.arange(pars['stimStart'], pars['timeMax'], 1000.0/pars['stimHz'])
    alphaTrainP=trainAlpha(samplingTimes=sampleTimes, pulseTimes=preAPs[1:], tauTrain=pars['tauTrain'],downAccel=1.0)
    alphaTrainX=trainAlpha(samplingTimes=sampleTimes, pulseTimes=preAPs, tauTrain=pars['tauTrain'],downAccel=1.0)
    alphaTrainP= alphaTrainP/alphaTrainP.max()
    alphaTrainX= alphaTrainX/alphaTrainX.max()
    pars['alphaTrainP']= alphaTrainP
    pars['alphaTrainX']= alphaTrainX
    pars['spikesCa']=lambda t: sc.interp(t,sampleTimes, alphaTrainP)
    pars['spikesX']=lambda t: sc.interp(t,sampleTimes, alphaTrainX)
    nStim= len(preAPs)
    peakInds= sc.zeros(nStim,'int32')
    for nn in range(0,nStim):
        peakInds[nn] = 1+gr.find(sampleTimes< preAPs[nn]+pars['tauTrain']).max()
    pars['sampleTimes']=sampleTimes
    pars['preAPs']=preAPs
    pars['peakInds']=peakInds
    pars['nStim']= nStim
    return pars









# ------------------------------------------------------------
# Phase plane for p and x alone
pars={'timeMin':0.0, 'timeMax': 1000.0, 'timeStep':1e-3,
      'tauX':10.0,  'tauP':20.0, 
      'asynX':0.9, 'asynP':0.3, 'asynXPm':'ko',
      'kP':0.0, 'kX':1.0, 
      'jumpP':0.4, 'tauTrain':1.0, 
      'stimHz':20.0, 'stimStart':10.0}
pars= stimTrains(pars)
#dxCurve2, dpCurve2, dxpCurve2=stFD([xRange,pRange, pRange*xRange],0,pars)

# Phase plane for p
if 0: 
    pRange=sc.arange(0,1.0001,1e-3)
    dpCurve1= (pRange**pars['kP'])*(pars['asynP']-pRange)/pars['tauP'] 
    dpCurve1AP= dpCurve1 + (1-pRange) * pars['jumpP']
    pPP= plt.figure()
    plt.ioff()
    ax=pPP.add_subplot(111)
    ax.plot([pRange.min(), pRange.max()], [0,0], 'k', alpha=0.5, lw=1)
    ax.plot(pRange, dpCurve1, 'b', lw=2, alpha=1.0, label='Between APs')
    ax.plot(pRange, dpCurve1AP, 'b--', lw=1, alpha=1, label='AP')
    ax.plot([pars['asynP'], ], [0,], 'ko', alpha=1, ms=2)
    pInfAP=(pars['asynP']+ pars['jumpP']*pars['tauP'])/(1 + pars['jumpP']*pars['tauP'])
    ax.plot([pInfAP, ], [0,], 'ko', alpha=1, ms=2)
    ax.set_xlabel(r'$p$')
    ax.set_ylabel(r'$\partial_t p$')
    ax.annotate(r"$p_{\infty}$",
                xy=(pars['asynP'],0), xycoords='data',
                xytext=(pars['asynP'], 0.1), textcoords='data',
                arrowprops=dict(arrowstyle="->",
                                connectionstyle="arc3"))
    ax.annotate(r"$\frac{p_{\infty} - h \tau_p}{1+h \tau_p}$",
                xy=(pInfAP,0), xycoords='data',
                xytext=(pInfAP, 0.1), textcoords='data',
                arrowprops=dict(arrowstyle="->",
                                connectionstyle="arc3"))
    str0=r'$p_\infty$=%.1f, $\tau_p$=%.0f, $h$=%.1f'%(pars['asynP'], pars['tauP'], pars['jumpP'])
    ax.text(0.1,-0.1,str0)
    ax.set_xlim(0,1)
    ax.set_ylim(-0.2,0.4)
    plt.ion(); plt.draw()

# ------------------------------------------------------------
# Solution to linear equation
if 0: 
    pars['tauP']=10.0
    pars['jumpP']=0.9
    p0= (1-pars['asynP']) * pars['jumpP'] ; 
    pars['kP']=0.0
    pRange=sc.arange(0,1.0001,1e-3)
    pSol1Curve=(pars['asynP'] - (pars['asynP']-p0) * sc.exp(-(pars['sampleTimes']) /pars['tauP'] ))
    b=sc.int32(0.8*pars['stimPeriod']/pars['timeStep'])
    pSolFig= plt.figure()
    plt.ioff()
    ax=pSolFig.add_subplot(111)
    ax.plot([pars['sampleTimes'][0], pars['stimStart'] +pars['stimPeriod'] +pars['sampleTimes'][b]], [0,0], 'k', alpha=0.5, lw=1)
    ax.plot([pars['sampleTimes'][0], pars['stimStart'] +pars['stimPeriod'] +pars['sampleTimes'][b]], [pars['asynP'], pars['asynP']], 'k', alpha=0.5, lw=1)
    ax.plot([pars['sampleTimes'][0], pars['stimStart']+pars['stimPeriod']], [pars['asynP'], pars['asynP']], 'b', alpha=1, lw=2)
    ax.plot(pars['stimStart'] +pars['stimPeriod'] + pars['sampleTimes'][:b], pSol1Curve[:b], 'b', lw=2, alpha=1.0)
    ax.annotate(r"$p_{\infty}$", xy=(pars['stimStart'], pars['asynP']), xycoords='data', xytext=(5, 0.5*pars['asynP']), textcoords='data', arrowprops=dict(arrowstyle="->", connectionstyle="arc3"))
    ax.annotate(r"$p_1$", xy=(pars['stimStart'] +pars['stimPeriod'] ,p0), xycoords='data', xytext=(pars['stimPeriod'], p0), textcoords='data', arrowprops=dict(arrowstyle="->", connectionstyle="arc3"))
    str0=r'$p_\infty$=%.1f, $\tau_p$=%.0f, $h$=%.1f'%(pars['asynP'], pars['tauP'], pars['jumpP'])
    ax.text(30,0.1,str0)
    str1=r'APs : $t_0$=%.f, $t_1$=%.f'%(pars['stimStart'], pars['stimStart']+pars['stimPeriod'])
    ax.text(30,0.2,str1)
    ax.set_xlabel('Milliseconds')
    ax.set_ylabel('P(release)')
    ax.set_ylim(-0.01,1.3*sc.maximum(pars['asynP'], p0))
    plt.ion(); plt.draw()

# ------------------------------------------------------------
# Phase plane for x
if 0: 
    p0=0.5
    pars['tauX']=0.9*pars['asynX']/p0
    xRange=sc.arange(0,1.0001,1e-3)
    dxCurve1= (xRange**pars['kX'])*(pars['asynX']-xRange)/pars['tauX'] 
    dxCurve1AP= dxCurve1 - xRange * p0
    xPP= plt.figure()
    plt.ioff()
    ax=xPP.add_subplot(111)
    ax.plot([xRange.min(), xRange.max()], [0,0], 'k', alpha=0.5, lw=1)
    ax.plot(pRange, dxCurve1, 'b', lw=2, alpha=1.0, label='Between APs')
    ax.plot(pRange, dxCurve1AP, 'b--', lw=1, alpha=1, label='AP')
    ax.plot([pars['asynX'], ], [0,], 'ko', alpha=1, ms=2)
    xInfAP=(pars['asynX'] - p0*pars['tauX'])
    ax.plot([xInfAP, ], [0,], 'ko', alpha=1, ms=2)
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$\partial_t x$')
    ax.annotate(r"$x_{\infty}$",
                xy=(pars['asynX'],0), xycoords='data',
                xytext=(pars['asynX'], 0.3), textcoords='data',
                arrowprops=dict(arrowstyle="->",
                                connectionstyle="arc3"))
    ax.annotate(r"$x_{\infty}-p \tau_x$",
                xy=(xInfAP,0), xycoords='data',
                xytext=(xInfAP, 0.3), textcoords='data',
                arrowprops=dict(arrowstyle="->",
                                connectionstyle="arc3"))
    str0=r'$x_\infty$=%.1f, $\tau_x$=%.0f, $p$=%.1f'%(pars['asynX'], pars['tauX'], p0)
    ax.text(0.1,-0.1,str0)
    ax.set_xlim(0,1)
    ax.set_ylim(-0.2,0.4)
    plt.ion(); plt.draw()

# ------------------------------------------------------------
# Solution to logistic equation
if 0: 
    p0=0.2
    pars['tauX']=10.0
    xRange=sc.arange(0,1.0001,1e-3)
    x0= p0*pars['asynX']
    b0=sc.int32(pars['stimPeriod']/pars['timeStep'])
    denom0=(x0 + (pars['asynX']-x0) * sc.exp(-pars['sampleTimes'][:b0] * pars['asynX']/pars['tauX'] ))
    xSol0Curve=x0*pars['asynX'] / denom0
    x1= p0*xSol0Curve[-1]
    b1=sc.int32(0.8*pars['stimPeriod']/pars['timeStep'])
    denom1=(x1 + (pars['asynX']-x1) * sc.exp(-pars['sampleTimes'][:b1] * pars['asynX']/pars['tauX'] ))
    xSol1Curve=x1*pars['asynX'] / denom1
    xSol= plt.figure()
    plt.ioff()
    ax=xSol.add_subplot(111)
    ax.plot([pars['sampleTimes'][0], pars['stimStart']+1.8*pars['stimPeriod']], [0,0], 'k', alpha=0.5, lw=1)
    ax.plot([pars['sampleTimes'][0], pars['stimStart']+1.8*pars['stimPeriod']], [pars['asynX'], pars['asynX']], 'k', alpha=0.5, lw=1)
    ax.plot([pars['sampleTimes'][0], pars['stimStart']], [pars['asynP'], pars['asynP']], 'b', alpha=1, lw=2)
    ax.plot(pars['stimStart']+pars['sampleTimes'][:b0], xSol0Curve[:b0], 'b', lw=2, alpha=1.0, label='Between APs')
    ax.plot(pars['stimStart']+pars['stimPeriod']+pars['sampleTimes'][:b1], xSol1Curve[:b1], 'b', lw=2, alpha=1.0, label='Between APs')
    ax.annotate(r"$x_{\infty}$", xy=(20, pars['asynX']), xycoords='data', xytext=(10, 0.8*pars['asynX']), textcoords='data', arrowprops=dict(arrowstyle="->", connectionstyle="arc3"))
    ax.annotate(r"$x_0$", xy=(pars['stimStart'],x0), xycoords='data', xytext=(1.2*pars['stimStart'], 0.05), textcoords='data', arrowprops=dict(arrowstyle="->", connectionstyle="arc3"))
    ax.annotate(r"$x_1$", xy=(pars['stimStart']+pars['stimPeriod'],x1), xycoords='data', xytext=(1.2*pars['stimStart']+pars['stimPeriod'], 0.05), textcoords='data', arrowprops=dict(arrowstyle="->", connectionstyle="arc3"))
    str0=r'$x_\infty$=%.1f, $\tau_x$=%.0f, $p$=%.1f'%(pars['asynX'], pars['tauX'], p0)
    ax.text(30,0.3,str0)
    str1=r'$\tau_x/x_\infty $=%.1f ms'%(pars['tauX']/pars['asynX'])
    ax.text(30,0.2,str1)
    str2=r'AP: $t_0$=%.f'%(pars['stimStart'] )
    ax.text(30,0.1,str2)
    ax.set_xlabel('Milliseconds')
    ax.set_ylabel('Occupancy')
    ax.set_ylim(-0.1,1.2*pars['asynX'])
    plt.ion(); plt.draw()

# Check plot sim
if 0:
    pars={'timeMin':0.0, 'timeMax': 1000.0, 'timeStep':1e-3,
      'tauX':75.0,  'tauP':100.0, 
      'asynX':0.9, 'asynP':0.3, 'asynXPm':'ko',
      'kP':0.0, 'kX':1.0, 
      'jumpP':0.4, 'tauTrain':1.0, 
      'stimHz':20.0, 'stimStart':10.0}
    pars= stimTrains(pars)
    orbit=simulateFD(pars)
    r=1; c=1
    f0=gr.figure(figsize=(13,5))
    ax1=f0.add_subplot(r,c,1)
    xpNorm=graphSimulation(ax=ax1, pa=pars, orbit=orbit,saveFig=0)
    ax1.set_ylim(-0.1,1.1*xpNorm.max())
    figName='taux=%.0f_taup=%.0f_h=%.1f_xa=%.1f_pa=%.1f'%(pars['tauX'], pars['tauP'], pars['jumpP'], pars['asynX'], pars['asynP'])
    figStr1=r'$\tau_x$=%.0f, $\tau_p=%.0f, $h$=%.1f, $x_{\infty}$=%.1f, $p_{\infty}$=%.1f'%(pars['tauX'], pars['tauP'], pars['jumpP'], pars['asynX'], pars['asynP'])
    if saveFig:
        f0.savefig('./figuresSynapse/'+figName+'.png')
    print (figName)

# .................................................................
# Ajuste de parametros
# .................................................................
if 1:
    pars1={'timeMin':0.0, 'timeMax': 500.0, 'timeStep':1e-3,
           'tauX':12.75,  'tauP':300.0, 
           'asynX':0.8, 'asynP':0.27, 'asynXPm':'ko',
           'kP':1.0, 'kX':1.0, 
           'jumpP':0.9, 'tauTrain':1.0, 
           'stimHz':20.0, 'stimStart':10.0}
    pars1= stimTrains(pars1)
    sim1=simulateFD(pars1)
    dataVals = subeBaja(x=sc.arange(1,sim1['nStim']+1), a=2.1348, b=3.109, t1=1.6353, t2=1.8831, y0=0.3986)
    #dataVals = subeBaja(x=sc.arange(0,450,1), a=3.09, b=4.7, t1=158.2, t2=37.7, c=-18) #STF Lesion
    f0= ajusteParametros(simData=sim1, expeData=dataVals)




# Check plot sim function
if 0: 
    pars={'timeMin':0.0, 'timeMax': 400.0, 'timeStep':1e-3,
      'tauX':75.0,  'tauP':150.0, 
      'asynX':0.9, 'asynP':0.3, 'asynXPm':'ko',
      'kP':0.0, 'kX':1.0, 
      'jumpP':0.4, 'tauTrain':1.0, 
      'stimHz':20.0, 'stimStart':10.0}

    pars= stimTrains(pars)
    xpFactors=sc.arange(0,10,1.0)
    taus=sc.arange(10,100+1e-7, 10.0)
    r=1; c=1
    f0=gr.figure(figsize=(13,5))
    ax1=f0.add_subplot(111)
    for m in range(len(taus)): 
        for n in range(len(xpFactors)):
            pars['tauX']= taus[m]
            pars['tauP']= xpFactors[n]*pars['tauX']         
            simulation=simulateFD(pars)
            r=1; c=1
            xpNorm=graphSimulation(ax=ax1, pa=simulation, orbit=orbit)
            f0.clf()
        
    

# Roll taus, jumps, and pinf
if 0:
    synapses=list()
    hArray=sc.arange(0.1,1,10*0.15)
    pinfArray=sc.arange(0.1,0.5,10*0.1)
    pars['asynX']=1.0
    pars['jumpP']=0.9
    pars['timeMax']=1000.0
    tauTr=0.1
    pars= stimTrains(pars, tauTrain=tauTr)
    for thisH in hArray:
        pars['jumpP']=thisH
        for thisP in pinfArray:
            pars['asynP']=thisP
            orbit=simulateFD(pars)
            pars['xOrbit']=orbit[0]; 
            pars['pOrbit']=orbit[1]; 
            pars['xpOrbit']=orbit[2]
            rollOverTaus(pars, tauMin=0.0, tauMax=300.0, tauStep=15.0)  
            synapses.append(pars)


# ------------------------------------------------------------
# Generate figures with  of (x,p) in the phase plane for different taus, 
# --------------------------------------------------
if 0:
    pars={'timeMin':0.0, 'timeMax': 500.0, 'timeStep':1e-3,
      'tauX':15.0,  'tauP':30.0, 
      'asynX':0.9, 'asynP':0.3, 'asynXPm':'ko',
      'kP':1.0, 'kX':1.0, 
      'jumpP':0.3, 'tauTrain':1.0, 
      'stimHz':20.0, 'stimStart':10.0}
    startTime=time.time()
    tau=sc.arange(5,75.001,5.0)
    jumps=sc.arange(0.1,1.0,0.1)
    orbits=list()
    pars
    for kP in [0,1.0]:
        pars['kP']=kP
        for kX in [0,1.0]:
            pars['kX']=kX
            for h in jumps:
                pars['jumpP']=h
                for m in range(len(tau)): 
                    pars['tauX']=tau[m]
                    for n in range(len(tau)): 
                        pars['tauP']=tau[n]
                        pars= stimTrains(pars)
                        simulation= simulateFD(pars)
                        fig=phase2Dynamics(simData=simulation,saveDir='jjj')
                        plt.close('all')
                        totTime =time.time()-startTime
                        totTimeMins=totTime/60.0
                        totTimeHrs=totTimeMins/60.0
                        #print "So far this job has taken %f seconds (%f minutes, %f hours) "%(totTime,totTimeMins,totTimeHrs)

    totTime =time.time()-startTime
    totTimeMins=totTime/60.0
    totTimeHrs=totTimeMins/60.0
    #print "Took %f seconds (%f minutes, %f hours) to finish the job"%(totTime,totTimeMins,totTimeHrs)


# -----------------------------------
def phase2DynamicsMovie(simData, pa, msFramePoints=3.1*20, msJumpPoints=10):
    peakInds = pa['peakInds']
    firstPeakInd=peakInds[0]
    relAmp= orbit[2][firstPeakInd]
    xpNorm=orbit[2][peakInds]/relAmp
    figStr1=r'$\tau_x$=%.0f, $\tau_p$=%.0f, $k_x$=%.1f, $k_p$=%.1f, $h$=%.1f, $x_{\infty}$=%.1f, $p_{\infty}$=%.1f'%(pa['tauX'], pa['tauP'], pa['kX'], pa['kP'], pa['jumpP'], pa['asynX'], pa['asynP'])
    figName='taux=%.0f_taup=%.0f_kx=%.1f_kp=%.1f_h=%.2f_xa=%.2f_pa=%.2f_maxT=%d'%(pa['tauX'], pa['tauP'], pa['kX'], pa['kP'], pa['jumpP'], pa['asynX'], pa['asynP'],pa['timeMax'])
    newPath="jjj" #+figName
    #if not os.path.exists(newPath): 
    #        os.makedirs(newPath,mode=0774)
    # Assuming the variable simData is a list 
    nVars=len(orbit)
    nSimPts=len(pa['sampleTimes'])
    #
    nPtsPerFrame= sc.int32(msFramePoints/pa['timeStep'])
    nPtsJump= sc.int32(msJumpPoints/pa['timeStep'])
    nFrames=1+nSimPts/nPtsJump
    lines=list()
    fig1 = plt.figure(figsize=(11,5))
    ax0=fig1.add_subplot(121)
    ax1=fig1.add_subplot(122)
    l0, = ax0.plot([], [], 'b.',ms=1)
    l1, = ax1.plot([], [], 'k.',ms=1,alpha=0.5)
    l2, = ax1.plot([], [], 'k',ms=1,lw=1,alpha=0.5)
    l3, = ax1.plot([], [], 'wo',ms=3)
    ax0.plot([pa['asynX'],], [pa['asynP'],], pa['asynXPm'], ms=4)
    # Comment this line if using with other variables
    ax0.set_xlim(-0.1,1.1);  ax0.set_ylim(-0.1,1.1)
    ax0.set_xlabel(r'$x$'); ax0.set_ylabel(r'$p$');
    ax1.set_xlabel(r'ms'); ax1.set_ylabel(r'$px$');
    ax1.set_xlim(0,pa['timeMax']);  ax1.set_ylim(-0.1,1.1*xpNorm.max())
    for k in range(0,nFrames):
        b=(k+1)*nPtsJump
        a= sc.maximum(b-nPtsPerFrame,0)
        if (b-1>nSimPts):
            b=nSimPts
        l0.set_data(simData[0][a:b], simData[1][a:b])
        l0.set_data(simData[0][a:b], simData[0][a:b]*simData[1][a:b])
        l1.set_data(pa['sampleTimes'][a:b], simData[0][a:b]*simData[1][a:b])
        l2.set_data(pa['sampleTimes'][peakInds], xpNorm)
        l3.set_data(pa['sampleTimes'][peakInds], xpNorm)
        fig1.canvas.draw()
        fig1.savefig(newPath+"/"+figName+"_%03d.png"%k)

    l0.set_data(simData[0], simData[1])
    l0.set_data(simData[0], simData[0]*simData[1])
    l1.set_data(pa['sampleTimes'], simData[0]*simData[1])
    fig1.canvas.draw()
    fig1.suptitle(figStr1, fontsize=10)
    fig1.savefig(newPath+"/"+figName+"_%03d.png"%k)
    print 'Making movie %s.mpg - this make take a while'%movieName
    os.system("mencoder 'mf://%s/%s*.png' -mf type=png:fps=5 -ovc lavc -lavcopts vcodec=wmv2 -oac copy -o %s/%s.mpg"%(newPath,figName,newPath, figName))
    for k in range(0,nFrames):
        os.system("rm '%s/%s*.png'"%(newPath,figName))
    plt.close('all')
    return figName


# --------------------------------------------------
# Movie from simulation of solutions (x,p,xp)(t)
# ------------------------------------------------------------
# Generate movies of (x,p) in the phase plane for different taus, 
# ------------------------------------------------------------
if 0:
    pars={'timeMin':0.0, 'timeMax': 1000.0, 'timeStep':1e-3,
      'tauX':15.0,  'tauP':30.0, 
      'asynX':0.9, 'asynP':0.3, 'asynXPm':'ko',
      'kP':1.0, 'kX':1.0, 
      'jumpP':0.3, 'tauTrain':1.0, 
      'stimHz':20.0, 'stimStart':10.0}

if 0: 
    startTime=time.time()
    tau=sc.arange(5,31,5.0)
    for kP in [0,1.0]:
        pars['kP']=kP
        for kX in [0,1.0]:
            pars['kX']=kX
            for h in sc.arange(0.1,1.0,0.1):
                pars['jumpP']=h
                for m in range(len(tau)): 
                    pars['tauX']=tau[m]
                    for n in range(len(tau)): 
                        pars['tauP']=tau[n]
                        pars= stimTrains(pars)
                        orbit= simulateFD(pars)
                        figName=phase2DynamicsMovie(simData=orbit[0:2], pars=pars, msFramePoints=1.5*20, msJumpPoints=10)
    totTime =time.time()-startTime
    totTimeMins=totTime/60.0
    totTimeHrs=totTimeMins/60.0
    print "Took %f seconds (%f minutes, %f hours) to finish the job"%(totTime,totTimeMins,totTimeHrs)


# -----------------------------------------
def phase1DynamicsMovie(simData, pars, msFramePoints=3.1*20, msJumpPoints=10, movieName='phase1Dynamics'):
    figStr1=r'$\tau_x$=%.0f, $\tau_p$=%.0f, $h$=%.1f, $x_{\infty}$=%.1f, $p_{\infty}$=%.1f'%(pars['tauX'], pars['tauP'], pars['jumpP'], pars['asynX'], pars['asynP'])
    figName='taux=%.0f_taup=%.0f_h=%.1f_xa=%.1f_pa=%.1f'%(pars['tauX'], pars['tauP'], pars['jumpP'], pars['asynX'], pars['asynP'])
    newPath="jjj"+figName
    if len(movieName)<1:
        movieName=figName
    if not os.path.exists(newPath): 
            os.makedirs(newPath,mode=0774)
    # Assuming the variable data is a list 
    nVars=len(orbit)
    nSimPts=len(pars['sampleTimes'])
    #
    nPtsPerFrame= sc.int32(msFramePoints/pars['timeStep'])
    nPtsJump= sc.int32(msJumpPoints/pars['timeStep'])
    nFrames=1+nSimPts/nPtsJump
    lines=list()
    fig1 = plt.figure(figsize=(5,5))
    ax=fig1.add_subplot(111)
    l, = ax.plot([], [], 'b.',ms=1)
    ax.plot([pars['asynX'],], [pars['asynP'],], pars['asynXPm'], ms=4)
    # Comment this line if using with other variables
    ax.set_xlim(-0.1,1.1);  ax.set_ylim(-0.1,1.1)
    ax.set_xlabel(r'$x$'); ax.set_ylabel(r'$p$');
    for k in range(0,nFrames):
        b=(k+1)*nPtsJump
        a= sc.maximum(b-nPtsPerFrame,0)
        if (b-1>nSimPts):
            b=nSimPts
        l.set_data(simData[0][a:b], simData[1][a:b])
        fig1.canvas.draw()
        fig1.savefig(newPath+"/"+figName+"_%03d.png"%k)

    l.set_data(simData[0], simData[1])
    fig1.canvas.draw()
    fig1.suptitle(figStr1, fontsize=10)
    fig1.savefig(newPath+"/"+figName+"_%03d.png"%k)
    print 'Making movie %s.mpg - this make take a while'%movieName
    os.system("mencoder 'mf://%s/%s*.png' -mf type=png:fps=5 -ovc lavc -lavcopts vcodec=wmv2 -oac copy -o %s/%s.mpg"%(newPath,figName,newPath, movieName))
    plt.close('all')
    return figName





def dynamicsMovie(simData, pars, msFramePoints=3.1*20, msJumpPoints=10,cm=['b.','g.','k.']):
    figStr1=r'$\tau_x$=%.0f, $\tau_p$=%.0f, $h$=%.1f, $x_{\infty}$=%.1f, $p_{\infty}$=%.1f'%(pars['tauX'], pars['tauP'], pars['jumpP'], pars['asynX'], pars['asynP'])
    figName='taux=%.0f_taup=%.0f_h=%.1f_xa=%.1f_pa=%.1f'%(pars['tauX'], pars['tauP'], pars['jumpP'], pars['asynX'], pars['asynP'])
    newPath="jjj"+figName
    if not os.path.exists(newPath): 
            os.makedirs(newPath,mode=0774)
    # Assuming the variable data is a list 
    nVars=len(orbit)
    nSimPts=len(pars['sampleTimes'])
    #
    nPtsPerFrame= sc.int32(msFramePoints/pars['timeStep'])
    nPtsJump= sc.int32(msJumpPoints/pars['timeStep'])
    nFrames=1+nSimPts/nPtsJump
    lines=list()
    fig1 = plt.figure()
    ax=fig1.add_subplot(111)
    strLabels=[r'$(t,x)$', r'$(t,p)$', r'$(t,xp)$']
    strLabels=['x', r'$(t,p)$', r'$(t,xp)$']
    for n in range(nVars):
        l, = ax.plot([], [], cm[n],ms=1,label=strLabels[n])
        lines.append(l)
        
    ax.set_xlim(pars['timeMin'],pars['timeMax'])
    ax.set_ylim(-0.1,1.1)
    ax.set_xlabel('time (ms)')
    ax.legend()
    for k in range(0,nFrames):
        b=(k+1)*nPtsJump
        a= sc.maximum(b-nPtsPerFrame,0)
        if (b-1>nSimPts):
            b=nSimPts
        for n in range(nVars):
            lines[n].set_data(pars['sampleTimes'][a:b], simData[n][a:b])
        fig1.canvas.draw()
        fig1.savefig("jjj/fig%03d.png"%k)

    for n in range(nVars):
        lines[n].set_data(pars['sampleTimes'], simData[n])

    fig1.canvas.draw()
    fig1.savefig("jjj/fig%03d.png"%k)
    newPath=r'./jjj'
    if not os.path.exists(newPath): 
            os.makedirs(newPath,mode=0774)
    print 'Making movie animation.mpg - this make take a while'
    os.system("mencoder 'mf://%s/*.png' -mf type=png:fps=5 -ovc lavc -lavcopts vcodec=wmv2 -oac copy -o %s/xpDynamics.mpg"%(newPath,newPath))
    return 

if 0: 
    dynamicsMovie(pars=pars, simData= orbit, msFramePoints=3.1*20, msJumpPoints=10)






# Creating an uncompressed file from all the PNG files in the current directory:
# os.system("mencoder 'mf://*.png' -mf w=800:h=600:fps=10:type=png -ovc raw -oac copy -o output.avi")
# Note: width must be integer multiple of 4, it is a limitation of the RAW RGB AVI format.
# Creating a Motion PNG (MPNG) file from all the PNG files in the current directory:
# mencoder mf://*.png -mf w=800:h=600:fps=25:type=png -ovc copy -oac copy -o output.avi


# -----------------------------------------
# Notes about movies:
# -----------------------------------------
# Creating an MPEG-4 file from all the JPEG files in the current directory:
# mencoder mf://*.jpg -mf w=800:h=600:fps=25:type=jpg -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o output.avi
# Creating an MPEG-4 file from some JPEG files in the current directory:
# mencoder mf://frame001.jpg,frame002.jpg -mf w=800:h=600:fps=25:type=jpg -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o output.avi
# Creating an MPEG-4 file from explicit list of JPEG files (list.txt in current directory contains the list of files to use as source, one per line):
# mencoder mf://@list.txt -mf w=800:h=600:fps=25:type=jpg -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o output.avi
# You can mix different types of images, regardless of the method you use individual filenames, wildcard or file with a list provided of course they have the same dimensions. So you can e.g. take title frame from PNG file, and then put a slideshow of your JPEG photos.
# Creating a Motion JPEG (MJPEG) file from all the JPEG files in the current directory:
# mencoder mf://*.jpg -mf w=800:h=600:fps=25:type=jpg -ovc copy -oac copy -o output.avi
# Creating an uncompressed file from all the PNG files in the current directory:
# mencoder mf://*.png -mf w=800:h=600:fps=25:type=png -ovc raw -oac copy -o output.avi
# Note: width must be integer multiple of 4, it is a limitation of the RAW RGB AVI format.
# Creating a Motion PNG (MPNG) file from all the PNG files in the current directory:
# mencoder mf://*.png -mf w=800:h=600:fps=25:type=png -ovc copy -oac copy -o output.avi
# Creating a Motion TGA (MTGA) file from all the TGA files in the current directory:
# mencoder mf://*.tga -mf w=800:h=600:fps=25:type=tga -ovc copy -oac copy -o output.avi

