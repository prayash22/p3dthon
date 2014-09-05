from matplotlib.collections import LineCollection                                                                                                          
from matplotlib.colors import Normalize

whichflg = 'll' 
#diag_type = 'single'
diag_type = 'multi'

if 'fig' not in locals():
    fig = figure()
fig.clf()


if diag_type == 'multi':
    ax= fig.add_subplot(211)
    ax.imshow(CR['tiparav']-mean(CR['tiparav'][0,:]),origin='low',extent=extent,cmap='bwr',vmin=-2.,vmax=2.)
    ax.autoscale(False)
    for p in range(TR._npart):
        print 'adding trajectory %d...' % p

        x = TR.r[0,p,::100]
        y = TR.r[1,p,::100]

        points = np.array([x, y]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)

        vmin = .08*p 
        vmax = 2.8*(p+1) 

        lc = LineCollection(segments,
                            cmap='cool')
                            #norm=Normalize(vmin,vmax))
        lc.set_array(x) 
        lc.set_linewidth(1)
        ax.add_collection(lc)

    ax = fig.add_subplot(223)
    p=0
    v_mag = sqrt(TR.v[0,p,::100]**2+TR.v[1,p,::100]**2+TR.v[2,p,::100]**2)
    timevec = arange(TR.tstart,TR.tend+TR.dt,100*TR.dt)
    #ax.plot(timevec,v_mag)

    cVD = VD['i'][1]
    Hmasked = np.ma.masked_where(cVD[1]==0,cVD[1])
    pcm = pcolormesh(cVD[2],cVD[3],Hmasked)
    ax.autoscale(False)
    ax.plot(arange(-200,200),0.*arange(-200,200),'k--')
    ax.plot(0.*arange(-200,200),arange(-200,200),'k--')
    ax.set_title(cVD[0])
    for p in range(TR._npart):
        ax.text(TR.v[0,p,int(TR.t0/TR.dt)],TR.v[2,p,int(TR.t0/TR.dt)],str(p))
        ax.set_xlim([-6,6])
        ax.set_ylim([-6,6])


    ax = fig.add_subplot(224)
    cVD = VD['i'][2]
    Hmasked = np.ma.masked_where(cVD[1]==0,cVD[1])
    pcm = pcolormesh(cVD[2],cVD[3],Hmasked)
    ax.autoscale(False)
    ax.plot(arange(-200,200),0.*arange(-200,200),'k--')
    ax.plot(0.*arange(-200,200),arange(-200,200),'k--')
    ax.set_title(cVD[0])
    for p in range(TR._npart):
        ax.text(TR.v[1,p,int(TR.t0/TR.dt)],TR.v[2,p,int(TR.t0/TR.dt)],str(p))
        ax.set_xlim([-6,6])
        ax.set_ylim([-6,6])

#c#         ax.plot(TR.r[0,p,:],TR.r[1,p,:])


if diag_type == 'single':
    which_part = 0

    ax= fig.add_subplot(321)
    ax.imshow(CR['tiparav']-mean(CR['tiparav'][0,:]),origin='low',extent=extent,cmap='bwr',vmin=-2.,vmax=2.)
    ax.autoscale(False)
#ax.plot(TR.r[0,0,:],TR.r[1,0,:],'b')
    ax.plot(TR.r[0,which_part,:],TR.r[1,which_part,:],'g')

    ax = fig.add_subplot(323)
    ax.imshow(CR['tiparav']-mean(CR['tiparav'][0,:]),origin='low',extent=[extent[0],extent[1],extent[2]-extent[3]+10.,10.],alpha=0.0)
    ax.autoscale(False)
    ax.plot(TR.r[0,which_part,:],TR.r[2,which_part,:],'g')
#ax.plot(TR.r[0,0,:],TR.r[2,0,:],'b')
    ax.set_title('x vx z')

    ax = fig.add_subplot(325)
    ax.plot(TR.r[1,which_part,:],TR.r[2,which_part,:],'g')
#ax.plot(TR.r[1,0,:],TR.r[2,0,:],'b')
    ax.set_title('y vx z')

    ax = fig.add_subplot(322)
    timevec = arange(TR.tstart,TR.tend+TR.dt,TR.dt)
    ax.plot(timevec,TR.r[0,which_part,:],'k')
    ax.plot(timevec,TR.r[1,which_part,:],'r')
    ax.plot(timevec,TR.r[2,which_part,:],'b')
    ax.set_title('x (k), y(r), z(b)')

    ax = fig.add_subplot(324)
    timevec = arange(TR.tstart,TR.tend+TR.dt,TR.dt)
    ax.plot(timevec,TR.v[0,which_part,:],'k')
    ax.plot(timevec,TR.v[1,which_part,:],'r')
    ax.plot(timevec,TR.v[2,which_part,:],'b')
    ax.set_title('v_x (k), v_y(r), v_z(b)')

    ax = fig.add_subplot(326)
    timevec = arange(TR.tstart,TR.tend+TR.dt,TR.dt)
    vmag = sqrt(sum(TR.v[:,which_part,:]**2,0))
    vpar = abs(sum(TR.v[:,which_part,:]*TR.Bp[:,which_part,:],0)/sqrt(sum(TR.Bp[:,which_part,:]**2,0)))
    vperp = sqrt(vmag**2 - vpar**2)
    ax.plot(timevec,.5*vmag**2,'k')
    ax.plot(timevec,.5*vpar**2,'r')
    ax.plot(timevec,.5*vperp**2,'b')
    ax.set_title('|v|^2 (k), v_par^2(r), v_perp^2(b)')



#c# if whichflg == 'll':
#c#     ax= fig.add_subplot(321)
#c#     ax.imshow(CR['tiparav']-mean(CR['tiparav'][0,:]),origin='low',extent=extent,cmap='bwr',vmin=-2.,vmax=2.)
#c#     ax.autoscale(False)
#c#     #ax.plot(TRlu.r[0,0,:],TRlu.r[1,0,:],'b')
#c#     ax.plot(TRll.r[0,0,:],TRll.r[1,0,:],'g')
#c# 
#c#     ax = fig.add_subplot(323)
#c#     ax.imshow(CR['tiparav']-mean(CR['tiparav'][0,:]),origin='low',extent=[extent[0],extent[1],extent[2]-extent[3]+10.,10.],alpha=0.0)
#c#     ax.autoscale(False)
#c#     ax.plot(TRll.r[0,0,:],TRll.r[2,0,:],'g')
#c#     #ax.plot(TRlu.r[0,0,:],TRlu.r[2,0,:],'b')
#c#     ax.set_title('x vx z')
#c# 
#c#     ax = fig.add_subplot(325)
#c#     ax.plot(TRll.r[1,0,:],TRll.r[2,0,:],'g')
#c#     #ax.plot(TRlu.r[1,0,:],TRlu.r[2,0,:],'b')
#c#     ax.set_title('y vx z')
#c# 
#c#     ax = fig.add_subplot(322)
#c#     timevec = arange(TRll.tstart,TRll.tend+TRll.dt,TRll.dt)
#c#     ax.plot(timevec,TRll.r[0,0,:],'k')
#c#     ax.plot(timevec,TRll.r[1,0,:],'r')
#c#     ax.plot(timevec,TRll.r[2,0,:],'b')
#c#     ax.set_title('x (k), y(r), z(b)')
#c# 
#c#     ax = fig.add_subplot(324)
#c#     timevec = arange(TRll.tstart,TRll.tend+TRll.dt,TRll.dt)
#c#     ax.plot(timevec,TRll.v[0,0,:],'k')
#c#     ax.plot(timevec,TRll.v[1,0,:],'r')
#c#     ax.plot(timevec,TRll.v[2,0,:],'b')
#c#     ax.set_title('v_x (k), v_y(r), v_z(b)')
#c# 
#c#     ax = fig.add_subplot(326)
#c#     timevec = arange(TRll.tstart,TRll.tend+TRll.dt,TRll.dt)
#c#     vmag = sqrt(sum(TRll.v[:,0,:]**2,0))
#c#     vpar = abs(sum(TRll.v[:,0,:]*TRll.Bp[:,0,:],0)/sqrt(sum(TRll.Bp[:,0,:]**2,0)))
#c#     vperp = sqrt(vmag**2 - vpar**2)
#c#     ax.plot(timevec,.5*vmag**2,'k')
#c#     ax.plot(timevec,.5*vpar**2,'r')
#c#     ax.plot(timevec,.5*vperp**2,'b')
#c#     ax.set_title('|v|^2 (k), v_par^2(r), v_perp^2(b)')
#c# 
#c# elif whichflg == 'lu':
#c#     ax= fig.add_subplot(321)
#c#     ax.imshow(CR['tiparav']-mean(CR['tiparav'][0,:]),origin='low',extent=extent,cmap='bwr',vmin=-2.,vmax=2.)
#c#     ax.autoscale(False)
#c#     ax.plot(TRlu.r[0,0,:],TRlu.r[1,0,:],'b')
#c#     #ax.plot(TRll.r[0,0,:],TRll.r[1,0,:],'g')
#c# 
#c#     ax = fig.add_subplot(323)
#c#     ax.imshow(CR['tiparav']-mean(CR['tiparav'][0,:]),origin='low',extent=[extent[0],extent[1],extent[2]-extent[3]+20.,20.],alpha=0.0)
#c#     ax.autoscale(False)
#c#     #ax.plot(TRll.r[0,0,:],TRll.r[2,0,:],'g')
#c#     ax.plot(TRlu.r[0,0,:],TRlu.r[2,0,:],'b')
#c#     ax.set_title('x vx z')
#c# 
#c#     ax = fig.add_subplot(325)
#c#     #ax.plot(TRll.r[1,0,:],TRll.r[2,0,:],'g')
#c#     ax.plot(TRlu.r[1,0,:],TRlu.r[2,0,:],'b')
#c#     ax.set_title('y vx z')
#c# 
#c# 
#c#     ax = fig.add_subplot(322)
#c#     timevec = arange(TRlu.tstart,TRlu.tend+TRlu.dt,TRlu.dt)
#c#     ax.plot(timevec,TRlu.r[0,0,:],'k')
#c#     ax.plot(timevec,TRlu.r[1,0,:],'r')
#c#     ax.plot(timevec,TRlu.r[2,0,:],'b')
#c#     ax.set_title('x (k), y(r), z(b)')
#c# 
#c#     ax = fig.add_subplot(324)
#c#     timevec = arange(TRlu.tstart,TRlu.tend+TRlu.dt,TRlu.dt)
#c#     ax.plot(timevec,TRlu.v[0,0,:],'k')
#c#     ax.plot(timevec,TRlu.v[1,0,:],'r')
#c#     ax.plot(timevec,TRlu.v[2,0,:],'b')
#c#     ax.set_title('v_x (k), v_y(r), v_z(b)')
#c# 
#c#     ax = fig.add_subplot(326)
#c#     timevec = arange(TRlu.tstart,TRlu.tend+TRlu.dt,TRlu.dt)
#c#     vmag = sqrt(sum(TRlu.v[:,0,:]**2,0))
#c#     vpar = abs(sum(TRlu.v[:,0,:]*TRlu.Bp[:,0,:],0)/sqrt(sum(TRlu.Bp[:,0,:]**2,0)))
#c#     vperp = sqrt(vmag**2 - vpar**2)
#c#     ax.plot(timevec,.5*vmag**2,'k')
#c#     ax.plot(timevec,.5*vpar**2,'r')
#c#     ax.plot(timevec,.5*vperp**2,'b')
#c#     ax.set_title('|v|^2 (k), v_par^2(r), v_perp^2(b)')

fig.show()






