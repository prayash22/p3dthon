fig.clf()

whichflg = 'll' 

if 'fig' not in locals():
    fig = figure()
if whichflg == 'll':
    ax= fig.add_subplot(321)
    ax.imshow(CR['tiparav']-mean(CR['tiparav'][0,:]),origin='low',extent=extent,cmap='bwr',vmin=-2.,vmax=2.)
    ax.autoscale(False)
    #ax.plot(TRlu.r[0,0,:],TRlu.r[1,0,:],'b')
    ax.plot(TRll.r[0,0,:],TRll.r[1,0,:],'g')

    ax = fig.add_subplot(323)
    ax.imshow(CR['tiparav']-mean(CR['tiparav'][0,:]),origin='low',extent=[extent[0],extent[1],extent[2]-extent[3]+10.,10.],alpha=0.0)
    ax.autoscale(False)
    ax.plot(TRll.r[0,0,:],TRll.r[2,0,:],'g')
    #ax.plot(TRlu.r[0,0,:],TRlu.r[2,0,:],'b')
    ax.set_title('x vx z')

    ax = fig.add_subplot(325)
    ax.plot(TRll.r[1,0,:],TRll.r[2,0,:],'g')
    #ax.plot(TRlu.r[1,0,:],TRlu.r[2,0,:],'b')
    ax.set_title('y vx z')

    ax = fig.add_subplot(322)
    timevec = arange(TRll.tstart,TRll.tend+TRll.dt,TRll.dt)
    ax.plot(timevec,TRll.r[0,0,:],'k')
    ax.plot(timevec,TRll.r[1,0,:],'r')
    ax.plot(timevec,TRll.r[2,0,:],'b')
    ax.set_title('x (k), y(r), z(b)')

    ax = fig.add_subplot(324)
    timevec = arange(TRll.tstart,TRll.tend+TRll.dt,TRll.dt)
    ax.plot(timevec,TRll.v[0,0,:],'k')
    ax.plot(timevec,TRll.v[1,0,:],'r')
    ax.plot(timevec,TRll.v[2,0,:],'b')
    ax.set_title('v_x (k), v_y(r), v_z(b)')

    ax = fig.add_subplot(326)
    timevec = arange(TRll.tstart,TRll.tend+TRll.dt,TRll.dt)
    vmag = sqrt(sum(TRll.v[:,0,:]**2,0))
    vpar = abs(sum(TRll.v[:,0,:]*TRll.Bp[:,0,:],0)/sqrt(sum(TRll.Bp[:,0,:]**2,0)))
    vperp = sqrt(vmag**2 - vpar**2)
    ax.plot(timevec,.5*vmag**2,'k')
    ax.plot(timevec,.5*vpar**2,'r')
    ax.plot(timevec,.5*vperp**2,'b')
    ax.set_title('|v|^2 (k), v_par^2(r), v_perp^2(b)')

elif whichflg == 'lu':
    ax= fig.add_subplot(321)
    ax.imshow(CR['tiparav']-mean(CR['tiparav'][0,:]),origin='low',extent=extent,cmap='bwr',vmin=-2.,vmax=2.)
    ax.autoscale(False)
    ax.plot(TRlu.r[0,0,:],TRlu.r[1,0,:],'b')
    #ax.plot(TRll.r[0,0,:],TRll.r[1,0,:],'g')

    ax = fig.add_subplot(323)
    ax.imshow(CR['tiparav']-mean(CR['tiparav'][0,:]),origin='low',extent=[extent[0],extent[1],extent[2]-extent[3]+20.,20.],alpha=0.0)
    ax.autoscale(False)
    #ax.plot(TRll.r[0,0,:],TRll.r[2,0,:],'g')
    ax.plot(TRlu.r[0,0,:],TRlu.r[2,0,:],'b')
    ax.set_title('x vx z')

    ax = fig.add_subplot(325)
    #ax.plot(TRll.r[1,0,:],TRll.r[2,0,:],'g')
    ax.plot(TRlu.r[1,0,:],TRlu.r[2,0,:],'b')
    ax.set_title('y vx z')


    ax = fig.add_subplot(322)
    timevec = arange(TRlu.tstart,TRlu.tend+TRlu.dt,TRlu.dt)
    ax.plot(timevec,TRlu.r[0,0,:],'k')
    ax.plot(timevec,TRlu.r[1,0,:],'r')
    ax.plot(timevec,TRlu.r[2,0,:],'b')
    ax.set_title('x (k), y(r), z(b)')

    ax = fig.add_subplot(324)
    timevec = arange(TRlu.tstart,TRlu.tend+TRlu.dt,TRlu.dt)
    ax.plot(timevec,TRlu.v[0,0,:],'k')
    ax.plot(timevec,TRlu.v[1,0,:],'r')
    ax.plot(timevec,TRlu.v[2,0,:],'b')
    ax.set_title('v_x (k), v_y(r), v_z(b)')

    ax = fig.add_subplot(326)
    timevec = arange(TRlu.tstart,TRlu.tend+TRlu.dt,TRlu.dt)
    vmag = sqrt(sum(TRlu.v[:,0,:]**2,0))
    vpar = abs(sum(TRlu.v[:,0,:]*TRlu.Bp[:,0,:],0)/sqrt(sum(TRlu.Bp[:,0,:]**2,0)))
    vperp = sqrt(vmag**2 - vpar**2)
    ax.plot(timevec,.5*vmag**2,'k')
    ax.plot(timevec,.5*vpar**2,'r')
    ax.plot(timevec,.5*vperp**2,'b')
    ax.set_title('|v|^2 (k), v_par^2(r), v_perp^2(b)')

fig.show()
