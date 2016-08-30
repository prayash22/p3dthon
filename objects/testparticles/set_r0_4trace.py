npart  = 10
charge = 1.
mass   = 1.0
t0     = 5.
tend   = 40.
v0     = 2.3
v0     = zeros(npart*3).reshape(3,npart)
for c,vv in enumerate(linspace(0.,5.,10)):
    v0[:,c] = array([vv,0.,0.])

dr0    = [.0,.0]
dv0    = 0.0*ones(3)
dt     = .025
r0y    = 149.
r0x    = 72.
r0y    = 60.
r0x    = 324.

dx = CR['xx'][1] - CR['xx'][0]
dy = CR['yy'][1] - CR['yy'][0]
ip = int( round( (r0x - CR['xx'][0])/dx ) )
jp = int( round( (r0y - CR['yy'][0])/dy ) )
bb = array([CR['bxav'][jp-.1/dy:jp+.1/dy,ip-1/dx:ip+1/dx],CR['byav'][jp-.1/dy:jp+.1/dy,ip-1/dx:ip+1/dx],CR['bzav'][jp-.1/dy:jp+.1/dy,ip-1/dx:ip+1/dx]])
ee = array([CR['exav'][jp-.1/dy:jp+.1/dy,ip-1/dx:ip+1/dx],CR['eyav'][jp-.1/dy:jp+.1/dy,ip-1/dx:ip+1/dx],CR['ezav'][jp-.1/dy:jp+.1/dy,ip-1/dx:ip+1/dx]])
exb = cross(ee,bb,axis=0)/sum(bb**2,axis=0)
bbh = bb/sqrt(sum(bb**2,axis=0))

#sys.exit()

exb = mean(exb,axis=(1,2))
bbh = mean(bbh,axis=(1,2))

r0 = [r0x,r0y]

TR = TPRun(CR,
           npart=npart,
           charge=charge, 
           mass=mass, 
           tstart=0., 
           tend=tend,
           t0=t0,
           r0=r0, 
           dr0=dr0, 
           v0= v0,
           #v0= 0.*bbh,
           dv0=dv0,
           dt=dt)

TR.move()

