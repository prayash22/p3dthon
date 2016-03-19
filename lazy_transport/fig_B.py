import matplotlib.gridspec as gridspec
from plot_defs import *
from collections import OrderedDict
from paper_fig_defs import *

plt.ioff()

# We need to build the figure
fig = figure(1); fig.clf()
fig.set_size_inches(8.0,10.5)

#============================================
# First we make the 1D plots on the left hand side

gs1 = gridspec.GridSpec(9, 1)
gs1.update(left=0.05, right=0.48, bottom=.035, top=.95, hspace=0.15)
axL = [plt.subplot(gs1[c, :]) for c in range(9)]


load_data = (raw_input('Load data? ') is 'y')
if  load_data:
    xp = 26.275
    yp = 37.875
    CR = readsave('./tave100reconn414u.dat')
    CR['xx'] = CR['xx'] - xp
    CR['yy'] = CR['yy'] - yp
    calc_vars(CR)

# We will comeback and do the distro later
vars_1D = [OrderedDict([(k , CR[k]) for k in ('bxav', '-bzav', 'byav')]),
           OrderedDict([(k , CR[k]) for k in ('eyav', '-vixby')]),
           OrderedDict([(k , CR[k]) for k in ('deniav', 'deneav')]),
           OrderedDict([(k , CR[k]) for k in ('teparav', 'teperp1av')]),
           OrderedDict([(k , CR[k]) for k in ('vexav', '-vezav', 'veyav')]),
           OrderedDict([(k , CR[k]) for k in ('vixav', '-vizav', 'viyav')]),
           OrderedDict([(k , CR[k]) for k in ('-vexbyav', 'vdgvyav',
                                              'gpynav', 'eyav')]),
           OrderedDict([(k , CR[k]) for k in ('piyyav', 'peyyav', 'pb', 'ptot')]),
           OrderedDict([(k , CR[k]) for k in ('texxav', 'tezzav', 'teyyav')])]


labels = {'bxav'     :'$B_{L}$',
          'byav'     :'$B_{N}$',
          '-bzav'    :'$B_{M}$',
          'eyav'     :'$E_{N}$',
          '-vixby'   :'$-(v_{i}\\times B)_{N}$',
          'deniav'   :'$n_{i}$',
          'deneav'   :'$n_{e}$',
          'teparav'  :'$T_{e\parallel}$',
          'teperp1av':'$T_{e\perp}$',
          'vexav'    :'$v_{eL}$',
          'veyav'    :'$v_{eN}$',
          '-vezav'   :'$v_{eM}$',
          'vixav'    :'$v_{iL}$',
          'viyav'    :'$v_{iN}$',
          '-vizav'   :'$v_{eM}$',
          '-vexbyav'  :'$-(v_{e}\\times B)_{N}$',
          'vdgvyav'  :'$-m_e v_e \\cdot \\nabla v_{eN}$',
          'gpynav'   :'$-(\\nabla P_e)_N/n_e $',
          'eyav'     :'$E_N$',
          'piyyav'   :'$P_{iNN}$',
          'peyyav'   :'$P_{eNN}$',
          'pb'       :'$P_{B}$',
          'ptot'     :'$P_{tot}$',
          'texxav'   :'$P_{iLL}$',
          'teyyav'   :'$P_{iNN}$',
          'tezzav'   :'$P_{iMM}$'}
           
xp = 26.275
yp = 37.875

lns = []
lgs = []
lgargs = {'loc':'best', 
          'ncol':1,
          'frameon':False,
          'prop':{'size':6}}

for c,(ax,v) in enumerate(zip(axL,vars_1D)):
    if c in [2, 4]: 
        lgargs['loc'] = 'lower left' 
    else: lgargs['loc'] = 'upper left' 

    if c == 6: lgargs['ncol'] = 2
    elif c == 7: lgargs['ncol'] = 4
    else: lgargs['ncol'] = 1

    _ = plot_1D(ax,CR,v,
                dir='y',
                loc=0.,
                #label=labels[c],
                lgargs=lgargs,
                no_lines=True,
                extent=(-6.,+6.))

    lns.append(_[0])
    lgs.append(_[1])

# Make the text bigger
[lg.texts[0].set_size(10) for lg in lgs]

# Cleaning up
[ax.invert_xaxis() for ax in axL]
[ax.set_title('',loc='right') for ax in axL]
[plt.setp(ax.get_xticklabels(), visible=False) for ax in axL[:-1]]
for lg in lgs:
    for txt in lg.texts:
        txt.set_text(labels[txt.get_text()[1:-1]])

for ax in axL[:-1]:
    ax.set_xlabel('') 
axL[6].set_ylim(-4.,8.)
axL[2].set_ylim(.05,1.5)
axL[2].set_yscale('log')


#============================================
# Now we handel the distros

gs2 = gridspec.GridSpec(9, 1)
gs2.update(left=0.54, right=0.97, bottom=.035, top=.9575, hspace=0.15)
axR = []
axR.append(plt.subplot(gs2[0, :]))

#===================
# Do the 2d box plot

if load_data:
    _CR = readsave('./tave100reconn414u.dat')
    roll_run(_CR)
    _CR['psi'] = calc_psi(_CR)
    calc_EvxB(_CR)

xp = 0.
yp = 0.

ip = int(round((xp - _CR['xx'][0])/(_CR['xx'][1] - _CR['xx'][0])))
jp = int(round((yp - _CR['yy'][0])/(_CR['yy'][1] - _CR['yy'][0])))

psi_levels = [_CR['psi'][jp,ip] - 2.5*z for z in range(-4,4)]

extent = [-30.,30.,yp-15./2.,yp+15./2.]
txtloc = [xp-1., yp+3.]

mp = get_mpline(_CR)
mip = ((mp - _CR['yy'][0])/\
       (_CR['xx'][1] - _CR['xx'][0])\
       ).round().astype(int).tolist()

k = [_CR['eyav'],  'seismic', '$E_N$']
im,cm = plot_2D(axR[0], _CR, k[0],
                extent=extent,
                cmap=k[1],
                psi_lvls=psi_levels,
                cut_locs=None)

axR[0].plot(_CR['xx'],mp,'k--')

ft = axR[0].text(txtloc[0],
            txtloc[1],
            r'%s'%k[2],
            fontsize=12,
            weight='bold',
            horizontalalignment='center',
            verticalalignment='center')

# Set title and stuff
im.set_clim(-1.*max(abs(array(im.get_clim()))),
                max(abs(array(im.get_clim()))))
axR[0].set_title('')
axR[0].set_xlabel(r'$L\ (d_i)$')
axR[0].set_ylabel(r'$N\ (d_i)$')
axR[0].yaxis.labelpad=-5.
#Code to draw boxes

#===================
# Do the distros

gs3 = gridspec.GridSpec(4, 2)
gs3.update(left=0.54, right=0.97, bottom=.035, top=.825, hspace=0.15)
axR.append(plt.subplot(gs3[0, 0]))
axR.append(plt.subplot(gs3[0, 1]))
axR.append(plt.subplot(gs3[1, 0]))
axR.append(plt.subplot(gs3[1, 1]))
axR.append(plt.subplot(gs3[2, 1]))
axR.append(plt.subplot(gs3[3, 1]))


### Lets handel the distro stuff later
plt.draw()
plt.show()
#fig.savefig('fig_B.eps')
sys.exit()

if raw_input('Load data? ') == 'y':
    CC=p3d_run('ar414')
    CR = CC.load_movie('all',99,movie_num=001)

# setting up some basic distro values
# determin the range and box size
bins = 75 
rgi = [[-7.5,7.5],[-7.5,7.5]]
rge = [[-15,15],[-15,15]]

# set the pizza origin and width
z0 = 0.  
dz = 6. # note this is for e, i -> dz/5

###################################### 1

x0   = 26.275
y0r  = [37.125, 37.125 +.5,37.125 + 1.]

dx = array([1.,.5])

fig = []
axs = []
ima = []

plt.ioff() # Trun off interactive mode

count = 1
for y0 in y0r:
    for LMN in [True,False]:
    
        f1,ax = set_f1(count)
        fig.append(f1)
        axs.append(ax)

        ima.append(\
        plot_pizza(CC,
                   ax=ax,
                   r0=[x0,y0],
                   dx=dx,
                   range=rge,
                   bins=bins,
                   dump_num=002,
                   par=not LMN,
                   pizza=True,
                   z0=z0,
                   dz=dz,
                   wax=0,
                   LMN=LMN,
                   sp ='e'))
        
        fname = 'fig_B_%i.eps'%count
        print 'Saving ',fname
        #f1.savefig(fname)
        count += 1

imax = max([c.get_clim()[1] for c in ima])
[c.set_clim(.1,imax) for c in ima]
print '\nMax = %1.1f\n'%imax

for c,f1 in enumerate(fig):
    fname = 'fig_B_%i.eps'%c
    print 'Saving ',fname
    f1.savefig(fname)



