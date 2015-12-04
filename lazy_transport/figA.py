from plot_defs import *
from paper_fig_defs import *

if raw_input('Load data? ') == 'y':
    CR = readsave('./tave100reconn414u.dat')
    roll_run(CR)
    CR['psi'] = calc_psi(CR)
    calc_EvxB(CR)

plt.ioff()

fig = figure(1);fig.clf()
fig.set_size_inches(8.,10.5)
ax = [fig.add_subplot(9,2,c+1) for c in range(9*2)]
fig.subplots_adjust(left  =.08,
                    right =.95,
                    bottom=.15,
                    top   =.99)

xp = 26.275 
yp = 37.875
xp = 0.0
yp = 0.0

ip = int(round((xp - CR['xx'][0])/(CR['xx'][1] - CR['xx'][0])))
jp = int(round((yp - CR['yy'][0])/(CR['yy'][1] - CR['yy'][0])))

psi_levels = [CR['psi'][jp,ip] - 2.5*z for z in range(-4,4)]

extent = [-30.,30.,yp-15./2.,yp+15./2.]
txtloc = [xp-1., yp+3.]

mp = get_mpline(CR)
mip = ((mp - CR['yy'][0])/\
       (CR['xx'][1] - CR['xx'][0])\
       ).round().astype(int).tolist()

#        Varaible        cmap       title
var = [[CR['vixav'],  'seismic', '$v_{iL}$'],
[-1.*CR['bzav'],      'seismic', '$B_{M}$'],
    [CR['vexav'],     'seismic', '$v_{eL}$'],
    [CR['eyav'],      'seismic', '$E_{N}$'],
[-1.*CR['vezav'],     'seismic', '$v_{iM}$'],
    [CR['EviBy'],     'seismic', '$(E + v_i\\times B)_{N}$'],
    [CR['teparav'],   'afmhot' , '$T_{e\parallel}$'],
    [CR['EveBy'],     'seismic', '$(E+ v_e\\times B)_{N}$'],
    [CR['teperp1av'], 'afmhot' , '$T_{e\perp}$'],
    [CR['De'],        'seismic', '$D_{e}$'],
[-1.*CR['jzav'],      'seismic', '$J_{M}$'],
    [CR['Ae'],        'afmhot' , '$A_{e}$'],
    [CR['deniav'],    'afmhot' , '$n_{i}$'],
    [CR['bx'],        'seismic', '$B_{L}$'],
    [CR['jxav'],      'seismic', '$J_{L}$'],
[-1.*CR['EviBz'],     'seismic', '$(E + v_i\\times B)_{M}$'],
[-1.*CR['vezav'],     'seismic', '$v_{eM}$'],
[-1.*CR['EveBz'],     'seismic', '$(E + v_e\\times B)_{M}$']]
      
ima = []
cma = []
for c,k in enumerate(var):
    
    _ = plot_2D(ax[c], CR, k[0],
                extent=extent,
                cmap=k[1],
                psi_lvls=psi_levels,
                cut_locs=None)

    ima.append(_[0])
    cma.append(_[1])

    ax[c].plot(CR['xx'],mp,'k--')

    ft = ax[c].text(txtloc[0],
               txtloc[1],
               r'%s'%k[2],
               fontsize=12,
               weight='bold',
               horizontalalignment='center',
               verticalalignment='center')
    if k[2] is '$T_{e\parallel}$' or \
       k[2] is '$T_{e\perp}$':
        ft.set_color('white')


# Now we clean up a bit
for c,a in enumerate(ax):
    a.set_title('')
    a.set_xlabel('')
    if c%2 == 0:
        a.set_ylabel(r'$N\ (d_{i})$')
    else:
        a.set_ylabel('')
    if c < 16:
        a.set_xlabel('')
    else:
        a.set_xlabel(r'$L\ (d_{i})$')

    if ima[c].get_cmap().name is 'seismic':
        cl = ima[c].get_clim()
        if abs(cl[0]) > abs(cl[1]):
            ima[c].set_clim(cl[0],abs(cl[0]))
        else:
            ima[c].set_clim(-1.*cl[1],cl[1])


plt.show()
plt.draw()

if raw_input('Save fig? ') is 'y':
    fname = 'figA.eps'
    print 'Saving ',fname
    fig.savefig(fname,rasterized=True,dpi=100)
