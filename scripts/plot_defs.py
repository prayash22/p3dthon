from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl
import matplotlib.pyplot as plt
from collections import OrderedDict
from sub import *

#==============================================

def set_fig():
    mpl.rcParams.update({'figure.autolayout': True})
    fig = plt.figure(1,figsize=(8.5, 11.))
    fig.clf()
    axs = [fig.add_subplot(6,2,x+1) for x in range(12)]
    return fig,axs

#==============================================

def close_fig(pdf,nosave=False):
    if nosave:
        plt.show()
        plt.draw()
        #plt.tight_layout()
        raw_input('nosave flag set. Return to see next figure: ')
        plt.clf()
    else:
        #plt.tight_layout()
        pdf.savefig()  # saves the current figure into a pdf page
        plt.close()

#==============================================

def plot_2D(ax,CR,var,
            extent='',
            psi_lvls=None,
            cut_locs=None,
            **kwargs):

    ext = [CR['xx'][0],
           CR['xx'][-1],
           CR['yy'][0],
           CR['yy'][-1]]

    if type(var) is str: 
        title = var
        var = CR[var]
    else:
        if 'title' in kwargs:
            title = kwargs.pop('title')
        else:
            title = 'NO TITLE SET!!!'
    
    im = ax.imshow(var,extent=ext,origin='low',**kwargs)
    
    if type(extent) is not str:
        ax.set_xlim(extent[:2])
        ax.set_ylim(extent[2:])
    ax.autoscale(False)
    
    ax.set_xlabel(r'$X\ (d_i)$',size=8)
    ax.set_ylabel(r'$Y\ (d_i)$',size=8)
    ax.set_title(title+': %1.3f, %1.3f'%(var.min(),var.max()),size=8)
    
    lsz = 6
    ax.xaxis.set_tick_params(which='both',labelsize=lsz)
    ax.yaxis.set_tick_params(which='both',labelsize=lsz)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", "3%", pad="1.5%")
    plt.colorbar(im, cax=cax)

    lsz = 6
    cax.xaxis.set_tick_params(which='both',labelsize=lsz)
    cax.yaxis.set_tick_params(which='both',labelsize=lsz)
    plt.sca(ax)
    plt.minorticks_on()

    if psi_lvls is not None:

        cs = ax.contour(CR['xx'],
                        CR['yy'],
                        CR['psi'],
                        levels = psi_lvls,
                        colors='k',
                        linewidths=.5)

        [c.set_linestyle('solid') for c in cs.collections]

    if cut_locs is not None:
        for cosa in cut_locs:
            draw_line(ax,'y',cosa)

    return im,cax

#==============================================

def plot_1D(ax,CR,var,
            dir='y',
            loc=0.,
            extent=None,
            **kwargs):

    lgargs = kwargs.pop('lgargs',None)
    no_lines = kwargs.pop('no_lines',False)
    c_itter = kwargs.pop('c_itter',['r','b','g','k'])

    if dir == 'y':

        ip = abs(CR['xx'] - loc).argmin()
        
        if type(var) in [dict,OrderedDict]:
            ln = []
            for key in var:
                lab ='$'+key+'$'
                ln.append(\
                ax.plot(CR['yy'],
                        var[key][:,ip],
                        color=c_itter.pop(0),
                        label=r'%s'%lab,
                        **kwargs)[0])

            ax.set_xlabel(r'$Y (d_i)$',fontsize=8)

            if lgargs is not None:
                lg = ax.legend(**lgargs)
            else:
                lg = ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), 
                               loc=3,
                               ncol=len(var),
                               borderaxespad=0.,
                               prop={'size':6})

        else:
            ln = ax.plot(CR['yy'],var[:,ip],**kwargs)

    if extent is not None:
        ax.set_xlim(extent[:2])
        if np.size(extent) > 2:
            ax.set_ylim(extent[2:])
    else:
        if dir == 'y':
            ax.set_xlim(CR['yy'][[0,-1]])
        else:
            ax.set_xlim(CR['xx'][[0,-1]])
    lsz = 6
    ax.xaxis.set_tick_params(which='both',labelsize=lsz)
    ax.yaxis.set_tick_params(which='both',labelsize=lsz)
    ax.set_title('cut @ x = %1.2f'%loc,size=8,loc='right')

    ax.autoscale(False)

    vcuts = locate_seps(CR,ip)
    
    if not no_lines:
        for offset in vcuts:
            draw_line(ax,cut='y',offset=offset)
    else:
# Draw line at mind plane
        xmp = CR['yy'][abs(CR['bxav'][:,ip]).argmin()]
        draw_line(ax,cut='x',offset=0.)
        draw_line(ax,cut='y',offset=xmp)
# or just draw it at the new 0
        #draw_line(ax,cut='y',offset=0.)

    plt.sca(ax)
    plt.minorticks_on()

    return ln,lg
 

#==============================================

def draw_line(ax,cut='x',offset=0.):
    ax.autoscale(False)
    bgarr = np.array([-10000,10000])

    if cut == 'x':
        ax.plot(bgarr, np.zeros(2)+offset, 'k--')

    elif cut == 'y':
        ax.plot(np.zeros(2)+offset, bgarr, 'k--')

    else: print 'What?'
        


#==============================================

def calc_extra_vars(CR):
    print 'Calculating Extra Variables to be included...'

    if 'ni' not in CR:
        CR['ni'] = CR['deni']
    if 'ne' not in CR:
        CR['ne'] = CR['dene']

    CR['psi'] = calc_psi(CR)
    CR['rho'] = CR['ni'] - CR['ne']
    CR['|b|'] = np.sqrt(CR['bx']**2+CR['by']**2+CR['bz']**2)
    for v in ['x','y','z']:
        for q,s in zip([1.,-1.],['i','e']):
            CR['v'+s+v] = q*CR['j'+s+v]/CR['n'+s]
    for v in ['xx','yy','zz','xy','xz','yz']:
        for s in ['i','e']:
            CR['t'+s+v] = CR['p'+s+v]/CR['n'+s]

    rotate_ten(CR,'ti',av='')
    rotate_ten(CR,'te',av='')

    bclip = CR['|b|'].clip(.05)

    CR['exbx'] = (CR['ey']*CR['bz'] - CR['ez']*CR['by'])/bclip**2
    CR['exby'] = (CR['ez']*CR['bx'] - CR['ex']*CR['bz'])/bclip**2
    CR['exbz'] = (CR['ex']*CR['by'] - CR['ey']*CR['bx'])/bclip**2
    CR['|exb|'] = np.sqrt(CR['exbx']**2+CR['exby']**2+CR['exbz']**2)

    CR['epar'] = (CR['ex']*CR['bx'] +
                  CR['ey']*CR['by'] +
                  CR['ez']*CR['bz'])/CR['|b|']
    CR['jpar'] = (CR['jx']*CR['bx'] +
                  CR['jy']*CR['by'] +
                  CR['jz']*CR['bz'])/CR['|b|']
    CR['phin'] = CR['tepar']*np.log(CR['ne']/np.mean(CR['ne'][[0,-1],:]))
    CR['phit'] = CR['tepar'] - np.mean(CR['tepar'][[0,-1],:])


    CR['vipar'] = (CR['vix']*CR['bx'] +
                   CR['viy']*CR['by'] +
                   CR['viz']*CR['bz'])/CR['|b|']
    CR['vepar'] = (CR['vex']*CR['bx'] +
                   CR['vey']*CR['by'] +
                   CR['vez']*CR['bz'])/CR['|b|']

    CR['viperp'] = (CR['vix']*CR['exbx'] +
                    CR['viy']*CR['exby'] +
                    CR['viz']*CR['exbz'])/CR['|exb|']
    CR['veperp'] = (CR['vex']*CR['exbx'] +
                    CR['vey']*CR['exby'] +
                    CR['vez']*CR['exbz'])/CR['|exb|']

    CR['pb']   = CR['|b|']**2/2.
    CR['ptot'] = CR['pb'] + CR['piyy'] + CR['peyy']


#==============================================

def cut_n_cont_locs(CR):
# Frist we need to find the xline
    jpx = np.shape(CR['psi'])[0]/2
    if CR['yy'][0] < 1.0:                 # If lower half
        ipx = CR['psi'][jpx,:].argmin()
    else:                           # Else upper half
        ipx = CR['psi'][jpx,:].argmax()

# Now we are going to determin the cut locations in real
# space and let the plot method figure the rest out
    xp = CR['xx'][ipx]
    yp = CR['yy'][jpx]

    lx = CR['xx'][-1] + CR['xx'][0] 
    ly = CR['yy'][-1] + CR['yy'][1] - 2.*CR['yy'][0]

# We will have a cut a 2.5 and 
# then double the lenght each time
# But this has now been replaced with
# the while loop for a better spread
# of cuts
    #ncuts = int(np.log(lx/2./2.5)/np.log(2.)) 
    #cut_locs = [2.5*2.**abs(n) for n in range(ncuts+1)]
    
    _dx = 5.0
    _ct = 0 
    _mx = 2
    cut_locs = [0.0]
    while cut_locs[-1] < lx/2.0:
        cut_locs += [(cut_locs[-1] + _dx)]
        _ct += 1
        if _ct >= _mx:
            _ct = 0
            _mx += 2
            _dx *= 2.

    cut_locs[-1] = lx/2.0
    cut_locs = [-1.0*x for x in cut_locs[::-1]] + \
               [0.] + cut_locs + [lx/2.]

    cut_locs = [(xp + x)%lx for x in cut_locs]
    cut_locs.sort()

    psi0 = CR['psi'][-1,ipx]
    dpsi = (CR['psi'][jpx,ipx] - psi0)/5.
    psi_lvls = np.arange(psi0,10*dpsi + psi0,dpsi)

    return cut_locs,psi_lvls


def locate_seps(CR,ip):
    """ This Function take
        CR: a standered simulation and
        ip: the x index for the y cut location

        and returns

        ypl: the real location of the left 
             speratrix
        ypm: the real location of the 
             midplane
        ypr: the real location of the right 
             speratrix
     """
# Frist we need to find the xline
    jpx = np.shape(CR['psi'])[0]/2
    if CR['yy'][0] < 1.0:                 # If lower half
        ipx = CR['psi'][jpx,:].argmin()
    else:                           # Else upper half
        ipx = CR['psi'][jpx,:].argmax()
    
    psi0 = CR['psi'][jpx,ipx]
    ypl = abs(CR['psi'][:jpx,ip] - psi0).argmin()
    ypr = abs(CR['psi'][jpx:,ip] - psi0).argmin()

    ypl = CR['yy'][ypl]
    ypm = CR['yy'][jpx]
    ypr = CR['yy'][jpx + ypr]

    return [ypl,ypm,ypr]
