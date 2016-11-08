import time
import operator
import numpy as np
import matplotlib.pyplot as plt
from scipy.io.idl import readsav 
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import AutoMinorLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable
from p3d_runs import p3d_run

#======================================================
def set_local(IDL_restore,lcl):
    #lcl = locals()
    print 'Setting B ...'
    if ('bxav' in IDL_restore) and ('byav' in IDL_restore) and ('bzav' in IDL_restore): lcl['B'] = np.array([IDL_restore['bxav'],IDL_restore['byav'],IDL_restore['bzav']])
    elif 'bx' in IDL_restore and 'by' in IDL_restore and 'bz' in IDL_restore: lcl['B'] = np.array([IDL_restore['bx'],IDL_restore['by'],IDL_restore['bz']])
    print 'Setting E ...'
    if 'exav' in IDL_restore and 'eyav' in IDL_restore and 'ezav' in IDL_restore: lcl['E'] = np.array([IDL_restore['exav'],IDL_restore['eyav'],IDL_restore['ezav']])
    elif 'ex' in IDL_restore and 'ey' in IDL_restore and 'ez' in IDL_restore: lcl['E'] = np.array([IDL_restore['ex'],IDL_restore['ey'],IDL_restore['ez']])
    print 'Setting J ...'
    if 'jxav' in IDL_restore and 'jyav' in IDL_restore and 'jzav' in IDL_restore: lcl['J'] = np.array([IDL_restore['jxav'],IDL_restore['jyav'],IDL_restore['jzav']])
    elif 'jx' in IDL_restore and 'jy' in IDL_restore and 'jz' in IDL_restore: lcl['J'] = np.array([IDL_restore['jx'],IDL_restore['jy'],IDL_restore['jz']])
    print 'Setting Ji ...'
    if 'jixav' in IDL_restore and 'jiyav' in IDL_restore and 'jizav' in IDL_restore: 
        lcl['Ji'] = np.array([IDL_restore['jixav'],IDL_restore['jiyav'],IDL_restore['jizav']])
        if 'J' in lcl and 'Ji' in lcl: lcl['Je'] = lcl['J'] - lcl['Ji']
    elif 'jix' in IDL_restore and 'jiy' in IDL_restore and 'jiz' in IDL_restore: 
        lcl['J'] = np.array([IDL_restore['jix'],IDL_restore['jiy'],IDL_restore['jiz']])
        if 'J' in lcl and 'Ji' in lcl: lcl['Je'] = lcl['J'] - lcl['Ji']
    print 'Setting deni ...'
    if 'deniav' in IDL_restore: lcl['deni'] = IDL_restore['deniav']
    elif 'deni' in IDL_restore: lcl['deni'] = IDL_restore['deni']
    print 'Setting dene ...'
    if 'deneav' in IDL_restore: lcl['dene'] = IDL_restore['deneav']
    elif 'dene' in IDL_restore: lcl['dene'] = IDL_restore['dene']
    print 'Setting P ...'
    if TestTen('Pi',lcl,'av'): lcl['Pi'] = np.array([[IDL_restore['pixxav'],IDL_restore['pixyav'],IDL_restore['pixzav']],
                                                  [IDL_restore['pixyav'],IDL_restore['piyyav'],IDL_restore['piyzav']],
                                                  [IDL_restore['pixzav'],IDL_restore['piyzav'],IDL_restore['pizzav']]])
    elif TestTen('Pi',lcl): lcl['Pi'] = np.array([[IDL_restore['pixx'],IDL_restore['pixy'],IDL_restore['pixz']],
                                               [IDL_restore['pixy'],IDL_restore['piyy'],IDL_restore['piyz']],
                                               [IDL_restore['pixz'],IDL_restore['piyz'],IDL_restore['pizz']]])
    if TestTen('Pe',lcl,'av'): lcl['Pe'] = np.array([[IDL_restore['pexxav'],IDL_restore['pexyav'],IDL_restore['pexzav']],
                                                  [IDL_restore['pexyav'],IDL_restore['peyyav'],IDL_restore['peyzav']],
                                                  [IDL_restore['pexzav'],IDL_restore['peyzav'],IDL_restore['pezzav']]])
    elif TestTen('Pe',lcl): lcl['Pe'] = np.array([[IDL_restore['pexx'],IDL_restore['pexy'],IDL_restore['pexz']],
                                               [IDL_restore['pexy'],IDL_restore['peyy'],IDL_restore['peyz']],
                                               [IDL_restore['pexz'],IDL_restore['peyz'],IDL_restore['pezz']]])


#======================================================
def TestTen(var,lcl,av=''):
    if var+'xx'+av in lcl and var+'yy'+av in lcl and var+'zz'+av in lcl and var+'xy'+av in lcl and var+'yz'+av in lcl and var+'yz'+av in lcl:
        return True
    else: 
        return False

#======================================================
def rotate_ten(CR,
               var='pi',
               av='av',
               overwrite=False,
               full_rotate=False):

    if var+'par'+av in CR and not overwrite:
        print 'Warning: %sparav was found in the' \
              'restored data: nothing will be rotated!!!!'
        pass

        
    elif full_rotate:
# e1 -> \hat{B} 
# e2 -> \hat{ExB}
# e3 -> \hat{Bx(ExB)}
        e1 = np.array([CR['bxav'],
                       CR['byav'],
                       CR['bzav']])

        e2 = np.cross(np.array([CR['exav'],
                                CR['eyav'],
                                CR['ezav']]), e1,axis=0)

        e1 = e1/np.sqrt(np.sum(e1**2,axis=0))
        e2 = e2/np.sqrt(np.sum(e2**2,axis=0))
        e3 = np.cross(e1,e2,axis=0)

        T = np.array([[CR[var+'xx'+av],CR[var+'xy'+av],CR[var+'xz'+av]],
                      [CR[var+'xy'+av],CR[var+'yy'+av],CR[var+'yz'+av]],
                      [CR[var+'xz'+av],CR[var+'yz'+av],CR[var+'zz'+av]]])

        Te1 = np.array([np.sum(T[0,:,:,:]*e1,axis=0),
                        np.sum(T[1,:,:,:]*e1,axis=0),
                        np.sum(T[2,:,:,:]*e1,axis=0)])
        Te2 = np.array([np.sum(T[0,:,:,:]*e2,axis=0),
                        np.sum(T[1,:,:,:]*e2,axis=0),
                        np.sum(T[2,:,:,:]*e2,axis=0)])
        Te3 = np.array([np.sum(T[0,:,:,:]*e3,axis=0),
                        np.sum(T[1,:,:,:]*e3,axis=0),
                        np.sum(T[2,:,:,:]*e3,axis=0)])

        Tpar  = np.sum(e1*Te1,axis=0)
        Tperp = (T[0,0,:,:] + T[1,1,:,:] + T[2,2,:,:] - Tpar)/2.0

# The diaganal is easy, we pick perp1 = perp2
        CR[var+'par'+av]   = Tpar
        CR[var+'perp1'+av] = Tperp
        CR[var+'perp2'+av] = Tperp

        a = np.sum(e2*Te1,axis=0)
        b = np.sum(e3*Te1,axis=0)
        c = np.sum(e2*Te2,axis=0)
        d = np.sum(e3*Te2,axis=0)
        e = np.sum(e3*Te3,axis=0)
        
        x = (c - e)/d/2.
        ct = np.sqrt(1. + 1./np.sqrt((1.+ (x)**2)))/np.sqrt(2)
        st = x/(np.sqrt(2.)*np.sqrt(x**2.+1.)*\
                np.sqrt(1./np.sqrt(x**2.+1.)+1.))

        T11 = Tpar
        T12 = a*ct - b*st
        T13 = a*st + b*ct
        T22 = c*ct**2 + e*st**2 - 2.*d*st*ct
        T23 = (c - e)*st*ct + d*(ct**2 - st**2)
        T33 = c*st**2 + e*ct**2 + 2.*d*st*ct

        CR[var+'11'+av] = T11
        CR[var+'12'+av] = T12
        CR[var+'13'+av] = T13
        CR[var+'22'+av] = T22
        CR[var+'23'+av] = T23
        CR[var+'33'+av] = T23

        # Now Agyrotropy code!
        CR[var+'agy'+av] = 2.*np.sqrt((T12**2 + T13**2 + T23**2)) \
                                 /(T11 + T22 + T33) 
    else:
# This was the old way, and it was very simple

        bmag = np.sqrt( CR['bx'+av]**2+
                        CR['by'+av]**2+
                        CR['bz'+av]**2)
        bbx = CR['bx'+av]/bmag
        bby = CR['by'+av]/bmag
        bbz = CR['bz'+av]/bmag

        CR[var+'par'+av] = (bbx*(bbx*CR[var+'xx'+av] + 
                                 bby*CR[var+'xy'+av] + 
                                 bbz*CR[var+'xz'+av])+
                            bby*(bbx*CR[var+'xy'+av] + 
                                 bby*CR[var+'yy'+av] + 
                                 bbz*CR[var+'yz'+av])+
                            bbz*(bbx*CR[var+'xz'+av] + 
                                 bby*CR[var+'yz'+av] + 
                                 bbz*CR[var+'zz'+av]))

        CR[var+'perp1'+av] = (CR[var+'xx'+av] + 
                              CR[var+'yy'+av] +
                              CR[var+'zz'+av] - 
                              CR[var+'par'+av])/2.

        CR[var+'perp2'+av] = CR[var+'perp1'+av]

#======================================================

def ims(fdic,
        key,
        ax=None,
        extent=None,
        cbar=None,
        cont=None,
        no_draw=None,
        ctargs={},
        **kwargs):
    """
    A wrapper function for imshow to do most 
    tedious stuff for my simulations
    """
    old_ax = plt.gca() # Get Current Axis
    if ax is None: 
        ax = old_ax
    else:
        plt.sca(ax)    # Set Current Axis

# Use the dict values of xx and yy to set extent
    ext = [fdic['xx'][0],
           fdic['xx'][-1],
           fdic['yy'][0],
           fdic['yy'][-1]]

    if type(key) is str: plt_val = fdic[key]
    else               : plt_val = key

    if kwargs.has_key('cmap'): cmap=kwargs.pop('cmap')
    else:                      cmap='PuOr'

    im = ax.imshow(plt_val,
                           origin='low',
                           extent=ext,
                           cmap=cmap,            # I just love this color map
                           **kwargs)

    if extent is not None:
        ax.set_xlim(extent[:2])
        ax.set_ylim(extent[2:])

    ax.autoscale(False)

    ax.set_xlabel(r'$X (d_i)$',size=8)
    ax.set_ylabel(r'$Y (d_i)$',size=8)

    ax.xaxis.set_tick_params(which='both',labelsize=8)
    #minorLocator = AutoMinorLocator()           # Note the second call is so that the minor x ticks are not
    #ax.xaxis.set_minor_locator(minorLocator)    # the same as the y ticks

    ax.yaxis.set_tick_params(which='both',labelsize=8)
    #minorLocator = AutoMinorLocator()
    #ax.yaxis.set_minor_locator(minorLocator)

    plt.minorticks_on()

    return_tup = im,
    # Code to implement for a cbar
    if cbar:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", "3%", pad="1.5%")
        plt.colorbar(im, cax=cax)

        cax.xaxis.set_tick_params(which='both',labelsize=8)
        cax.yaxis.set_tick_params(which='both',labelsize=8)

        return_tup += (cax,)

    if cont:
        if 'psi' in fdic:
            psi = fdic['psi']
        else:
            psi = calc_psi(fdic)
        if 'colors' not in ctargs: ctargs['colors'] = 'k'
        if 'linestyles' not in ctargs: ctargs['linestyles'] = 'solid'
        
        cts = ax.contour(fdic['xx'],fdic['yy'],psi,**ctargs)

        return_tup += (cts,)

    if not no_draw:
        plt.draw()
    plt.sca(old_ax)

    if len(return_tup) == 1:
        return return_tup[0]
    else:
        return return_tup

#======================================================

def find_xpt(d):
    psi = calc_psi(d)
    jp = int(np.round(len(d['yy'])/2.))
    yp = d['yy'][jp]

    if 'bxav' in d: av = 'av'
    else: av = ''
    lBm = np.mean(d['bx'+av][:jp,:])
    uBm = np.mean(d['bx'+av][jp:,:])

    if lBm > uBm: #Upper
        ip = psi[jp,:].argmax()
        print 'Finding upper max'
    else:         #Lower
        ip = psi[jp,:].argmin()
        print 'Finding lower min'
    xp = d['xx'][ip]

    return ip,jp,xp,yp

#======================================================
    
def ims_subplot(d,var,ax,window,**kwargs):
    ix = lambda x,xx: np.abs(xx - x).argmin()
    iwd = [ix(w,d[v]) for w,v in zip(window,['xx','yy','xx','yy'])]
    
    if 'bxav' in d: av = 'av'
    else: av = ''
    td = {'xx':d['xx'][iwd[0]:iwd[2]],
          'yy':d['yy'][iwd[1]:iwd[3]],
          'bx'+av:d['bx'+av][iwd[1]:iwd[3], iwd[0]:iwd[2]],
          'by'+av:d['by'+av][iwd[1]:iwd[3], iwd[0]:iwd[2]],
          'bz'+av:d['bz'+av][iwd[1]:iwd[3], iwd[0]:iwd[2]]}

    rst = ims(td, var[iwd[1]:iwd[3], iwd[0]:iwd[2]], ax, **kwargs)

    return rst 

#======================================================

def shift_to_xpt_frame(d):
    ip,jp,xp,yp = find_xpt(d)
    d['xx'] = d['xx'] - xp  
    d['yy'] = d['yy'] - yp  

#======================================================

def var_at(fdic,key,r0,ordflg='idl'):
    delx = fdic['xx'][1] - fdic['xx'][0] 
    dely = fdic['yy'][1] - fdic['yy'][0] 

    xind = (np.floor((r0[0]-delx/2.0)/delx)).astype(int)
    yind = (np.floor((r0[1] - fdic['yy'][0])/dely)).astype(int)

    wx = (r0[0]-delx/2.0)%delx
    wy = (r0[1]-fdic['yy'][0])%dely

    if ordflg =='idl':
        var = wx     *wy     *fdic[key][yind    ,xind    ] + \
              (1.-wx)*wy     *fdic[key][yind    ,(xind+1)] + \
              wx     *(1.-wy)*fdic[key][(yind+1),xind    ] + \
              (1.-wx)*(1.-wy)*fdic[key][(yind+1),(xind+1)] 
    else:
        var = wx     *wy     *fdic[key][xind    ,yind] + \
              (1.-wx)*wy     *fdic[key][(xind+1),yind] + \
              wx     *(1.-wy)*fdic[key][xind    ,(yind+1)] + \
              (1.-wx)*(1.-wy)*fdic[key][(xind+1),(yind+1)] 

    return var

#======================================================

def load_movie(**kwargs):

    return p3d_run('local').load_movie('all')

#======================================================

def show_energy(fname=None):
    if fname is None:
        fname = raw_input('Enter p3d.stdout file: ')

    f = open(fname, 'r')
    eng = []
    for lines in f:
        if lines.find('ENERGY') > -1 and lines.find('ENERGY:') < 0:
            eng.append(lines.split()[1:4])
    f.close()

    return np.array(eng)

#======================================================

def calc_psi(CR):
# Calculating Psi                                                                                   
    if 'bxav' in CR and 'byav' in CR:
        bx = CR['bxav']
        by = CR['byav']
    else:
        bx = CR['bx']
        by = CR['by']

    psi = 0.0*bx
    psi[0,1:] = np.cumsum(by[0,1:])*(CR['yy'][2] - CR['yy'][1])
    psi[1:,:] = psi[0,:] - np.cumsum(bx[1:,:],axis=0)*(CR['xx'][2] - CR['xx'][1])
    return psi

#======================================================

def plot_line(itcpt,dir='y',ax=None,**kwargs):
    if ax is None:
        ax = plt.gca()
    
    xarr = array(ax.get_xlim())
    yarr = array(ax.get_ylim())

    if dir == 'x':
        yarr = yarr*0.0 + itcpt
    elif dir == 'y':
        xarr = xarr*0.0 + itcpt
    else:
        print 'I dont understand what direction ' + str(dir) + \
              'is! Nothing plotted.'
        return None

    ax.plot(xarr,yarr,**kwargs)

#======================================================

def avg_movie(fname=None, 
              param=None,
              mov=None,
              ntimes=None,
              save_full=None):
    
    way = 'slow'

    if fname is None:
        fname = raw_input('Enter Save file name base: ')

    if way == 'fast':
        CC = p3d_run('local',param=param)
        print 'Loading time %i'%0
        CR = CC.load_movie('all',0)

        ntimes = CC.movie.num_of_times

        keys = CR.keys()
        keys.pop(keys.index('xx'))
        keys.pop(keys.index('yy'))

        t = time.time()
        for k in keys:
            _ = CC.load_movie(k,range(1,ntimes),mov)
            CR[k] += np.sum(_[k],axis=0) # sum along time
            CR[k+'av'] = CR[k]/1.0/ntimes
            CR[k] = _[k][-1,:,:]

        print 'TOTAL TIME: %f'%(time.time() - t)

    else:
    ## First way I tried, may be slow?
        CC = p3d_run('local',param=param)
        print 'Loading time %i'%0
        CR = CC.load_movie('all',0,mov)
        
        if ntimes is None:
            ntimes = CC.movie.num_of_times

        keys = CR.keys()
        keys.pop(keys.index('xx'))
        keys.pop(keys.index('yy'))

        t = time.time()
        for cosa in range(1,ntimes):
            print '\n==================\n' \
                    'Loading time %i' \
                  '\n==================\n'%cosa
            _ = CC.load_movie('all',cosa)
            for k in keys:
                CR[k] += _[k]
                if cosa == ntimes -1:
                    CR[k+'av'] = 1.0*CR[k]/ntimes
                    CR[k] = _[k]

        rotate_ten(CR,'pi')
        rotate_ten(CR,'pe')
        CR['tiparav'] = CR['piparav']/CR['niav']
        CR['tiperp1av'] = CR['piperp1av']/CR['niav']
        CR['tiperp2av'] = CR['piperp2av']/CR['niav']
        CR['teparav'] = CR['peparav']/CR['neav']
        CR['teperp1av'] = CR['peperp1av']/CR['neav']
        CR['teperp2av'] = CR['peperp2av']/CR['neav']

        print 'TOTAL TIME: %f'%(time.time() - t)

        CRL = {}
        CRU = {}
        ylen = len(CR['yy'])
        for k in CR:
            if k != 'yy' and k != 'xx':
                CRL[k] = np.squeeze(CR[k])[:ylen/2,:]
                CRU[k] = np.squeeze(CR[k])[ylen/2:,:]
            elif k == 'yy':
                CRL[k] = CR[k][:ylen/2]
                CRU[k] = CR[k][ylen/2:]
            else:
                CRL[k] = CR[k]
                CRU[k] = CR[k]
        
        print 'Saving lower data...'
        np.save(fname+'_lower',CRL)
        print 'Saving upper data...'
        np.save(fname+'_upper',CRU)

        if save_full:
            np.save(fname,CR)

    return CR


#======================================================

def roll_run(CR,sx=None):
    """ Roll every variable in a simulation (CR)
        in the x direction by length in indexspace (sx)
    """
    klst = ['rho','jx','jy','jz','bx','by','bz',
            'ex','ey','ez','ne','jex','jey','jez',
            'pexx','peyy','pezz','pexy','peyz','pexz',
            'ni','jix','jiy','jiz','vix','viy','viz',
            'vex','vey','vez','dene','deni',
            'pixx','piyy','pizz','pixy','piyz','pixz',
            'tepar','teperp1','tipar','tiperp1']

    kavlst = [k+'av' for k in klst]
    if sx is None:
        if CR['yy'][0] < 1.0: 
            print 'Gonna roll RIGHT!!!!!!!!!!!!'
            sx = -1*np.size(CR['xx'])/4
        else: 
            print 'Gonna roll LEFT!!!!!!!!!!!!'
            sx = np.size(CR['xx'])/4
    
    for key in CR.keys():
# Old way not super smart
#        if key.rfind('av') == len(key)-2 and len(key) > 2:
#            print 'Rolling ',key
#            CR[key] = np.roll(CR[key],sx,axis=1)
# New way a little bit smarter
        if key in klst or key in kavlst:
            print 'Rolling ',key
            CR[key] = np.roll(CR[key],sx,axis=1)

#======================================================

def readsave(restore_fname):
    if restore_fname[restore_fname.rfind('.'):] == '.npy':
        cr = np.load(restore_fname).all()
        for v in ['ni', 'ne', 'niav', 'neav']:
            if v in cr:
                cr['de'+v] = cr[v]
        return cr
    else:
        return readsav(restore_fname)

#======================================================

def colby_ct():
    incs = [0., 1.5, 2., 3., 4., 4.5, 6. ]
    incs = [c/6. for c in incs]
    cdict = {'red':  ((incs[0], 0.0, 0.0),
                      (incs[1], 0.0, 0.0),
                      (incs[2], 0.0, 0.0),
                      (incs[3], 1.0, 1.0),
                      (incs[4], 1.0, 1.0),
                      (incs[5], 1.0, 1.0),
                      (incs[6], 0.3, 0.0)),
     
             'green':((incs[0], 0.0, 0.0),
                      (incs[1], 0.0, 0.0),
                      (incs[2], 1.0, 1.0),
                      (incs[3], 1.0, 1.0),
                      (incs[4], 1.0, 1.0),
                      (incs[5], 0.0, 0.0),
                      (incs[6], 0.0, 0.0)),
     
             'blue': ((incs[0], 0.0, 0.3),
                      (incs[1], 1.0, 1.0),
                      (incs[2], 1.0, 1.0),
                      (incs[3], 1.0, 1.0),
                      (incs[4], 0.0, 0.0),
                      (incs[5], 0.0, 0.0),
                      (incs[6], 0.0, 0.0))} 

    mjet =  {'red' : ((0.0, 0.5, 0.5),
                       (0.11, 1, 1), 
                       (0.34, 1, 1), 
                       (0.50, 1, 1), 
                       (0.65, 0, 0), 
                       (1, 0, 0)),
            'green' : ((0.0, 0, 0), 
                       (0.125, 0, 0), 
                       (0.375, 1, 1), 
                       (0.50, 1, 1), 
                       (0.64, 1, 1), 
                       (0.91, 0, 0), 
                       (1, 0, 0)),
              'blue' : ((0.0, 0, 0), 
                       (0.35, 0, 0), 
                       (0.50, 1, 1), 
                       (0.66, 1, 1), 
                       (0.89, 1, 1), 
                       (1, 0.5, 0.5))}

    new_jet = {'blue' : ((0.0, 0.5, 0.5),
                   (0.11, 1, 1),
                   (0.34, 1, 1),
                   (0.50, 1, 1),
                   (0.65, 0, 0),
                   (1, 0, 0)),
        'green' : ((0.0, 0, 0),
                   (0.125, 0, 0),
                   (0.375, 1, 1),
                   (0.50, 1, 1),
                   (0.64, 1, 1),
                   (0.91, 0, 0),
                   (1, 0, 0)),
          'red' : ((0.0, 0, 0),
                   (0.35, 0, 0),
                   (0.50, 1, 1),
                   (0.66, 1, 1),
                   (0.89, 1, 1),
                   (1, 0.5, 0.5))}


    cdict = {'red' : ((0.0,   1.0, 1.0),
                      (0.03,  1.0, 1.0),
                      (0.215, 1.0, 1.0),
                      (0.4,   0.0, 0.0),
                      (0.5,   1.0, 1.0),
                      (0.586, 0.0, 0.0),
                      (0.77,  0.0, 0.0),
                      (0.954, 1.0, 1.0),
                      (1.0,   1.0, 1.0)),
 
        'green' : ((0.0,   0.0, 0.0),
                   (0.03,  0.0, 0.0),
                   (0.215, 1.0, 1.0),
                   (0.4,   1.0, 1.0),
                   (0.5,   1.0, 1.0),
                   (0.586, 1.0, 1.0),
                   (0.77,  0.0, 0.0),
                   (0.954, 0.0, 0.0),
                   (1.0,   0.0, 0.0)),
 
         'blue' : ((0.0,   0.16, 0.16),
                   (0.03,  0.0 , 0.0),
                   (0.215, 0.0 , 0.0),
                   (0.4,   0.0 , 0.0),
                   (0.5,   1.0, 1.0),
                   (0.586, 1.0 , 1.0),
                   (0.77,  1.0 , 1.0),
                   (0.954, 1.0 , 1.0),
                   (1.0,   0.75, 0.75))}

    g_modr = {'red' : ((0.0,   1.0, 1.0), #The one mike likes
                      (0.046, 1.0, 1.0), # gist_rainbow mod
                      (0.23,  0.0, 0.0),
                      (0.414, 0.0, 0.0),
                      (0.5,   1.0, 1.0),
                      (0.6,   0.0, 0.0),
                      (0.785, 1.0, 1.0),
                      (0.97,  1.0, 1.0),
                      (1.0,   1.0, 1.0)),
 
        'green' : ((0.0,   0.0, 0.0),
                   (0.046, 0.0, 0.0),
                   (0.23,  0.0, 0.0),
                   (0.414, 1.0, 1.0),
                   (0.5,   1.0, 1.0),
                   (0.6,   1.0, 1.0),
                   (0.785, 1.0, 1.0),
                   (0.97,  0.0, 0.0),
                   (1.0,   0.0, 0.0)),
 
         'blue' : ((0.0,   0.75, 0.75),
                   (0.046, 1.0 , 1.0 ),
                   (0.23,  1.0 , 1.0 ),
                   (0.414, 1.0 , 1.0 ),
                   (0.5,   1.0 , 1.0 ),
                   (0.6,   0.0 , 0.0 ),
                   (0.785, 0.0 , 0.0 ),
                   (0.97,  0.0 , 0.0 ),
                   (1.0,   0.16, 0.16))}


    cdict = {'red' : ((0.0,   1.0, 1.0),
                      (0.046, 1.0, 1.0),
                      (0.23,  0.0, 0.0),
                      (0.414, 0.0, 0.0),
                      (0.5,   1.0, 1.0),
                      (0.6,   0.0, 0.0),
                      (0.785, 1.0, 1.0),
                      (0.97,  1.0, 1.0),
                      (1.0,   1.0, 1.0)),
 
        'green' : ((0.0,   0.0, 0.0),
                   (0.046, 0.0, 0.0),
                   (0.23,  0.0, 0.0),
                   (0.414, 1.0, 1.0),
                   (0.5,   1.0, 1.0),
                   (0.6,   1.0, 1.0),
                   (0.785, 1.0, 1.0),
                   (0.97,  0.0, 0.0),
                   (1.0,   0.0, 0.0)),
 
         'blue' : ((0.0,   0.75, 0.75),
                   (0.046, 1.0 , 1.0 ),
                   (0.23,  1.0 , 1.0 ),
                   (0.414, 1.0 , 1.0 ),
                   (0.5,   1.0 , 1.0 ),
                   (0.6,   0.0 , 0.0 ),
                   (0.785, 0.0 , 0.0 ),
                   (0.97,  0.0 , 0.0 ),
                   (1.0,   0.16, 0.16))}



    ccct = LinearSegmentedColormap('ccct', cdict)
    plt.register_cmap(cmap=ccct)

    new_jet = LinearSegmentedColormap('new_jet', new_jet)
    plt.register_cmap(cmap=new_jet)

    g_modr = LinearSegmentedColormap('g_modr', g_modr)
    plt.register_cmap(cmap=g_modr)


def calc_pdf(ar,min=99999,max=99999,weight=100,inc=0,ax=0):
   if len(ar) == 0:
      print 'No array provided! Exiting!'
      return
   if min == 99999:
      min=ar.min()
   if max == 99999:
      max=ar.max()
   # If PDF of increment, then increment the array
   if inc > 0:
      ar = ar - np.roll(ar,inc,axis=ax)
   # Find the total length of data set
   arsize=reduce(operator.mul, np.shape(ar),1)
   # Find the RMS value of data set and normalize to it.
   rmsval = np.sqrt(np.mean(ar**2))
   if rmsval != 0:
      ar = ar/rmsval
   # Reshape the array to 1D & sort it.
   arr=np.reshape(ar,arsize)
   np.ndarray.sort(arr,kind='heapsort')
   # Empty arrays for output
   bins=int(arsize/weight); pdf=np.zeros(bins); binvals=np.zeros(bins)
   # Fill the bins 
   for i in range(bins):
      start=i*weight
      binvals[i] = np.mean(arr[start:start+weight])
      pdf[i] = weight/(arr[start:start+weight].max()-arr[start:start+weight].min())
   pdf = pdf/arsize
   return binvals,pdf

def rs3d(arr):
    """ Reshape an array as a 3D array (for Tulasi's stupid code)"""
    return arr.reshape(arr.shape + (1,))

def date_file_prefix():
    import datetime
    return datetime.date.today().strftime('%y.%m.%d')
