import numpy as np
import matplotlib.pyplot as plt
from scipy.io.idl import readsav
from matplotlib.ticker import AutoMinorLocator
from p3d_runs import p3d_run

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


def TestTen(var,lcl,av=''):
    if var+'xx'+av in lcl and var+'yy'+av in lcl and var+'zz'+av in lcl and var+'xy'+av in lcl and var+'yz'+av in lcl and var+'yz'+av in lcl:
        return True
    else: 
        return False

#======================================================

def ims(fdic,key,ax=None,ordflg='idl',**kwargs):
    """
    A wrapper function for imshow to do most tedious stuff for my simulations
    """
    if ax is None: ax = plt.gca()

    if ordflg == 'idl': plt_val = fdic[key]
        
    else: plt_val = fdic[key].T
# Use the dict values of xx and yy to set extent
    extent = [fdic['xx'][0],
              fdic['xx'][-1],
              fdic['yy'][0],
              fdic['yy'][-1]]

    if kwargs.has_key('cmap'): cmap=kwargs.pop('cmap')
    else:                      cmap='PuOr'

    return_ims = ax.imshow(plt_val,
                           origin='low',
                           extent=extent,
                           cmap=cmap,            # I just love this color map
                           **kwargs)

    # Giveing the plot minor tick marks
    minorLocator = AutoMinorLocator()
    ax.yaxis.set_minor_locator(minorLocator)
    minorLocator = AutoMinorLocator()           # Note the second call is so that the minor x ticks are not
    ax.xaxis.set_minor_locator(minorLocator)    # the same as the y ticks

    plt.draw()
    return return_ims

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


def load_movie(**kwargs):
    print 'not coded yet...'
    return p3d_run('local').load_movie('all')


