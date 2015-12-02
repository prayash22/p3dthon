import numpy as np
from sub import *

def set_f1(fnum):

    f1 = figure(fnum)
    f1.clf()
    f1.set_size_inches(1.75,1.75)
    ax = f1.add_subplot(111)
    f1.subplots_adjust(left  =.2,
                       right =.95,
                       bottom=.15,
                       top   =.85)

    ax.xaxis.labelpad=0.
    ax.yaxis.labelpad=-7.

    return f1,ax

def calc_vars(CR):
    print 'Calculating Extra Variables to be included...'

    if 'ni' not in CR:
        CR['ni'] = CR['deniav']
    if 'ne' not in CR:
        CR['ne'] = CR['deneav']

    CR['psi'] = calc_psi(CR)
    CR['rho'] = CR['ni'] - CR['ne']
    CR['|b|'] = np.sqrt(CR['bxav']**2+CR['byav']**2+CR['bzav']**2)
    for v in ['xav','yav','zav']:
        for q,s in zip([1.,-1.],['i','e']):
            CR['v'+s+v] = q*CR['j'+s+v]/CR['n'+s]
    for v in ['xxav','yyav','zzav','xyav','xzav','yzav']:
        for s in ['i','e']:
            CR['t'+s+v] = CR['p'+s+v]/CR['n'+s]

    rotate_ten(CR,'ti')
    rotate_ten(CR,'te')

    bclip = CR['|b|'].clip(.05)


    #CR['exbx'] = (CR['ey']*CR['bz'] - CR['ez']*CR['by'])/bclip**2
    #CR['exby'] = (CR['ez']*CR['bx'] - CR['ex']*CR['bz'])/bclip**2
    #CR['exbz'] = (CR['ex']*CR['by'] - CR['ey']*CR['bx'])/bclip**2
    #CR['|exb|'] = np.sqrt(CR['exbx']**2+CR['exby']**2+CR['exbz']**2)
    dx = CR['xx'][1] - CR['xx'][0]
    CR['-vixby'] = -1.*(CR['vizav']*CR['bxav'] - CR['vixav']*CR['bzav'])
    CR['-vexbyav'] = -1.*(CR['vezav']*CR['bxav'] - CR['vexav']*CR['bzav'])
    CR['vdgvyav'] = -1./25.*(CR['vexav']*(np.roll(CR['veyav'],-1,axis=1) - 
                                                  CR['veyav'])/dx +
                             CR['veyav']*(np.roll(CR['veyav'],-1,axis=0) - 
                                                   CR['veyav'])/dx)

    CR['gpynav'] = -1.0/dx*((np.roll(CR['pexyav'],-1,axis=1) - CR['pexyav']) +
                   (np.roll(CR['peyyav'],-1,axis=0) - CR['peyyav']))/CR['deneav']

    CR['-tezzav'] = -1.0*CR['tezzav']
    CR['-bzav'] = -1.0*CR['bzav']
    CR['-vizav'] = -1.0*CR['vizav']
    CR['-vezav'] = -1.0*CR['vezav']

    #CR['epar'] = (CR['ex']*CR['bx'] +
    #              CR['ey']*CR['by'] +
    #              CR['ez']*CR['bz'])/CR['|b|']
    #CR['jpar'] = (CR['jx']*CR['bx'] +
    #              CR['jy']*CR['by'] +
    #              CR['jz']*CR['bz'])/CR['|b|']

    CR['pb']   = CR['|b|']**2/2.
    CR['ptot'] = CR['pb'] + CR['piyyav'] + CR['peyyav']

    # Flipping values!!!!

def get_mpline(CR):
    mp = 0.0*CR['bxav'][0,:]
    for c in range(np.size(CR['bxav'][0,:])):
        mp[c] = CR['yy'][abs(CR['bxav'][:,c]).argmin()]

    return mp

def calc_EvxB(CR):

    CR['vixav'] = CR['jixav']/CR['deniav'] 
    CR['viyav'] = CR['jiyav']/CR['deniav'] 
    CR['vizav'] = CR['jizav']/CR['deniav'] 
    CR['vexav'] = -1.*CR['jexav']/CR['deneav'] 
    CR['veyav'] = -1.*CR['jeyav']/CR['deneav'] 
    CR['vezav'] = -1.*CR['jezav']/CR['deneav'] 
    CR['EviBy'] = CR['eyav'] + \
                 (CR['vizav']*CR['bxav'] - CR['vixav']*CR['bzav'])
    CR['EveBy'] = CR['eyav'] + \
                 (CR['vezav']*CR['bxav'] - CR['vexav']*CR['bzav'])
    CR['De']    = 0.0*CR['bzav']
    CR['Ae']    = CR['De']
    CR['EviBz'] = CR['ezav'] + \
                 (CR['vixav']*CR['byav'] - CR['viyav']*CR['bxav'])
    CR['EveBz'] = CR['ezav'] + \
                 (CR['vexav']*CR['byav'] - CR['veyav']*CR['bxav'])

    xp = 26.275 + (CR['xx'][-1] + CR['xx'][-0])/4.0
    yp = 37.875 

    CR['xx'] -= xp
    CR['yy'] -= yp


