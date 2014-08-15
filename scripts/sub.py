import numpy as np

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

