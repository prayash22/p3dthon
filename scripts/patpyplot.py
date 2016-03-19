from plot_defs import *
from matplotlib.backends.backend_pdf import PdfPages
from sub import calc_psi


# Assume Movie Dat has been loaded
# and is called CR

# Calc extra variables
def make_pat_plot(CR,fname=None):

    plt.ioff()
    _CR = dict(CR)
    for k in _CR.keys():
        if len(k) > 2 and k.rfind('av') > 0: 
            _CR[k[:-2]] = _CR[k]
    calc_extra_vars(_CR)
    cut_locs,psi_lvls = cut_n_cont_locs(_CR)

    vars_2D = []
    pg_cmap = []

# Page #1
    vars_2D.append(['ni','rho',
                    'bx','ex',
                    'by','ey',
                    'bz','ez',
                    '|b|','jx',
                    'jy','jz',])
    pg_cmap.append(['gist_heat','gist_heat',
                    'bwr','bwr',
                    'bwr','bwr',
                    'bwr','bwr',
                    'gist_heat','bwr',
                    'bwr','bwr'])
# Page #2
    vars_2D.append(['jix','jex',
                    'jiy','jey',
                    'jiz','jez',
                    'vix','vex',
                    'viy','vey',
                    'viz','vez'])
    pg_cmap.append(['bwr','bwr',
                    'bwr','bwr',
                    'bwr','bwr',
                    'bwr','bwr',
                    'bwr','bwr',
                    'bwr','bwr'])
# Page #3
    vars_2D.append(['tixx','texx',
                    'tiyy','teyy',
                    'tizz','tezz',
                    'tixy','texy',
                    'tixz','texz',
                    'tiyz','teyz'])
    pg_cmap.append(['gist_heat','gist_heat',
                    'gist_heat','gist_heat',
                    'gist_heat','gist_heat',
                    'gist_heat','gist_heat',
                    'gist_heat','gist_heat',
                    'gist_heat','gist_heat'])
# Page #4
#    vars_2D.append(['vipar','vepar', # This was the old set up
#                    'viperp','veperp',
    vars_2D.append(['tipar','tepar',
                    'tiperp1','teperp1',
                    'viperp','vepar',
                    'exbx','exby',
                    'exbz','|exb|',
                    'epar', 'jpar'])
    pg_cmap.append(['gist_heat','gist_heat',
                    'gist_heat','gist_heat',
                    'bwr','bwr',
                    'bwr','bwr',
                    'bwr','gist_heat',
                    'bwr','bwr'])

# This is bascily the dict that that coresponds to each location 
# of the 1D cuts
    vars_1D = [{k : _CR[k] for k in ('bx', 'by', 'bz', '|b|')},
               {k : _CR[k] for k in ('ni', 'ne', 'rho')},
               {k : _CR[k] for k in ('ex', 'ey', 'ez')},
               {k : _CR[k] for k in ('piyy', 'peyy', 'pb', 'ptot')},
               {k : _CR[k] for k in ('exbx', 'exby', 'exbz', '|exb|')},
               {k : _CR[k] for k in ('jx', 'jy', 'jz')},
               {k : _CR[k] for k in ('vix', 'viy', 'viz')},
               {k : _CR[k] for k in ('vex', 'vey', 'vez')},
               {k : _CR[k] for k in ('tixx', 'tiyy', 'tizz')},
               {k : _CR[k] for k in ('texx', 'teyy', 'tezz')},
               {k : _CR[k] for k in ('tipar', 'tiperp1')},
               {k : _CR[k] for k in ('tepar', 'teperp1')}]

    pg = 0
    plt.ioff()
    if fname is None: fname = raw_input('Please enter run name: ')
    with PdfPages(fname+'_patplots.pdf') as pdf:
# First we plot the 2D's
        for dvar,cmaps in zip(vars_2D,pg_cmap):
            pg += 1
            print 'Plotting pg%i...'%(pg)
            fig,axs = set_fig()
            for ax,v,cmap in zip(axs,dvar,cmaps):
                _ = plot_2D(ax, _CR, v,
                            cmap=cmap,
                            psi_lvls=psi_lvls,
                            cut_locs=cut_locs)

                if cmap is 'bwr':
                    imax = np.max(np.abs(np.array(_[0].get_clim())))
                    _[0].set_clim(-1.*imax,imax)

            close_fig(pdf)

# Now 1D plots
        for loc in cut_locs:
            pg += 1
            print 'Plotting pg%i...'%(pg)
            fig,axs = set_fig()
            for (ax,d) in zip(axs,vars_1D):
                plot_1D(ax,_CR,d,loc=loc)

            close_fig(pdf)
        print 'All done, just finishing up...'


#    # We can also set the file's metadata via the PdfPages object:
#    d = pdf.infodict()
#    d['Title'] = 'Multipage PDF Example'
#    d['Author'] = u'Jouni K. Sepp\xe4nen'
#    d['Subject'] = 'How to create a multipage pdf file and set its metadata'
#    d['Keywords'] = 'PdfPages multipage keywords author title subject'
#    d['CreationDate'] = datetime.datetime(2009, 11, 13)
#    d['ModDate'] = datetime.datetime.today()


