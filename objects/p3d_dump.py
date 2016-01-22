#######################################################################
#######################################################################
#                                                                     #
#                  Python Progs :  p3d_dump.py                        #
#                  Aruthor      :  Colby Haggerty                     #
#                  Date         :  2014.08.15                         #
#                                                                     #
#                                                                     #
#######################################################################
### Discription:
#

import os
import sys 
import datetime
import numpy as np
import struct
import glob
import pdb
from scipy.io.idl import readsav
from scipy.ndimage import gaussian_filter

class p3d_dump(object):
    """p3d_run object 

        p3d_movie class is a movie object to be handled by p3d_run. This object
        keeps all of the information for a  given movie. We have seperated this from the
        p3d_object to keep things neatter and cleaner. The plan is then apply this same
        idea to the dump files. NOTE there may be a time trade off issue between this
        and doing every thing together, but i dont think that will be an issue.  #
    """

    def __init__(self, dump_path, param_dict, dump_num=-1): 
        """ Initilazition Routine for the p3d_run object

        Fill in method discription

        @param dump_path: the path to the set of movies you are interested in
        @param dump_num    : The number identive which set of movies you want to look at (which restart of the run)
        @param param_dict   : The dictonary of the param file associated with this run. It should be passed from the p3d_run object

        @return: @todo

        Example  :

        Creation
        :
        2013-05-08
        11:21:19.671182
        """
        self._set_local_flag = False

        if dump_path.strip()[-1] == '/': self.dump_path = dump_path[:-1].strip()
        else: self.dump_path = dump_path
        self.param_dict = param_dict
        if dump_num < 0: dump_num = self._dump_num_options()
        self.dump_num = self._num_to_ext(dump_num) # Code _movie_num_opt

        if 'stat_ions' in self.param_dict: self.species = ['i','e'] 
        else: self.species = ['i','e'] #all species contained within dump files. will need to recode for multi species

# this method assumes that particles are already laoded
    def moments(self):
        num_es = np.size(self.particles['e'])
        box_area = self._dx[0]*self._dx[1]
        dx = self.param_dict['lx']/self.param_dict['nx']/self.param_dict['pex']
        dy = self.param_dict['ly']/self.param_dict['ny']/self.param_dict['pey']
 #sloopying coding, kmalakit uses n1 and marc uses n0
        n0 = float(self.param_dict.get('n0',self.param_dict['n1']))
        ppdi = self.param_dict['ppg']/dx/dy/n0

        ne_loc = num_es/box_area/ppdi
        
        uex_loc = np.mean(self.particles['e']['vx'])
        uey_loc = np.mean(self.particles['e']['vy'])
        uez_loc = np.mean(self.particles['e']['vz'])

        pexx_loc = np.mean((self.particles['e']['vx'] - uex_loc)**2)
        peyy_loc = np.mean((self.particles['e']['vy'] - uey_loc)**2)
        pezz_loc = np.mean((self.particles['e']['vz'] - uez_loc)**2)


        print 'local ne    : ',ne_loc
        print 'local uex   : ',uex_loc
        print 'local uey   : ',uey_loc
        print 'local uez   : ',uez_loc
        print 'local pexx  : ',pexx_loc
        print 'local peyy  : ',peyy_loc
        print 'local pezz  : ',pezz_loc
        pass
        

    def read_dump_file(self,dump_index=-1,fields=False,verbose=False):
        """
#-----------------------------------------------------------------------------
#   Method      : read_dump_file
# 
#   Args        : dump_index (is the middle 3 digets of the dump file. This 
#                 orseponds to the processors)
#               : dump_num (is the last 3 digets of the dump file. This 
#                 spesifies which run)
# 
#   BIG NOTE    : This assumes that we are looking at 2 D dump files! so file 1 
#                 is field file 
# 
#   Comment     : 
#----------------------------------------------------------------------------
        """
        
        return_dict = {}

        if fields:
            dump_index='001'

        dump_index = self._num_to_ext(dump_index)
        fname = '%s/p3d-%s.%s'%(self.dump_path,dump_index,self.dump_num)
        #qc#print 'fle name = '+fname

# Make sure that the file exsists and suguset alternetives to fix it
        try:
            self.open_dump_file = open(fname, "rb")
        except IOError as e:
            print "I/O error({0}): {1}".format(e.errno, e.strerror)
            print "ERROR: Could not open file. " + fname
            print "Possible Reasons:"
            print "     File does not exist"
            print "     File extention is wrong"
            print "     run_list.dat is not set properly (i.e. looking in the wrong directory)"
            return -1


        self._pad_headder()
# First we need to deal with special case of if it is a field file
# If the dump file we are trying to read is the first, we need to first read the fields and then the particles 
        if dump_index == '001':
            return_dict = self._get_fields()
            self.fields = return_dict
# Now we open 
        if not fields:
            return_dict = dict(self._get_particles().items() + return_dict.items())
            
        self.open_dump_file.close()
        return return_dict


    def _get_fields(self):
        #qc#print 'Reading fields: ex, ey, ez, bx, by, bz'
        self.dump_field_dict={}
        for field_var in ['ex','ey','ez','bx','by','bz']:
            dump_dat=[]
            for pey_current in range(self.py):
                pad_head = struct.unpack('<i',self.open_dump_file.read(4))[0]
                dump_dat.append(np.fromfile(self.open_dump_file,dtype='float64',count=pad_head/8)) # Float size times number of floats
                pad_butt = struct.unpack('<i',self.open_dump_file.read(4))[0]
            dump_dat = np.concatenate(dump_dat)
            dump_dat.shape = (self.py,self.px)
            self.dump_field_dict[field_var] = dump_dat
        return self.dump_field_dict
# You should relly wtright in a check to compare pex from param to what you find here
# Like they do in p3d

    def _pad_headder(self):
        #qc#print "Reading Header from dump file"
        header_chunksize=44 # Size of the header in bytes
        header_binary = self.open_dump_file.read(header_chunksize)
        #if sys.byteorder == 'little': ic = '<' # Littel Endian (do I even need to do this?) 
        #else ic = '>' # Big Endian
        header = struct.unpack('<idd6i', header_binary) #Int, Double Double, 6 Ints
        #qc#print "              : "+str(header)
        self.time = header[1]
        self.n_avg = header[2]
        self.px = header[3]
        self.py = header[4]
        self.pz = header[5]
        self.bufsize = header[6]
        self.nchannels = header[7]
        #qc#print 'time = %f'%(self.time)
        #qc#print 'n_avg = %f'%(self.n_avg)
        #qc#print 'px =  %i, %i'%(self.px,self.param_dict['nx']*self.param_dict['pex'])
        #qc#print 'py =  %i, %i'%(self.py,self.param_dict['ny']*self.param_dict['pey'])
        #qc#print 'pz =  %i, %i'%(self.pz,self.param_dict['nz']*self.param_dict['pez'])
        #qc#print 'bufsize = %i, %i'%(self.bufsize, self.param_dict['bufsize'])
        #qc#print 'nchannels = %i, %i\n'%(self.nchannels, self.param_dict['nchannels'])

    def _get_particles(self):
        data_type = np.dtype([('x', 'float32'), ('y', 'float32'), ('z', 'float32'),('vx', 'float32'), ('vy', 'float32'), ('vz', 'float32')])
        #all_particles = [] 
        all_particles = {} 
        for species in self.species: 
            pad_head = struct.unpack('<i',self.open_dump_file.read(4))[0]
            #c#verboseprint('Note sure about the point of this number? ' +str(struct.unpack('<i',self.open_dump_file.read(4))[0]))
            #print 'Not sure about the point of this number? ' +str(struct.unpack('<i',self.open_dump_file.read(4))[0])
            str(struct.unpack('<i',self.open_dump_file.read(4))[0]) # I think you need this comand so that it will
                                                                    # unpack what ever this value is
            pad_butt = struct.unpack('<i',self.open_dump_file.read(4))
            all_sub_species = []
            for current_sub_proc in range(self.param_dict['pey']*int(np.round(self.param_dict['pex']/self.param_dict['nchannels']))): 
                pad_head = struct.unpack('<i',self.open_dump_file.read(4))[0]
                #print 'pad_head = '+str(pad_head)
                number_of_part_on_pe = struct.unpack('<i',self.open_dump_file.read(4))[0]
                pad_butt = struct.unpack('<i',self.open_dump_file.read(4))
                dump_dat=[]
                bufsize_lastcase = number_of_part_on_pe % self.bufsize
                #c#verboseprint('Reading from proc number: '+str(current_sub_proc)+' Number of part on sub proc: '+str(number_of_part_on_pe))
                #vbct# print 'Reading from proc number: '+str(current_sub_proc)+' Number of part on sub proc: '+str(number_of_part_on_pe)
                if abs(1.0*number_of_part_on_pe/self.bufsize -  round(number_of_part_on_pe/self.bufsize)) < 1./self.bufsize:
                    rgn = number_of_part_on_pe/self.bufsize-1
                else:
                    rgn = number_of_part_on_pe/self.bufsize
                for current_sub_buffer in range(rgn): 
                    pad_head = struct.unpack('<i',self.open_dump_file.read(4))[0]
                    #print 'Reading Buffer number: '+str(current_sub_buffer)+' Size of next set of bytes: '+str(pad_head)
# Colby Maybe try switching these two to see if it runs faster. The two should be Equivelent
#   1:
#;#                all_particles_from_file = np.append(all_particles_from_file,np.fromfile(f,dtype=dt,count=pad_head/(4*6))) # Float size times number of floats
#   2:
#;#                dump_dat = np.fromfile(f,dtype=dt,count=pad_head/(4*6)) # Float size times number of floats
#;#                all_particles_from_file = np.append(all_particles_from_file,dump_dat)
#   3:
                    dump_dat.append(np.fromfile(self.open_dump_file,dtype=data_type,count=pad_head/(4*6))) # Float size times number of floats
#   end
                    pad_butt = struct.unpack('<i',self.open_dump_file.read(4))
                    #pad_head = struct.unpack('<i',self.open_dump_file.read(4))[0]
                    pad_head = struct.unpack('<i',self.open_dump_file.read(4))
                    #print 'pad_head = '+str(pad_head)
                    pad_head = pad_head[0]
                    #print 'Reading Buffer number: '+str(current_sub_buffer)+' Size of next set of bytes: '+str(pad_head)+' !SKIPPED!'
                    
# I think this run might not have any tags
                    np.fromfile(self.open_dump_file,dtype='int64',count=pad_head/(8)) # 1 int8 and we probobly dont need the tag
                    pad_butt = struct.unpack('<i',self.open_dump_file.read(4))
                # Special Case of the particles left over
                #   If you are looking hear to figure out an issue it could be that we do not check to make sure that
                #   this happens in reading the dump file. It is possible that the number of particles excatly fills
                #   up the buffer and you dont have this extra case. But I dont think this is likly to happen (1/ bufsize)
                #   First read the number of total particles on this PE, This will tell us how much our next byte size should be
                #print 'Number of particles in the final buffer: '+str(bufsize_lastcase)
                #   Next we read the data with the special size. NOTE this could be done in a cleaner way
                pad_head = struct.unpack('<i',self.open_dump_file.read(4))[0]
                #print 'Reading Buffer number: '+str(number_of_part_on_pe/bufsize+1)+' Size of next set of bytes: '+str(pad_head)
                temp_dat = np.fromfile(self.open_dump_file,dtype=data_type,count=pad_head/(4*6)) # Float size times number of floats
                pad_butt = struct.unpack('<i',self.open_dump_file.read(4))[0]
# Trim all of the extra zeros
                dump_dat.append(temp_dat[0:bufsize_lastcase])
# Appending temp dump_dat to all_particles
                #all_particles_from_file = np.append(all_particles_from_file,dump_dat)
# Now skip over the tags
                pad_head = struct.unpack('<i',self.open_dump_file.read(4))[0]
                #print 'Reading Buffer number: '+str(number_of_part_on_pe/bufsize+1)+' Size of next set of bytes: '+str(pad_head)+' !SKIPPED!'
                np.fromfile(self.open_dump_file,dtype='int64',count=pad_head/(8)) # 1 int8 and we probobly dont need the tag
                pad_butt = struct.unpack('<i',self.open_dump_file.read(4))
                all_sub_species.append(np.concatenate(dump_dat))
            all_particles[species] = all_sub_species
        return all_particles

    def _dump_num_options(self):
        choices = glob.glob(self.dump_path+'/p3d-001.*')
        if len(choices) == 0: 
            print '!!! WARNING: the direcotry we are looking in does not have any moive.bz.XXX so we are crashing'
            return -1
        for var in range(np.size(choices)):
            choices[var] = choices[var][-3:]
        print 'Select from the following possible dump numbers: \n'+str(choices),
        movie_num_int = raw_input()
        return int(movie_num_int)



    def _dump_num_options(self):
        choices = glob.glob(self.dump_path+'/p3d-001.*')
        if len(choices) == 0: 
            print '!!! WARNING: the direcotry we are looking in does not have any moive.bz.XXX so we are crashing'
            return -1
        for var in range(np.size(choices)):
            choices[var] = choices[var][-3:]
        print 'Select from the following possible dump numbers: \n'+str(choices),
        movie_num_int = raw_input()
        return int(movie_num_int)

    def vdist_2d(self,
                 r0=[0.5,0.5],
                 dx=[1.0,1.0],
                 par=False,
                 Bvec=False,
                 pitch=False,
                 pizza=False,
                 OneD=False,
                 **kwargs):

        """ Generates differnent 2-Dimensional histograms for particles 
            Please god colby add a doc string!!!
        """

        self._par   = par
        self._pizza = pizza
        self._Bvec  = Bvec

        if isinstance(dx,float):
            dx = [dx,dx] # Square box. Does this change dx?

# Turn this into a method
#%%%%%%%%%%%%%%%
        if OneD: # Make dy = ly
            r0[1] = self.param_dict['ly']/2.0
            dx = [dx[0], self.param_dict['ly']/2.0]
        if not hasattr(self,'_r0') or not hasattr(self,'particles'):
            self._r0 = r0[:]
            self._dx = dx[:]
# Calling get particles in box to make the vdist
            #qc#print 'Reading Ions and Electrons from the Dump File'
            self.particles = self._part_in_box(r0,dx)
        else:
            if (self._r0[0] == r0[0]) and \
               (self._r0[1] == r0[1]) and \
               (self._dx[0] == dx[0]) and \
               (self._dx[1] == dx[1]):

                print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
                print '~~~ Same r0 and dx found, using old particles ~~~'
                print '~~~ r0 = '+str(r0)+' _r0 = '+str(self._r0)  
                print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

            else:
                self._r0 = r0[:]
                self._dx = dx[:]
                #qc#print 'Reading Ions and Electrons from the Dump File'
                self.particles = self._part_in_box(r0,dx)

        #print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
        #print kwargs
        #if kwargs.pop('v0_frame',False):
        #    print 'THE SHIFFFFFFTTTTINGGGGGG!!!!!!'
        #    self._shift_frame()

        if not kwargs.has_key('bins'): kwargs['bins']=51
        if par or Bvec or pitch or pizza:
            #qc#print 'Reading in the Fields form the Dump File'
            self.dump_field_dict = self.read_dump_file(fields=True)
            if pitch:
                return_hist = self._vdist_pitch(par=par,**kwargs)
            elif pizza:
                return_hist = self._vdist_pizza(par,**kwargs)
            elif par and not pitch and not pizza:
                return_hist = self._vdist_2d_par(**kwargs)
            else:
                return_hist = self._vdist_2d(**kwargs)
                #qc#print 'Interpolating the Bfield at the given r0 value'
                return_hist['B'] = np.array([
                   self.interp_field(self.dump_field_dict['bx']),
                   self.interp_field(self.dump_field_dict['by']), 
                   self.interp_field(self.dump_field_dict['bz'])])
        else:
            return_hist = self._vdist_2d(**kwargs)

        self._box = [self._r0[0]-self._dx[0]/2.,
                     self._r0[0]+self._dx[0]/2.,
                     self._r0[1]-self._dx[1]/2.,
                     self._r0[1]+self._dx[1]/2.]

        return_hist['box'] = self._box

        return return_hist





#wut?    def _vdist_pitch(self,**kwargs):
#wut?        bx_interp = self.interp_field(self.dump_field_dict['bx'],self.param_dict['lx'],self.param_dict['ly'],self._r0)
#wut?        by_interp = self.interp_field(self.dump_field_dict['by'],self.param_dict['lx'],self.param_dict['ly'],self._r0)
#wut?        bz_interp = self.interp_field(self.dump_field_dict['bz'],self.param_dict['lx'],self.param_dict['ly'],self._r0)
#wut?        bmag_interp = (bx_interp**2+by_interp**2+bz_interp**2)**.5
#wut?
#wut?        if by_interp > 0.:
#wut?            b_perp1x = 0.
#wut?            b_perp1y = -1.*bz_interp/(bz_interp**2 + by_interp**2)**(.5)
#wut?            b_perp1z = by_interp/(bx_interp**2 + by_interp**2)**(.5)
#wut?        else:
#wut?            b_perp1x = 0.
#wut?            b_perp1y = bz_interp/(bz_interp**2 + by_interp**2)**(.5)
#wut?            b_perp1z = -1.*by_interp/(bx_interp**2 + by_interp**2)**(.5)
#wut?
#wut?        b_perpmag = (b_perp1x**2+b_perp1y**2+b_perp1z**2)**.5
#wut?        b_perp1x = b_perp1x/b_perpmag
#wut?        b_perp1y = b_perp1y/b_perpmag
#wut?        b_perp1z = b_perp1z/b_perpmag
#wut?
#wut?        b_perp2x = (by_interp*b_perp1z - bz_interp*b_perp1y)
#wut?        b_perp2y = (bz_interp*b_perp1x - bx_interp*b_perp1z)
#wut?        b_perp2z = (bx_interp*b_perp1y - by_interp*b_perp1x)
#wut?        b_perpmag = (b_perp2x**2+b_perp2y**2+b_perp2z**2)**.5
#wut?        b_perp2x = b_perp2x/b_perpmag
#wut?        b_perp2y = b_perp2y/b_perpmag
#wut?        b_perp2z = b_perp2z/b_perpmag
#wut?
#wut?        velo={}
#wut?        for species in self.species:
#wut?            velo[species] = {}
#wut?
#wut?            velo[species]['par']   = (bx_interp*self.particles[species]['vx']+by_interp*self.particles[species]['vy']+bz_interp*self.particles[species]['vz'])/bmag_interp
#wut?            velo[species]['perp1'] = self.particles[species]['vx']*b_perp1x+self.particles[species]['vy']*b_perp1y+self.particles[species]['vz']*b_perp1z
#wut?            velo[species]['perp2'] = self.particles[species]['vx']*b_perp2x+self.particles[species]['vy']*b_perp2y+self.particles[species]['vz']*b_perp2z
#wut?
#wut?        return_hist_dict = {}
#wut?        for species in self.species:
#wut?            return_hist_dict[species] = []
#wut?
#wut?        for species in velo.keys():
#wut?            H, xedges, yedges = np.histogram2d(velo[species]['par'],velo[species]['perp1'],**kwargs)
#wut?# H needs to be rotated and flipped
#wut?            H = np.rot90(H)
#wut?            H = np.flipud(H)
#wut?            return_hist_dict[species].append(['Parallel vs Perp 1 (+zy)',H,xedges,yedges])
#wut?
#wut?            H, xedges, yedges = np.histogram2d(velo[species]['par'],velo[species]['perp2'],**kwargs)
#wut?            H = np.rot90(H)
#wut?            H = np.flipud(H)
#wut?            return_hist_dict[species].append(['Parallel vs Perp 2',H,xedges,yedges])
#wut?
#wut?            H, xedges, yedges = np.histogram2d(velo[species]['perp1'],velo[species]['perp2'],**kwargs)
#wut?            H = np.rot90(H)
#wut?            H = np.flipud(H)
#wut?            return_hist_dict[species].append(['Perp 1 (+zy) vs Perp 2',H,xedges,yedges])
#wut?
#wut?# Mask zeros
#wut?         #Hmasked = np.ma.masked_where(H==0,H)
#wut?
#wut?        return return_hist_dict


    def _vdist_pitch(self,par=False,energy=False,**kwargs): 
        """ Pitch party!!!!
        """
# I have a lot of code that can return a 2D distrubution function for a givn pitch angle range
# but I don't think that this will ever be usefull... I should remove it
        pa  = kwargs.pop('pa',90.)
        dpa = kwargs.pop('dpa',5.)
        wax = kwargs.pop('wax',0)
        v0_frame = kwargs.pop('v0_frame',None)

        if energy: par = True; kwargs['normed'] = True
        self._wax = wax
##############################################

        if par:
            b_interp = np.array([ self.interp_field(self.dump_field_dict['bx'],self.param_dict['lx'],self.param_dict['ly'],self._r0),
                                  self.interp_field(self.dump_field_dict['by'],self.param_dict['lx'],self.param_dict['ly'],self._r0),
                                  self.interp_field(self.dump_field_dict['bz'],self.param_dict['lx'],self.param_dict['ly'],self._r0)])

            e_interp = np.array([ self.interp_field(self.dump_field_dict['ex'],self.param_dict['lx'],self.param_dict['ly'],self._r0),
                                  self.interp_field(self.dump_field_dict['ey'],self.param_dict['lx'],self.param_dict['ly'],self._r0),
                                  self.interp_field(self.dump_field_dict['ez'],self.param_dict['lx'],self.param_dict['ly'],self._r0)])

            exb = np.cross(e_interp,b_interp)/sum(b_interp**2)
            if abs(np.sqrt(sum(exb**2))) < np.spacing(10):
                exb = np.cross(np.array([0.,1.,1.]),b_interp)
                exb = exb/np.sqrt(sum(exb**2))
            else:
                exb = exb/np.sqrt(sum(exb**2))

            b_interp = b_interp/np.sqrt(sum(b_interp**2))

            bxexb = np.cross(b_interp,exb) 
            #print b_interp
            #print exb
            #print bxexb

        velo={}
        return_hist_dict = {}
        self.subpart = {}
        #print 'Mucking this up!!'
        #Npart = np.shape(self.particles['e']['vx'])
        #self.particles['e']['vx'] = 1.*np.random.normal(scale=2.5,size=Npart)
        #self.particles['e']['vy'] = 1.*np.random.normal(scale=2.5,size=Npart)
        #self.particles['e']['vz'] = 1.*np.random.normal(scale=2.5,size=Npart)
        #print 'Done Mucking this up!!!'
        for species in self.species:
       
            if par:
                v0 = (b_interp[0]*self.particles[species]['vx']+
                      b_interp[1]*self.particles[species]['vy']+
                      b_interp[2]*self.particles[species]['vz'])

                v1 = (exb[0]*self.particles[species]['vx']+
                      exb[1]*self.particles[species]['vy']+
                      exb[2]*self.particles[species]['vz'])

                v2 = (bxexb[0]*self.particles[species]['vx']+
                      bxexb[1]*self.particles[species]['vy']+
                      bxexb[2]*self.particles[species]['vz'])
# Muck ===============
#               print 'Mucking this up!!!'
#               Npart = 1000000
#               v0 = 1.*np.random.normal(scale=np.sqrt(50.),size=Npart)
#               v1 = 1.*np.random.normal(scale=np.sqrt(50.),size=Npart)
#               v2 = 1.*np.random.normal(scale=np.sqrt(50.),size=Npart)
# Muck ===============
            else:
                v0 = self.particles[species]['vx']
                v1 = self.particles[species]['vy']
                v2 = self.particles[species]['vz']

            if wax == 0: # wax is just which axis defines the plane we are looking at
                vp0 = v1
                vp1 = v2
                vax = v0
            elif wax == 1:
                vp0 = v0
                vp1 = v1
                vax = v2
            elif wax == 2:
                vp0 = v0
                vp1 = v1
                vax = v2
            else :
                print ''
                print 'The plane axis is out of bounds fo wax = ',wax
                print 'vdist_2d is crashing!!!!'
                print ''
            
            if v0_frame:
                if energy:
                    vax = vax - np.mean(vax)
                    vp0 = vp0 - np.mean(vp0)
                    vp1 = vp1 - np.mean(vp1)
                else:
                    print '~'*80
                    print 'Shifting vax by ',np.mean(vax)
                    vax = vax - np.mean(vax)

            #self.vax = vax
            #self.vp = np.sqrt(vp0**2+vp1**2)
            #self.exb = exb
            #self.bxexb = bxexb
            # This means field aligned is 0, perp is 90 and anti aligned is 180
            pitch_angle = -1.0*(np.arctan(vax/np.sqrt(vp0**2+vp1**2))/np.pi*180. - 90.)
            self.ptc = pitch_angle

# Not sure which of these two is right
            #pitch_angle = np.arccos(vpar/vmag)/np.pi*180.
            #pitch_angle = np.arctan(np.sqrt((vmag**2 - vpar**2)/2.)/vpar)/np.pi*180. + 90.

            subpartind = np.where(abs(pitch_angle - pa) < dpa/2.)

            #subpart = self.particles[species][subpartind]
            #self.subpart[species] = subpart
            #subpitch_angle = pitch_angle[subpartind] 

            #print 'Total %s in pitch angle range are: %i'%(species,len(subpart))
            return_hist_dict[species] = []

#This is silly just do the whole distro right here
            if energy:
                if species == 'e':
                    reltv = True
                    if reltv:
                        v2 = (self.particles[species]['vx']**2+
                              self.particles[species]['vy']**2+
                              self.particles[species]['vz']**2)

                        KE = self.param_dict['m_e']* self.param_dict['c_2']*\
                             (1./np.sqrt(1. - v2/self.param_dict['c_2']) - 1.)

                    else:
                        KE = self.param_dict['m_e']/2.0*\
                            (self.particles[species]['vx']**2+
                             self.particles[species]['vy']**2+
                             self.particles[species]['vz']**2)
# Muck ===============
#                   print 'More Mucking!!!'
#                   KE = self.param_dict['m_e']/2.0*(v0**2+
#                                                    v1**2+
#                                                    v2**2)
#                   
# Muck ===============
                else:
                    KE = 1.0/2.0*(self.particles[species]['vx']**2+
                                  self.particles[species]['vy']**2+
                                  self.particles[species]['vz']**2)
# Muck ===============
#                   print 'More Mucking!!!'
#                   KE = 1./2.0*(v0**2+
#                                v1**2+
#                                v2**2)
#                   
# Muck ===============
                #H,xedges = np.histogram(KE,**kwargs)
                #return_hist_dict[species].append(H)
                #return_hist_dict[species].append((xedges[:-1]+xedges[1:])/2.0)


                print 'pa_size = ',np.size(pitch_angle)
                print 'KE_size = ',np.size(KE)
                H,xedges,yedges = np.histogram2d(pitch_angle,KE,**kwargs)
                xx,yy = np.meshgrid(yedges,xedges)
                eng = (xx[1:,1:]+xx[:-1,:-1])/2.
                ynorm = (np.cos(yy[:-1,:-1]/180.*np.pi) - np.cos(yy[1:,1:]/180.*np.pi))
                ynorm = ynorm/sum(ynorm)
                if species == 'e':
                    self.ynorm=ynorm
                    self.pa = pitch_angle
                    self.KE = KE
                    self.eng = eng
                    self.H = H
                    self.yyy = yy

# We need to account for a geometric factor because
# E is basicly v^2, so we are binning circles?

                #H = H/np.sqrt(eng)
                #return_hist_dict[species].append(H/ynorm*eng**2.*np.size(pitch_angle))

                Ebar = eng/self.param_dict['m_e']/self.param_dict['c_2']
                rel_vel = np.sqrt((Ebar**2 + 2.*Ebar)/(Ebar**2 + 2.*Ebar + 1.))
                rel_vel = rel_vel*np.sqrt(self.param_dict['c_2'])

                return_hist_dict[species].append(H/ynorm*eng*rel_vel*np.size(pitch_angle))

                return_hist_dict[species].append(xedges)
                return_hist_dict[species].append(yedges)

                
            else:
                #Colby devide by bin size!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                H, xedges, yedges = np.histogram2d(vp0[subpartind],
                                                   vp1[subpartind],
                                                   **kwargs)
# Muck ===============
# Mucking!!!!!
#                H = gaussian_filter(H, sigma=1.)
# Muck ===============

                xx,yy = np.meshgrid(yedges,xedges)
                print pa,dpa

                #int_cone = 1./18.*(
                self.inc = 1./18.*(
                       (-xx[:-1,:-1]**3+6.*xx[:-1,:-1]*yy[:-1,:-1]*np.sqrt(xx[:-1,:-1]**2+yy[:-1,:-1]**2) + 
                        3.*yy[:-1,:-1]**3*np.log(np.sqrt(xx[:-1,:-1]**2+yy[:-1,:-1]**2) + xx[:-1,:-1] + np.spacing(4)) + 
                        3.*xx[:-1,:-1]**3*np.log(np.sqrt(xx[:-1,:-1]**2+yy[:-1,:-1]**2)+yy[:-1,:-1] + np.spacing(4))) -
                       (-xx[1:,1:]**3+6.*xx[1:,1:]*yy[:-1,:-1]*np.sqrt(xx[1:,1:]**2+yy[:-1,:-1]**2) + 
                        3.*yy[:-1,:-1]**3*np.log(np.sqrt(xx[1:,1:]**2+yy[:-1,:-1]**2) + xx[1:,1:] + np.spacing(4)) + 
                        3.*xx[1:,1:]**3*np.log(np.sqrt(xx[1:,1:]**2+yy[:-1,:-1]**2)+yy[:-1,:-1] + np.spacing(4))) -
                       (-xx[:-1,:-1]**3+6.*xx[:-1,:-1]*yy[1:,1:]*np.sqrt(xx[:-1,:-1]**2+yy[1:,1:]**2) + 
                        3.*yy[1:,1:]**3*np.log(np.sqrt(xx[:-1,:-1]**2+yy[1:,1:]**2) + xx[:-1,:-1] + np.spacing(4)) + 
                        3.*xx[:-1,:-1]**3*np.log(np.sqrt(xx[:-1,:-1]**2+yy[1:,1:]**2)+yy[1:,1:] + np.spacing(4))) +
                       (-xx[1:,1:]**3+6.*xx[1:,1:]*yy[1:,1:]*np.sqrt(xx[1:,1:]**2+yy[1:,1:]**2) + 
                        3.*yy[1:,1:]**3*np.log(np.sqrt(xx[1:,1:]**2+yy[1:,1:]**2) + xx[1:,1:] + np.spacing(4)) + 
                        3.*xx[1:,1:]**3*np.log(np.sqrt(xx[1:,1:]**2+yy[1:,1:]**2)+yy[1:,1:] + np.spacing(4))))

                self.inc = self._int_cone(pa,dpa,xedges,yedges)
                #H = H/self.inc
                H = np.rot90(H)
                H = np.flipud(H)
                return_hist_dict[species].append(H)
                return_hist_dict[species].append(xedges)
                return_hist_dict[species].append(yedges)


        return return_hist_dict

    def _vdist_pizza(self,
                     par=False,
                     z0=0.,
                     dz=1.,
                     wax=0,
                     v0_frame=False,
                     **kwargs): 
        """ Pizza Party!!!
            
            Also I am sooo sory if you ever have to
            use this code in this version... :(

            A little help:
            wax = which axis do we integarte over?
            wax = (0,1,2) -> (x,y,z) 
                or if in field line coordinates
            wax = (0,1,2) -> (||, ExB, bxExB)

        """

        self._z0  = z0
        self._dz  = dz
        self._wax = wax

        if par:
            b_interp = np.array([ 
                       self.interp_field(self.dump_field_dict['bx'],
                                         self.param_dict['lx'],
                                         self.param_dict['ly'],self._r0),
                       self.interp_field(self.dump_field_dict['by'],
                                         self.param_dict['lx'],
                                         self.param_dict['ly'],self._r0),
                       self.interp_field(self.dump_field_dict['bz'],
                                         self.param_dict['lx'],
                                         self.param_dict['ly'],self._r0)])

            e_interp = np.array([ 
                       self.interp_field(self.dump_field_dict['ex'],
                                         self.param_dict['lx'],
                                         self.param_dict['ly'],self._r0),
                       self.interp_field(self.dump_field_dict['ey'],
                                         self.param_dict['lx'],
                                         self.param_dict['ly'],self._r0),
                       self.interp_field(self.dump_field_dict['ez'],
                                         self.param_dict['lx'],
                                         self.param_dict['ly'],self._r0)])
            # Mike wanted me to do it this way
            fld = self.box_avg()

            #print 'Frist the interp fields:'
            #print 'B = ',b_interp
            #print 'E = ',e_interp
            b_interp = np.array([fld['bx'],fld['by'],fld['bz']])
            e_interp = np.array([fld['ex'],fld['ey'],fld['ez']])

            #print 'Then the avged fields:'
            #print 'B = ',b_interp
            #print 'E = ',e_interp

            exb = np.cross(e_interp,b_interp)/sum(b_interp**2)
            exb = exb/np.sqrt(sum(exb**2))

            b_interp = b_interp/np.sqrt(sum(b_interp**2))

            bxexb = np.cross(b_interp,exb) 
            #print b_interp
            #print exb
            #print bxexb

        velo={}
        return_hist_dict = {}
        self.subpart = {}
        #print 'Mucking this up!!'
        #Npart_2 = int(np.shape(self.particles['i']['vx'])[0]/2)

        #self.particles['i']['vx'][:Npart_2] = 1.*np.random.normal(scale=1.0,size=Npart_2)
        #self.particles['i']['vy'][:Npart_2] = 1.*np.random.normal(scale=1.0,size=Npart_2)
        #self.particles['i']['vz'][:Npart_2] = 1.*np.random.normal(scale=1.0,size=Npart_2)

        #self.particles['i']['vx'][Npart_2:] = 1.*np.random.normal(loc=10.,scale=1.0,size=Npart_2)
        #self.particles['i']['vy'][Npart_2:] = 1.*np.random.normal(loc=4.,scale=1.0,size=Npart_2)
        #self.particles['i']['vz'][Npart_2:] = 1.*np.random.normal(loc=4.,scale=1.0,size=Npart_2)
        ##print 'Done Mucking this up!!!'
        for species in self.species:
            return_hist_dict[species] = []
       
            if par:
# MuckMuckMuckMuck!!!!!
#                vxxx = self.particles[species]['vx'] + 1.0
# MuckMuckMuckMuck!!!!!

                #v0 = (b_interp[0]*vxxx+
                v0 = (b_interp[0]*self.particles[species]['vx']+
                      b_interp[1]*self.particles[species]['vy']+
                      b_interp[2]*self.particles[species]['vz'])

                #v1 = (exb[0]*vxxx+
                v1 = (exb[0]*self.particles[species]['vx']+
                      exb[1]*self.particles[species]['vy']+
                      exb[2]*self.particles[species]['vz'])

                #v2 = (bxexb[0]*vxxx+
                v2 = (bxexb[0]*self.particles[species]['vx']+
                      bxexb[1]*self.particles[species]['vy']+
                      bxexb[2]*self.particles[species]['vz'])
#Light Debuging
                self.bbb = b_interp
                self.exb = exb
                self.beb = bxexb
                #print 'Mucking this up!!!'
                #Npart = 100000
                #v0 = 1.*np.random.normal(scale=3.0,size=Npart)
                #v1 = 1.*np.random.normal(scale=1.0,size=Npart)
                #v2 = 1.*np.random.normal(scale=1.0,size=Npart)
            else:
                v0 = self.particles[species]['vx']
                v1 = self.particles[species]['vy']
                v2 = self.particles[species]['vz']

            if wax == 0: # wax is just which axis defines the plane we are looking at
                vp0 = v1
                vp1 = v2
                vax = v0
            elif wax == 1: 
                vp0 = v2
                vp1 = v0
                vax = v1
            elif wax == 2:
                vp0 = v0
                vp1 = v1
                vax = v2
            else :
                print ''
                print 'The plane axis is out of bounds fo wax = ',wax
                print 'vdist_2d is crashing!!!!'
                print ''

            if v0_frame:
                print '#'*80
                print 'Shifting vax by ',np.mean(vax)
                print '#'*80
                vax = vax - np.mean(vax)

            subpartind = np.where(abs(vax - z0) < dz/2.)
            print ''
            print 'The reduced number of {0} is {1}'.format(species,np.size(subpartind))
            print ''

            H, xedges, yedges = np.histogram2d(vp0[subpartind],vp1[subpartind],**kwargs)

            H = np.rot90(H)
            H = np.flipud(H)
            return_hist_dict[species].append(H)
            return_hist_dict[species].append(xedges)
            return_hist_dict[species].append(yedges)


        return return_hist_dict




    def _vdist_2d_par(self,**kwargs):
        """
        I dont know what is better
            To rotate all particles to a particular b feild
            or rotate all particles to their own b field?

            right now each part has its own b field
        """
        
        b_interp = np.array([ self.interp_field(self.dump_field_dict['bx'],self.param_dict['lx'],self.param_dict['ly'],self._r0),
                           self.interp_field(self.dump_field_dict['by'],self.param_dict['lx'],self.param_dict['ly'],self._r0),
                           self.interp_field(self.dump_field_dict['bz'],self.param_dict['lx'],self.param_dict['ly'],self._r0)])
        b_interp = b_interp/np.sqrt(sum(b_interp**2))

        e_interp = np.array([ self.interp_field(self.dump_field_dict['ex'],self.param_dict['lx'],self.param_dict['ly'],self._r0),
                           self.interp_field(self.dump_field_dict['ey'],self.param_dict['lx'],self.param_dict['ly'],self._r0),
                           self.interp_field(self.dump_field_dict['ez'],self.param_dict['lx'],self.param_dict['ly'],self._r0)])

        print 'b_interp'
        print b_interp
        self.trash1 = b_interp
        print 'e_interp'
        print e_interp
        self.trash2 = e_interp

        exb = np.cross(e_interp,b_interp) 
        exb = exb/np.sqrt(sum(exb**2))

        bxexb = np.cross(b_interp,exb) 

        if abs(sum(bxexb**2) - 1.0) > .001:
            print '###################### WARNING ######################'
            print '      your b and exb seem to not be perpandicular!!!'
            bxexb = bxexb/np.sqrt(sum(bxexb**2))

        #qc#print 'Rotating Velocties'
        velo={}
        for species in self.species:
            velo[species] = {}

            rotate_all_parts = False
            if rotate_all_parts == True:

                delx = self.param_dict['lx']*1.0/self.param_dict['pex']/self.param_dict['nx']
                dely = self.param_dict['ly']*1.0/(self.param_dict['pey']*self.param_dict['ny'])

                xind = (np.floor((self.particles[species]['x']-delx/2.0)/delx)).astype(int)
                yind = (np.floor((self.particles[species]['y']-dely/2.0)/dely)).astype(int)

                wx = (self.particles[species]['x']-delx/2.0)%delx
                wy = (self.particles[species]['y']-dely/2.0)%dely

                partbx = wx     *wy     *self.dump_field_dict['bx'][xind.tolist()    ,yind.tolist()] + \
                         (1.-wx)*wy     *self.dump_field_dict['bx'][(xind+1).tolist(),yind.tolist()] + \
                         wx     *(1.-wy)*self.dump_field_dict['bx'][xind.tolist()    ,(yind+1).tolist()] + \
                         (1.-wx)*(1.-wy)*self.dump_field_dict['bx'][(xind+1).tolist(),(yind+1).tolist()] 

                partby = wx     *wy     *self.dump_field_dict['by'][xind.tolist()    ,yind.tolist()] + \
                         (1.-wx)*wy     *self.dump_field_dict['by'][(xind+1).tolist(),yind.tolist()] + \
                         wx     *(1.-wy)*self.dump_field_dict['by'][xind.tolist()    ,(yind+1).tolist()] + \
                         (1.-wx)*(1.-wy)*self.dump_field_dict['by'][(xind+1).tolist(),(yind+1).tolist()] 

                partbz = wx     *wy     *self.dump_field_dict['bz'][xind.tolist()    ,yind.tolist()] + \
                         (1.-wx)*wy     *self.dump_field_dict['bz'][(xind+1).tolist(),yind.tolist()] + \
                         wx     *(1.-wy)*self.dump_field_dict['bz'][xind.tolist()    ,(yind+1).tolist()] + \
                         (1.-wx)*(1.-wy)*self.dump_field_dict['bz'][(xind+1).tolist(),(yind+1).tolist()] 

                (partbx,partby,partbz) = (partbx/np.sqrt(partbx**2+partby**2+partbz**2),
                                          partby/np.sqrt(partbx**2+partby**2+partbz**2),
                                          partbz/np.sqrt(partbx**2+partby**2+partbz**2))

                partpb1 =  0.0*partbx
                partpb2 = -1.*np.sign(by_interp)*partbz/(partbz**2 + partby**2)**(.5)
                partpb3 = 1.*np.sign(by_interp)*partby/(partbx**2 + partby**2)**(.5)

                (partpb1,partpb2,partpb3) = (partpb1/np.sqrt(partpb1**2+partpb2**2+partpb3**2),
                                             partpb2/np.sqrt(partpb1**2+partpb2**2+partpb3**2),
                                             partpb3/np.sqrt(partpb1**2+partpb2**2+partpb3**2))

                partpp1 =  (partby*partpb3 - partbz*partpb2)
                partpp2 =  (partbz*partpb1 - partbx*partpb3)
                partpp3 =  (partbx*partpb2 - partby*partpb1)

                velo[species]['par']   = (partbx*self.particles[species]['vx']+partby*self.particles[species]['vy']+partbz*self.particles[species]['vz'])
                velo[species]['perp1'] = (partpb1*self.particles[species]['vx']+partpb2*self.particles[species]['vy']+partpb3*self.particles[species]['vz'])
                velo[species]['perp2'] = (partpp1*self.particles[species]['vx']+partpp2*self.particles[species]['vy']+partpp3*self.particles[species]['vz'])

            else:
                velo[species]['par']   = (b_interp[0]*self.particles[species]['vx']+
                                          b_interp[1]*self.particles[species]['vy']+
                                          b_interp[2]*self.particles[species]['vz'])

                velo[species]['perp1'] = (exb[0]*self.particles[species]['vx']+
                                          exb[1]*self.particles[species]['vy']+
                                          exb[2]*self.particles[species]['vz'])

                velo[species]['perp2'] = (bxexb[0]*self.particles[species]['vx']+
                                          bxexb[1]*self.particles[species]['vy']+
                                          bxexb[2]*self.particles[species]['vz'])


        #velo['e']['par']   = (bx_interp*self.particles['e']['vx']+by_interp*self.particles['e']['vy']+bz_interp*self.particles['e']['vz'])/bmag_interp
        #velo['e']['perp1'] = self.particles['e']['vx']*b_perp1x+self.particles['e']['vy']*b_perp1y+self.particles['e']['vz']*b_perp1z
        #velo['e']['perp2'] = self.particles['e']['vx']*b_perp2x+self.particles['e']['vy']*b_perp2y+self.particles['e']['vz']*b_perp2z

        #qc#print 'Generating Histograms'

        return_hist_dict = {}
        for species in self.species:
            return_hist_dict[species] = []

        for species in velo.keys():
            H, xedges, yedges = np.histogram2d(velo[species]['par'],velo[species]['perp1'],**kwargs)
# H needs to be rotated and flipped
            H = np.rot90(H)
            H = np.flipud(H)
            return_hist_dict[species].append(['V_b vs V_exb',H,xedges,yedges])

            H, xedges, yedges = np.histogram2d(velo[species]['par'],velo[species]['perp2'],**kwargs)
            H = np.rot90(H)
            H = np.flipud(H)
            return_hist_dict[species].append(['V_b vs V_bxexb 2',H,xedges,yedges])

            H, xedges, yedges = np.histogram2d(velo[species]['perp1'],velo[species]['perp2'],**kwargs)
            H = np.rot90(H)
            H = np.flipud(H)
            return_hist_dict[species].append(['V_exb 1 vs V_bxexb',H,xedges,yedges])

# Mask zeros
         #Hmasked = np.ma.masked_where(H==0,H)

        return return_hist_dict


    def _vdist_2d(self,**kwargs):
        """
        """
        #qc#print 'Generating Histograms'

        return_hist_dict = {}
        return_hist_dict['i'] = []
        return_hist_dict['e'] = []
        if kwargs.has_key('gen_vec'):
            gvec = kwargs.pop('gen_vec')
            vec1 = np.array(gvec[0])
            vec2 = np.array(gvec[1])
            vec1 = vec1/np.sqrt(np.sum(vec1**2))
            vec2 = vec2/np.sqrt(np.sum(vec2**2))
            if abs(np.dot(vec1,vec2)) > 10e-4:
                print ''
                print 'The two vectors given are not orthognal!!!!'
                print 'The exception for handeling this has not been coded'
                print 'kindly fix your shit and try again ~Colby'
                print ''
                return -1
            else:
                vec3 = np.cross(vec1,vec2)

            velo = {}
            for species in self.species:

                velo[species] = {}

                velo[species]['v1']   = (vec1[0]*self.particles[species]['vx']+
                                         vec1[1]*self.particles[species]['vy']+
                                         vec1[2]*self.particles[species]['vz'])

                velo[species]['v2'] = (vec2[0]*self.particles[species]['vx']+
                                       vec2[1]*self.particles[species]['vy']+
                                       vec2[2]*self.particles[species]['vz'])

                velo[species]['v3'] = (vec3[0]*self.particles[species]['vx']+
                                       vec3[1]*self.particles[species]['vy']+
                                       vec3[2]*self.particles[species]['vz'])

                H, xedges, yedges = np.histogram2d(velo[species]['v1'],velo[species]['v2'],**kwargs)
# H needs to be rotated and flipped
                H = np.rot90(H)
                H = np.flipud(H)
                return_hist_dict[species].append(['V_1 vs V_2',H,xedges,yedges])

                H, xedges, yedges = np.histogram2d(velo[species]['v1'],velo[species]['v3'],**kwargs)
                H = np.rot90(H)
                H = np.flipud(H)
                return_hist_dict[species].append(['V_1 vs V_3',H,xedges,yedges])

                H, xedges, yedges = np.histogram2d(velo[species]['v2'],velo[species]['v3'],**kwargs)
                H = np.rot90(H)
                H = np.flipud(H)
                return_hist_dict[species].append(['V_2 vs V_3',H,xedges,yedges])
        else:
            for species in self.species:
                H, xedges, yedges = np.histogram2d(self.particles[species]['vx'],self.particles[species]['vy'],**kwargs)
# H needs to be rotated and flipped
                H = np.rot90(H)
                H = np.flipud(H)
                return_hist_dict[species].append(['V_X vs V_Y',H,xedges,yedges])

                H, xedges, yedges = np.histogram2d(self.particles[species]['vx'],self.particles[species]['vz'],**kwargs)
                H = np.rot90(H)
                H = np.flipud(H)
                return_hist_dict[species].append(['V_X vs V_Z',H,xedges,yedges])

                H, xedges, yedges = np.histogram2d(self.particles[species]['vy'],self.particles[species]['vz'],**kwargs)
                H = np.rot90(H)
                H = np.flipud(H)
                return_hist_dict[species].append(['V_Y vs V_Z',H,xedges,yedges])

# Mask zeros
         #Hmasked = np.ma.masked_where(H==0,H)

        return return_hist_dict

    #def get_part_in_box([location, width]):
    def _part_in_box(self,r0=[0.5,0.5],dx=[1.,1.]):
        """
        #--------------------------------------------------------------
        #   Method      : get_part_in_box
        #
        #   Discription : This method accepts a point and a width in 
        #                 simulation units (c/wpi) to define a box.
        #               : In that box we bin all of the particles to 
                          form the effective distrobution function
        #
        #   Args     r0 : location [x,y] ( where you want the center 
                          of you box to be located at)
        #            dx : width [x,y] (the width of the box to bin 
                           particles 
        #               : dump_num (This spesifies the particular 
                          runs dump file 
        #
        #   Comments    : It would be pretty easy and potential 
                          usefull to allow this to wrap around 
                          the edges
        #               : so we can pick zero as a boundry and the 
                          code will know what todo.
                          #--------------------------------------------------------------
        """
        x0 = r0[0]
        y0 = r0[1]
        if isinstance(dx,float) or isinstance(dx,float):
            dx = dx*1.0
            dy = dx       # Square box is assumed
        else:
            dy = dx[1]
            dx = dx[0]

        #qc#print 'r0 = [%f,%f] and dx = [%f,%f]'%(x0,y0,dx,dy)
# Figure out which set of processors we are on
        xlb = x0 - dx/2.
        xub = x0 + dx/2.
        ylb = y0 - dy/2.
        yub = y0 + dy/2.

# BIGNOTE: This seems iffy you should come back and double check this colby. It doesnt make sence that we should have to add 1
        #xproc_lb = (int(np.floor(1.0*self.param_dict['pex']*xlb/self.param_dict['lx']))+1)%int(self.param_dict['nchannels'])
        #xproc_lb = (int(np.floor(1.0*self.param_dict['pex']*xlb/self.param_dict['lx']))+1)%int(self.param_dict['nchannels'])
        if type(self.param_dict['nchannels']) == str:
# Man you should fix this colby!!!
            #self.param_dict['nchannels'] =  self.param_dict['pex'] # a lot of times we just set nchannels as pex
            self.param_dict['nchannels'] =  self.param_dict[self.param_dict['nchannels']]
            
        
        #xproc_lb = (int(np.floor(1.0*self.param_dict['pex']*xlb/self.param_dict['lx']))+1)%self.param_dict['nchannels']
# We are chaning this but it may not work for all stuff so you know 
        xproc_lb = (int(np.floor(1.0*self.param_dict['pex']*xlb/self.param_dict['lx'])))%self.param_dict['nchannels'] +1
# Colby this needs to be cooded better!!!
# if you have a diffent number of channels than pex you run into some shit
        yproc_lb = int(np.floor(1.0*self.param_dict['pey']*ylb/self.param_dict['ly']))
        #yproc_lb = int(np.floor(1.0*self.param_dict['pey']*ylb/self.param_dict['ly']))*2 #Uncomment for yishin
        #xproc_ub = (int(np.floor(1.0*self.param_dict['pex']*xub/self.param_dict['lx']))+1)%int(self.param_dict['nchannels'])


        if abs(xub - self.param_dict['lx']) < abs(np.spacing(2)): 
            xub = self.param_dict['lx'] - np.spacing(2)
        xproc_ub = (int(np.floor(1.0*self.param_dict['pex']*xub/self.param_dict['lx'])))%self.param_dict['nchannels'] +1
        yproc_ub = int(np.floor(1.0*self.param_dict['pey']*yub/self.param_dict['ly'])) 
        #yproc_ub = int(np.floor(1.0*self.param_dict['pey']*yub/self.param_dict['ly']))*2 #Uncomment for yishin

        if xproc_lb > xproc_ub:
            print 'Lower Bound greater than upper bound! That is obviously an issue!'
            return -1

        if xproc_lb < 1: xproc_lb = 1
        if xproc_ub > self.param_dict['nchannels']: xproc_ub = self.param_dict['nchannels'] 

#Colby! you can code this smarter but presently you are not doing that!
# To code smarter, insted of making this simple upper bound you should code
# to allow for nchanles to be a non multiple of pex
        max_yproc =  int(round(self.param_dict['pex']*self.param_dict['pey']/self.param_dict['nchannels']))
        if yproc_lb < 0: yproc_lb = 0
        if yproc_ub > max_yproc:
            yproc_ub = max_yproc


        xprocs_2load = range(xproc_lb,xproc_ub+1)
        yprocs_2load = range(yproc_lb,yproc_ub+1)
        
        print 'The x processors we will be loading (i.e. the dump files) are: {}'.format(xprocs_2load)
        print 'The y processors we will be loading (i.e. the sub arrays) are: {}'.format(yprocs_2load)
        #c# print 'That means we have {} processors to load and you can expect an aprox. {:.2f} min wait time'.format(
        #c#      len(xprocs_2load)*len(yprocs_2load),1.41*len(xprocs_2load)*len(yprocs_2load)*4./60.)

# Load in the appropriate Processors
        temp_dump_pruned = {}# List to hold all the particles
        for species in self.species:
            temp_dump_pruned[species] = []
        #first for loop over px
        for xprocs_index in xprocs_2load:
            dump_dat_dict = {}
            dump_index = self._num_to_ext(xprocs_index)
            #qc#print 'Loading in xproc number '+dump_index
            temp_dump_dat = self.read_dump_file(dump_index)

            #c# new_tdd = {}
# We need to#c#  throw out the extra data
            #c# if (int(np.floor(1.0*self.param_dict['pex']*xlb/self.param_dict['lx']))+1)/self.param_dict['nchannels'] < 1:
            #c#     new_tdd['i'] = temp_dump_dat['i'][:len(temp_dump_dat['i'])/2-1]
            #c#     new_tdd['e'] = temp_dump_dat['e'][:len(temp_dump_dat['e'])/2-1]
            #c# else:
            #c#     new_tdd['i'] = temp_dump_dat['i'][len(temp_dump_dat['i'])/2:]
            #c#     new_tdd['e'] = temp_dump_dat['e'][len(temp_dump_dat['e'])/2:]

            #c# temp_dump_dat = new_tdd


# We need to looop over Ions and Electrons
            for species in self.species:
                #if species == 'i':
                #    print '\tSelecting Ions'
                #else:
                #    print '\tSelecting Electrons'
# second for loop over the py
# also just doing electron for right now
                for yprocs_index in yprocs_2load:
                    self._debug = temp_dump_dat 
                    temp_dump_yproc = temp_dump_dat[species][yprocs_index]
### Lets try somthing new, that might be faster.  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    temp_index = np.where((temp_dump_dat[species][yprocs_index]['y'] - yub)**2 +
                    (temp_dump_dat[species][yprocs_index]['y'] - ylb)**2 < (yub-ylb)**2 )
                    temp_dump_xproc = temp_dump_dat[species][yprocs_index][temp_index]

                    temp_index = np.where((temp_dump_xproc['x'] - xub)**2 +
                    (temp_dump_xproc['x'] - xlb)**2 < (xub-xlb)**2 )

                    temp_dump_pruned[species].append(temp_dump_xproc[temp_index])


#test ### This is how we were sorting this But I think it is MUCH SLOWER
#test # You only need to sort if you are on the edge processors
#test                     if yprocs_index == yprocs_2load[0] or yprocs_index == yprocs_2load[-1]: 
#test                         #qc#print '\t\tSorting yproc number '+str(yprocs_index)
#test                         sorted_index = temp_dump_dat[species][yprocs_index].argsort(order='y')
#test                         temp_dump_yproc = temp_dump_dat[species][yprocs_index][sorted_index]
#test # Here we need kind of a complecated if structure to get all the poible cases since
#test # we are scaning over muliple processors.
#test # If you are on your first y processor then you need to find a lower boundry
#test                     if yprocs_index == yprocs_2load[0]: 
#test                         #qc#print '\t\t\tFinding lower yboundry index '
#test                         lower_yboundry_index = np.searchsorted(temp_dump_yproc['y'],ylb)
#test                     else:
#test                         lower_yboundry_index = 0#np.searchsorted(temp_dump_yproc['y'],ylb)
#test # If you are on your last y processor then you need to find a upper boundry
#test                     if yprocs_index == yprocs_2load[-1]: 
#test                         #qc#print '\t\t\tFinding upper yboundry index '
#test                         upper_yboundry_index = np.searchsorted(temp_dump_yproc['y'],yub)
#test                     else:
#test                         upper_yboundry_index = -1#np.searchsorted(temp_dump_yproc['y'],yub)
#test                     # You only need to sort if you are on the edge processors
#test                     temp_dump_xproc = temp_dump_yproc[lower_yboundry_index:upper_yboundry_index]
#test                     if xprocs_index == xprocs_2load[0] or xprocs_index == xprocs_2load[-1]: 
#test                         #qc#print '\t\tNow sorting x values for remaing data'
#test                         sorted_index = temp_dump_xproc.argsort(order='x')
#test                         temp_dump_xproc = temp_dump_xproc[sorted_index] 
#test # If you are on your first x processor then you need to find a lower boundry
#test                     if xprocs_index == xprocs_2load[0]: 
#test                         #qc#print '\t\t\tFinding lower xboundry index '
#test                         lower_xboundry_index = np.searchsorted(temp_dump_xproc['x'],xlb)
#test                     else:
#test                         lower_xboundry_index = 0#np.searchsorted(temp_dump_xproc['x'],xlb)
#test # If you are on your last x processor then you need to find a upper boundry
#test                     if xprocs_index == xprocs_2load[-1]: 
#test                         #qc#print '\t\t\tFinding upper xboundry index '
#test                         upper_xboundry_index = np.searchsorted(temp_dump_xproc['x'],xub)
#test                     else: 
#test                         upper_xboundry_index = -1#np.searchsorted(temp_dump_xproc['x'],xub)
#test                     temp_dump_pruned[species].append(temp_dump_xproc[lower_xboundry_index:upper_xboundry_index])
        for species in self.species:
            temp_dump_pruned[species] = np.concatenate(temp_dump_pruned[species])
            print 'Total %s in box are: %i'%(species,len(temp_dump_pruned[species]))
        return temp_dump_pruned



    def _num_to_ext(self,num):
        if type(num) is str: 
            return num
        else:
            if int(np.floor(num / 10)) > 0: 
                if int(np.floor(num / 100)) > 0: 
                    return str(num)
                else:
                    return '0'+str(num)
            else:
                return '00'+str(num)

    def _int_cone(self,pa,dpa,xedges,yedges):
        """
#---------------------------------------------------------------------------------------------------------------
#   Method      : _int_cone
#
#   Discription : This integrates a cone over a differental grid
#
#   Args        : xedges the x edges of the histogram
#               : yedges the y edges of the histogram
#
#   Comments    : I think this is working ok? 
#---------------------------------------------------------------------------------------------------------------
        """
        

        xx,yy = np.meshgrid(yedges,xedges)
        intcone = 1./18.*(
               (-xx[:-1,:-1]**3+6.*xx[:-1,:-1]*yy[:-1,:-1]*np.sqrt(xx[:-1,:-1]**2+yy[:-1,:-1]**2) + 
                3.*yy[:-1,:-1]**3*np.log(np.sqrt(xx[:-1,:-1]**2+yy[:-1,:-1]**2) + xx[:-1,:-1] + np.spacing(4)) + 
                3.*xx[:-1,:-1]**3*np.log(np.sqrt(xx[:-1,:-1]**2+yy[:-1,:-1]**2)+yy[:-1,:-1] + np.spacing(4))) -
               (-xx[1:,1:]**3+6.*xx[1:,1:]*yy[:-1,:-1]*np.sqrt(xx[1:,1:]**2+yy[:-1,:-1]**2) + 
                3.*yy[:-1,:-1]**3*np.log(np.sqrt(xx[1:,1:]**2+yy[:-1,:-1]**2) + xx[1:,1:] + np.spacing(4)) + 
                3.*xx[1:,1:]**3*np.log(np.sqrt(xx[1:,1:]**2+yy[:-1,:-1]**2)+yy[:-1,:-1] + np.spacing(4))) -
               (-xx[:-1,:-1]**3+6.*xx[:-1,:-1]*yy[1:,1:]*np.sqrt(xx[:-1,:-1]**2+yy[1:,1:]**2) + 
                3.*yy[1:,1:]**3*np.log(np.sqrt(xx[:-1,:-1]**2+yy[1:,1:]**2) + xx[:-1,:-1] + np.spacing(4)) + 
                3.*xx[:-1,:-1]**3*np.log(np.sqrt(xx[:-1,:-1]**2+yy[1:,1:]**2)+yy[1:,1:] + np.spacing(4))) +
               (-xx[1:,1:]**3+6.*xx[1:,1:]*yy[1:,1:]*np.sqrt(xx[1:,1:]**2+yy[1:,1:]**2) + 
                3.*yy[1:,1:]**3*np.log(np.sqrt(xx[1:,1:]**2+yy[1:,1:]**2) + xx[1:,1:] + np.spacing(4)) + 
                3.*xx[1:,1:]**3*np.log(np.sqrt(xx[1:,1:]**2+yy[1:,1:]**2)+yy[1:,1:] + np.spacing(4))))

        norm = 1./18.*(
               (-xx[0,0]**3+6.*xx[0,0]*yy[0,0]*np.sqrt(xx[0,0]**2+yy[0,0]**2) + 
                3.*yy[0,0]**3*np.log(np.sqrt(xx[0,0]**2+yy[0,0]**2) + xx[0,0] + np.spacing(4)) + 
                3.*xx[0,0]**3*np.log(np.sqrt(xx[0,0]**2+yy[0,0]**2)+yy[0,0] + np.spacing(4))) -
               (-xx[-1,-1]**3+6.*xx[-1,-1]*yy[0,0]*np.sqrt(xx[-1,-1]**2+yy[0,0]**2) + 
                3.*yy[0,0]**3*np.log(np.sqrt(xx[-1,-1]**2+yy[0,0]**2) + xx[-1,-1] + np.spacing(4)) + 
                3.*xx[-1,-1]**3*np.log(np.sqrt(xx[-1,-1]**2+yy[0,0]**2)+yy[0,0] + np.spacing(4))) -
               (-xx[0,0]**3+6.*xx[0,0]*yy[-1,-1]*np.sqrt(xx[0,0]**2+yy[-1,-1]**2) + 
                3.*yy[-1,-1]**3*np.log(np.sqrt(xx[0,0]**2+yy[-1,-1]**2) + xx[0,0] + np.spacing(4)) + 
                3.*xx[0,0]**3*np.log(np.sqrt(xx[0,0]**2+yy[-1,-1]**2)+yy[-1,-1] + np.spacing(4))) +
               (-xx[-1,-1]**3+6.*xx[-1,-1]*yy[-1,-1]*np.sqrt(xx[-1,-1]**2+yy[-1,-1]**2) + 
                3.*yy[-1,-1]**3*np.log(np.sqrt(xx[-1,-1]**2+yy[-1,-1]**2) + xx[-1,-1] + np.spacing(4)) + 
                3.*xx[-1,-1]**3*np.log(np.sqrt(xx[-1,-1]**2+yy[-1,-1]**2)+yy[-1,-1] + np.spacing(4))))
        
        self.normie=norm

        print 'norm = ',norm
        print 'otha = ',abs(1./np.tan((pa-dpa/2.)/180.*np.pi) - 1./np.tan((pa+dpa/2.)/180.*np.pi))

        #return intcone/norm#*abs(1./np.tan((pa-dpa/2.)/180.*np.pi) - 1./np.tan((pa+dpa/2.)/180.*np.pi))
        return intcone/intcone.min()#*abs(1./np.tan((pa-dpa/2.)/180.*np.pi) - 1./np.tan((pa+dpa/2.)/180.*np.pi))

# Ripped this code from else where nee to add it
    def box_avg(self):
        
        if not hasattr(self, 'dump_field_dict'):
            print 'Field data not loaded somthings wrong!!!!'
            return None
        
        r0 = self._r0
        dx = self._dx

        x0 = r0[0] - dx[0]/2.
        x1 = r0[0] + dx[0]/2.
        y0 = r0[1] - dx[1]/2.
        y1 = r0[1] + dx[1]/2.

        dx0 = 1.0*self.param_dict['lx']/(self.param_dict['nx']*self.param_dict['pex'])

        ip0 = int(np.ceil((x0 -  dx0/2.)/dx0))
        ip1 = int(np.floor((x1 - dx0/2.)/dx0))
        jp0 = int(np.ceil((y0 -  dx0/2.)/dx0))
        jp1 = int(np.floor((y1 - dx0/2.)/dx0))

        if ip1 < ip0 or jp1 < jp0:
            print 'No grid points found in box!!!!'
            return None
        #else:
        #    print 'DEBUG: ',ip0,ip1,jp0,jp1
       
        iprg, jprg = [],[]
        for c in range(ip0,ip1+1):
            for d in range(jp0,jp1+1):
                iprg.append(c)
                jprg.append(d)

        avg_flds = {}
        for k,fld in self.dump_field_dict.iteritems():
            avg_flds[k] = np.mean(fld[jprg,iprg])

        return avg_flds


    def _shift_frame(self):
        """ Shifts all the particles velocties into their fluid frame.
            This finds the mean ux uy and uz for e and i and shifts everying
            by thoes values
        """

        for s in self.species:
            for c in ['vx','vy','vz']:
                uu = np.mean(self.particles[s][c])
                print '{0}_{1} = {2}'.format(c,s,uu)
                self.particles[s][c] = self.particles[s][c] - uu


    def interp_field(self,field,lx=None,ly=None,r0=None):
        """
#----------------------------------------------------------------------
#        
#   Method      : interp_field
#
#   Discription : This method takes a field and a floating point, and returns the linear fit value 
#               : between the grid points
#
#   Args        : field  The field you are interpolating
#               : r0[0] The xpoint to interpolate at
#               : r0[1] The ypoint to interpolate at
#
#   Comments    : I think this is working ok? It would be smart to make
#                 this an object method that just reads
#               : the internally saved field. so CODE IN THE FUTURE
#----------------------------------------------------------------------
        """
        if lx is None:lx=self.param_dict['lx']
        if ly is None:ly=self.param_dict['ly']
        if r0 is None:r0=self._r0

        nx = len(field[0,:])
        ny = len(field[:,0])
        ip = int(np.floor(1.0*r0[0]/lx*nx))
        jp = int(np.floor(1.0*r0[1]/ly*ny))

        if ip + 1 > nx-1: ipp1 = 0
        else: ipp1 = ip+1

        if jp + 1 > ny-1: jpp1 = 0
        else: jpp1 = jp+1

       # print ip,jp,nx,ny

        wx = 1.0*r0[0]/lx*nx - np.floor(1.0*r0[0]/lx*nx)
        wy = 1.0*r0[1]/ly*ny - np.floor(1.0*r0[1]/ly*ny)

        return (1.-wx)*(1.-wy)*field[jp  , ip]+\
               (wx)   *(1.-wy)*field[jpp1, ip]+\
               (1.-wx)*   (wy)*field[jp  ,ipp1]+\
               (wx)   *   (wy)*field[jpp1,ipp1]
    


#Orphaned method that I use in load param but not quite sure how to fit it in
def convert(val):
    constructors = [int, float, str]
    for c in constructors:
        try:
            return c(val)
        except ValueError:
            pass
            

#c# This was a set of code I used to rotate each particle indviduly
#c# Mike insists this is silly, so we will not be using this
#c#            rotate_all_parts = False
#c#            if rotate_all_parts == True:
#c#
#c#                xind = (np.floor((self.particles[species]['x']-delx/2.0)/delx)).astype(int)
#c#                yind = (np.floor((self.particles[species]['y']-dely/2.0)/dely)).astype(int)
#c#
#c#                wx = (self.particles[species]['x']-delx/2.0)%delx
#c#                wy = (self.particles[species]['y']-dely/2.0)%dely
#c#
#c#                (1.-wx)*(1.-wy)
#c#                (1.-wx)*wy     
#c#                wx     *(1.-wy)
#c#                wx     *wy     
#c#
#c#                partbx = (1.-wx)*(1.-wy)*self.dump_field_dict['bx'][xind.tolist()    ,yind.tolist()] + \
#c#                         (1.-wx)*wy     *self.dump_field_dict['bx'][(xind+1).tolist(),yind.tolist()] + \
#c#                         wx     *(1.-wy)*self.dump_field_dict['bx'][xind.tolist()    ,(yind+1).tolist()] + \
#c#                         wx     *wy     *self.dump_field_dict['bx'][(xind+1).tolist(),(yind+1).tolist()] 
#c#
#c#                partby = (1.-wx)*(1.-wy)*self.dump_field_dict['by'][xind.tolist()    ,yind.tolist()] + \
#c#                         (1.-wx)*wy     *self.dump_field_dict['by'][(xind+1).tolist(),yind.tolist()] + \
#c#                         wx     *(1.-wy)*self.dump_field_dict['by'][xind.tolist()    ,(yind+1).tolist()] + \
#c#                         wx     *wy     *self.dump_field_dict['by'][(xind+1).tolist(),(yind+1).tolist()] 
#c#
#c#                partbz = (1.-wx)*(1.-wy)*self.dump_field_dict['bz'][xind.tolist()    ,yind.tolist()] + \
#c#                         (1.-wx)*wy     *self.dump_field_dict['bz'][(xind+1).tolist(),yind.tolist()] + \
#c#                         wx     *(1.-wy)*self.dump_field_dict['bz'][xind.tolist()    ,(yind+1).tolist()] + \
#c#                         wx     *wy     *self.dump_field_dict['bz'][(xind+1).tolist(),(yind+1).tolist()] 
#c#
#c#                vmag = np.sqrt(self.particles[species]['vx']**2+self.particles[species]['vy']**2+self.particles[species]['vz']**2)
#c#                vpar = (self.particles[species]['vx']*partbx+self.particles[species]['vy']*partby+self.particles[species]['vz']*partbz)/np.sqrt(partbx**2+partby**2+partbz**2)



        

