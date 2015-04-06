########################################################################################################################################################
########################################################################################################################################################
#                                                                                                                                                      #
#                                                 Python Progs :  p3d_dump.py                                                                         #
#                                                 Aruthor      :  Colby Haggerty                                                                       #
#                                                 Date         :  2014.08.15                                                                           #
#                                                                                                                                                      #
#                                                                                                                                                      #
########################################################################################################################################################
### Discription:
#
#       p3d_movie class is a movie object to be handled by p3d_run. This object keeps all of the information for a  
#       given movie. We have seperated this from the p3d_object to keep things neatter and cleaner. The plan is then
#       apply this same idea to the dump files. NOTE there may be a time trade off issue between this and doing every
#       thing together, but i dont think that will be an issue.
#

import os
import sys 
import datetime
import numpy as np
import struct
import glob
from scipy.io.idl import readsav

class p3d_dump(object):
    """p3d_run object """

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



    def read_dump_file(self,dump_index=-1,fields=False,verbose=False):
        """
#---------------------------------------------------------------------------------------------------------------
#   Method      : read_dump_file
# 
#   Args        : dump_index (is the middle 3 digets of the dump file. This corseponds to the processors)
#               : dump_num (is the last 3 digets of the dump file. This spesifies which run)
# 
#   BIG NOTE    : This assumes that we are looking at 2 D dump files! so file 1 is field file 
# 
#   Comment     : 
#---------------------------------------------------------------------------------------------------------------
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

    def vdist_2d(self,r0=[0.5,0.5],dx=[1.0,1.0],species='e',par=False,Bvec=False,OneD=False,**kwargs):
        """
        #---------------------------------------------------------------------------------------------------------------
        #   Method      : vdist_2d_par
        #
        #   Discription : This method accepts a point and a width in simulation units (c/wpi) to define a box.
        #               : In that box we bin all of the particles to form the effective distrobution function
        #
        #   Args        : location [x,y] ( where you want the center of you box to be located at)
        #               : width [x,y] (the width of the box to bin particles 
        #               : dump_num (This spesifies the particular runs dump file 
        #
        #   Returns     : return_hist_dict['axis'][i/e][hist_array]
        #               :   'axis'  slecets the data for the 2d projection of the distabution function 
        #                   'axis' = 'parperp1', 'parperp1', 'perp1perp2'
        #                   i/e is an int that selects either the electron or ion data. i = 0, e = 1
        #                   The rest is just an object that is returned directly from the histogram2d method
        #
        #
        #   Comments    : It would be pretty easy and potential usefull to allow this to wrap around the edges
        #               : so we can pick zero as a boundry and the code will know what todo.
        #---------------------------------------------------------------------------------------------------------------
        """

# Turn this into a method
#%%%%%%%%%%%%%%%
        if OneD: # Make dy = ly
            r0[1] = self.param_dict['ly']/2.0
            dx = [dx[0], self.param_dict['ly']/2.0]
        if not hasattr(self,'_r0') or not hasattr(self,'particles'):
            self._r0 = r0
            self._dx = dx
# Calling get particles in box to make the vdist
            #qc#print 'Reading Ions and Electrons from the Dump File'
            self.particles = self._part_in_box(r0,dx)
        else:
            if (self._r0[0] == r0[0]) and (self._r0[1] == r0[1]) and (self._dx == dx):
                print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
                print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Same r0 and dx found, using old particles ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
                print 'r0 = '+str(r0)+' _r0 = '+str(self._r0)
                print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
                print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
            else:
                self._r0 = r0
                self._dx = dx
                #qc#print 'Reading Ions and Electrons from the Dump File'
                self.particles = self._part_in_box(r0,dx)
#%%%%%%%%%%%%%%%

# Reading in fields to calculate vpar
# Now there are two ways to do this, a faster and a slower
# Faster: just take the average B field and and use that
#         for every particle
# Slower: Use an interpolated Bfield for each point
#
# Right now im only coding the faster one
        if not kwargs.has_key('bins'): kwargs['bins']=51
        if par or Bvec:
            #qc#print 'Reading in the Fields form the Dump File'
            self.dump_field_dict = self.read_dump_file(fields=True)
            if par:
                return_hist = self._vdist_2d_par(**kwargs)
            else:
                return_hist = self._vdist_2d(**kwargs)
                #qc#print 'Interpolating the Bfield at the given r0 value'
                return_hist['B'] = np.array([self.interp_field(self.dump_field_dict['bx']), \
                                             self.interp_field(self.dump_field_dict['by']), \
                                             self.interp_field(self.dump_field_dict['bz'])])
        else:
            return_hist = self._vdist_2d(**kwargs)

        return return_hist


    def _vdist_2d_par(self,**kwargs):
        #qc#print 'Calculating B field'
        # this version doent work 
        #bx_interp = self.interp_field(self.dump_field_dict['bx'],self.param_dict['lx'],self.param_dict['ly'],self._r0[0],self._r0[1])
        #by_interp = self.interp_field(self.dump_field_dict['by'],self.param_dict['lx'],self.param_dict['ly'],self._r0[0],self._r0[1])
        #bz_interp = self.interp_field(self.dump_field_dict['bz'],self.param_dict['lx'],self.param_dict['ly'],self._r0[0],self._r0[1])
        bx_interp = self.interp_field(self.dump_field_dict['bx'],self.param_dict['lx'],self.param_dict['ly'],self._r0)
        by_interp = self.interp_field(self.dump_field_dict['by'],self.param_dict['lx'],self.param_dict['ly'],self._r0)
        bz_interp = self.interp_field(self.dump_field_dict['bz'],self.param_dict['lx'],self.param_dict['ly'],self._r0)
        bmag_interp = (bx_interp**2+by_interp**2+bz_interp**2)**.5

        if by_interp > 0.:
            b_perp1x = 0.
            b_perp1y = -1.*bz_interp/(bz_interp**2 + by_interp**2)**(.5)
            b_perp1z = by_interp/(bx_interp**2 + by_interp**2)**(.5)
        else:
            b_perp1x = 0.
            b_perp1y = bz_interp/(bz_interp**2 + by_interp**2)**(.5)
            b_perp1z = -1.*by_interp/(bx_interp**2 + by_interp**2)**(.5)

        b_perpmag = (b_perp1x**2+b_perp1y**2+b_perp1z**2)**.5
        b_perp1x = b_perp1x/b_perpmag
        b_perp1y = b_perp1y/b_perpmag
        b_perp1z = b_perp1z/b_perpmag

        b_perp2x = (by_interp*b_perp1z - bz_interp*b_perp1y)
        b_perp2y = (bz_interp*b_perp1x - bx_interp*b_perp1z)
        b_perp2z = (bx_interp*b_perp1y - by_interp*b_perp1x)
        b_perpmag = (b_perp2x**2+b_perp2y**2+b_perp2z**2)**.5
        b_perp2x = b_perp2x/b_perpmag
        b_perp2y = b_perp2y/b_perpmag
        b_perp2z = b_perp2z/b_perpmag

        #print '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
        #print '%%%%%%%%%%%%%%%%%%%%%%%%%%% DEBUG %%%%%%%%%%%%%%%%%%%%%%%%%%%'
        #print 'bx,by,bz = %1.2f, %1.2f, %1.2f'%(bx_interp,by_interp,bz_interp)
        #print '1x,1y,1z = %1.2f, %1.2f, %1.2f'%(b_perp1x,b_perp1y,b_perp1z)
        #print '2x,2y,2z = %1.2f, %1.2f, %1.2f'%(b_perp2x,b_perp2y,b_perp2z)
        #print '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
        #print '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'

        #qc#print 'Rotating Velocties'
        velo={}
        for species in self.species:
            velo[species] = {}

            velo[species]['par']   = (bx_interp*self.particles[species]['vx']+by_interp*self.particles[species]['vy']+bz_interp*self.particles[species]['vz'])/bmag_interp
            velo[species]['perp1'] = self.particles[species]['vx']*b_perp1x+self.particles[species]['vy']*b_perp1y+self.particles[species]['vz']*b_perp1z
            velo[species]['perp2'] = self.particles[species]['vx']*b_perp2x+self.particles[species]['vy']*b_perp2y+self.particles[species]['vz']*b_perp2z


        #velo['e']['par']   = (bx_interp*self.particles['e']['vx']+by_interp*self.particles['e']['vy']+bz_interp*self.particles['e']['vz'])/bmag_interp
        #velo['e']['perp1'] = self.particles['e']['vx']*b_perp1x+self.particles['e']['vy']*b_perp1y+self.particles['e']['vz']*b_perp1z
        #velo['e']['perp2'] = self.particles['e']['vx']*b_perp2x+self.particles['e']['vy']*b_perp2y+self.particles['e']['vz']*b_perp2z

        #qc#print 'Generating Histograms'

        return_hist_dict = {}
        for species in self.species:
            return_hist_dict[species] = []

        for species in velo.keys():
            if kwargs.has_key('range'):
                H, xedges, yedges = np.histogram2d(velo[species]['par'],velo[species]['perp1'],**kwargs)
            else:
                H, xedges, yedges = np.histogram2d(velo[species]['par'],velo[species]['perp1'])
# H needs to be rotated and flipped
            H = np.rot90(H)
            H = np.flipud(H)
            return_hist_dict[species].append(['Parallel vs Perp 1 (+zy)',H,xedges,yedges])

            if kwargs.has_key('range'):
                H, xedges, yedges = np.histogram2d(velo[species]['par'],velo[species]['perp2'],**kwargs)
            else:
                H, xedges, yedges = np.histogram2d(velo[species]['par'],velo[species]['perp2'])
            H = np.rot90(H)
            H = np.flipud(H)
            return_hist_dict[species].append(['Parallel vs Perp 2',H,xedges,yedges])

            if kwargs.has_key('range'):
                H, xedges, yedges = np.histogram2d(velo[species]['perp1'],velo[species]['perp2'],**kwargs)
            else:
                H, xedges, yedges = np.histogram2d(velo[species]['perp1'],velo[species]['perp2'])
            H = np.rot90(H)
            H = np.flipud(H)
            return_hist_dict[species].append(['Perp 1 (+zy) vs Perp 2',H,xedges,yedges])

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
      
        for species in self.species:
            H, xedges, yedges = np.histogram2d(self.particles[species]['vx'],self.particles[species]['vy'],**kwargs)
# H needs to be rotated and flipped
            H = np.rot90(H)
            H = np.flipud(H)
            return_hist_dict[species].append(['X vs Y',H,xedges,yedges])

            H, xedges, yedges = np.histogram2d(self.particles[species]['vx'],self.particles[species]['vz'],**kwargs)
            H = np.rot90(H)
            H = np.flipud(H)
            return_hist_dict[species].append(['X vs Z',H,xedges,yedges])

            H, xedges, yedges = np.histogram2d(self.particles[species]['vy'],self.particles[species]['vz'],**kwargs)
            H = np.rot90(H)
            H = np.flipud(H)
            return_hist_dict[species].append(['Y vs Z',H,xedges,yedges])

# Mask zeros
         #Hmasked = np.ma.masked_where(H==0,H)

        return return_hist_dict

    #def get_part_in_box([location, width]):
    def _part_in_box(self,r0=[0.5,0.5],dx=[1.,1.]):
        """
        #---------------------------------------------------------------------------------------------------------------
        #   Method      : get_part_in_box
        #
        #   Discription : This method accepts a point and a width in simulation units (c/wpi) to define a box.
        #               : In that box we bin all of the particles to form the effective distrobution function
        #
        #   Args        : location [x,y] ( where you want the center of you box to be located at)
        #               : width [x,y] (the width of the box to bin particles 
        #               : dump_num (This spesifies the particular runs dump file 
        #
        #   Comments    : It would be pretty easy and potential usefull to allow this to wrap around the edges
        #               : so we can pick zero as a boundry and the code will know what todo.
        #---------------------------------------------------------------------------------------------------------------
        """
        x0 = r0[0]
        y0 = r0[1]
        if isinstance(dx,list):
            dy = dx[1]
            dx = dx[0]
        else: 
            dy = dx       # Square box is assumed

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

        xproc_lb = (int(np.floor(1.0*self.param_dict['pex']*xlb/self.param_dict['lx']))+1)%self.param_dict['nchannels']
# Colby this needs to be cooded better!!!
# if you have a diffent number of channels than pex you run into some shit
        yproc_lb = int(np.floor(1.0*self.param_dict['pey']*ylb/self.param_dict['ly']))
        #yproc_lb = int(np.floor(1.0*self.param_dict['pey']*ylb/self.param_dict['ly']))*2 #Uncomment for yishin
        #xproc_ub = (int(np.floor(1.0*self.param_dict['pex']*xub/self.param_dict['lx']))+1)%int(self.param_dict['nchannels'])
        xproc_ub = (int(np.floor(1.0*self.param_dict['pex']*xub/self.param_dict['lx']))+1)%self.param_dict['nchannels']
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
        for cosa in self.species:
            temp_dump_pruned[cosa] = []
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
                if species == 'i':
                    print '\tSelecting Ions'
                else:
                    print '\tSelecting Electrons'
# second for loop over the py
# also just doing electron for right now
                for yprocs_index in yprocs_2load:
                    self._debug = temp_dump_dat 
                    temp_dump_yproc = temp_dump_dat[species][yprocs_index]
# You only need to sort if you are on the edge processors
                    if yprocs_index == yprocs_2load[0] or yprocs_index == yprocs_2load[-1]: 
                        #qc#print '\t\tSorting yproc number '+str(yprocs_index)
                        sorted_index = temp_dump_dat[species][yprocs_index].argsort(order='y')
                        temp_dump_yproc = temp_dump_dat[species][yprocs_index][sorted_index]
# Here we need kind of a complecated if structure to get all the poible cases since
# we are scaning over muliple processors.
# If you are on your first y processor then you need to find a lower boundry
                    if yprocs_index == yprocs_2load[0]: 
                        #qc#print '\t\t\tFinding lower yboundry index '
                        lower_yboundry_index = np.searchsorted(temp_dump_yproc['y'],ylb)
                    else:
                        lower_yboundry_index = 0#np.searchsorted(temp_dump_yproc['y'],ylb)
# If you are on your last y processor then you need to find a upper boundry
                    if yprocs_index == yprocs_2load[-1]: 
                        #qc#print '\t\t\tFinding upper yboundry index '
                        upper_yboundry_index = np.searchsorted(temp_dump_yproc['y'],yub)
                    else:
                        upper_yboundry_index = -1#np.searchsorted(temp_dump_yproc['y'],yub)
                    # You only need to sort if you are on the edge processors
                    temp_dump_xproc = temp_dump_yproc[lower_yboundry_index:upper_yboundry_index]
                    if xprocs_index == xprocs_2load[0] or xprocs_index == xprocs_2load[-1]: 
                        #qc#print '\t\tNow sorting x values for remaing data'
                        sorted_index = temp_dump_xproc.argsort(order='x')
                        temp_dump_xproc = temp_dump_xproc[sorted_index] 
# If you are on your first x processor then you need to find a lower boundry
                    if xprocs_index == xprocs_2load[0]: 
                        #qc#print '\t\t\tFinding lower xboundry index '
                        lower_xboundry_index = np.searchsorted(temp_dump_xproc['x'],xlb)
                    else:
                        lower_xboundry_index = 0#np.searchsorted(temp_dump_xproc['x'],xlb)
# If you are on your last x processor then you need to find a upper boundry
                    if xprocs_index == xprocs_2load[-1]: 
                        #qc#print '\t\t\tFinding upper xboundry index '
                        upper_xboundry_index = np.searchsorted(temp_dump_xproc['x'],xub)
                    else: 
                        upper_xboundry_index = -1#np.searchsorted(temp_dump_xproc['x'],xub)
                    temp_dump_pruned[species].append(temp_dump_xproc[lower_xboundry_index:upper_xboundry_index])

        for species in self.species:
            temp_dump_pruned[species] = np.concatenate(temp_dump_pruned[species])
            print 'Total %s in box are: %i'%(species,len(temp_dump_pruned[species]))
        return temp_dump_pruned


########################################################################################################################
################################################### Note Coded Yet!!! ##################################################
########################################################################################################################
#c# #---------------------------------------------------------------------------------------------------------------
#c# #   Method      : vdist_1d
#c# #
#c# #   Discription : This method accepts a point and a width in simulation units (c/wpi) to define a box.
#c# #               : In that box we bin all of the particles to form the effective distrobution function
#c# #
#c# #   Args        : location [x,y] ( where you want the center of you box to be located at)
#c# #               : width [x,y] (the width of the box to bin particles 
#c# #               : dump_num (This spesifies the particular runs dump file 
#c# #
#c# #   Comments    : It would be pretty easy and potential usefull to allow this to wrap around the edges
#c# #               : so we can pick zero as a boundry and the code will know what todo.
#c# #---------------------------------------------------------------------------------------------------------------
#c#     #def vdist_1d([location, width,dump_index]):
#c#     def vdist_1d(self,r_simUnit=[76.8,76.8],x_width_simUnit=1.0,dump_num='000',bins=200,vflag=''):
#c# # First we handel the input
#c#         x_simUnit = r_simUnit[0]
#c#         y_simUnit = r_simUnit[0]
#c#         if isinstance(x_width_simUnit,list): # Square box is assumed
#c#             y_width_simUnit = x_width_simUnit[1]
#c#             x_width_simUnit = x_width_simUnit[0]
#c#         else:
#c#             y_width_simUnit = x_width_simUnit
#c#         if len(vflag) > 0 and vflag[0].lower() == 'v':
#c#             verbose = 1
#c#             vflag = vflag[1:]
#c#         else: 
#c#             verbose = 0
#c#         if verbose:
#c#             def verboseprint(*args):
#c#                 for arg in args:
#c#                     print arg,
#c#                 print
#c#         else:   
#c#             verboseprint = lambda *a: None      # do-nothing
#c# 
#c# # Calling get particles in box to make the vdist
#c#         verboseprint('Reading Ions and Electrons from the Dump File')
#c#         particles = self.get_part_in_box([r_simUnit,[x_width_simUnit,y_width_simUnit],dump_num,vflag])
#c# 
#c# # Here we are generating 3 histograms from the particles we read
#c#         verboseprint('Generating Histograms')
#c#         return_histx=[0.,0.]
#c#         return_histy=[0.,0.]
#c#         return_histz=[0.,0.]
#c#         return_hist_dict={}
#c#         centerx=[0.,0.]
#c#         centery=[0.,0.]
#c#         centerz=[0.,0.]
#c# 
#c# # the vx hist
#c#         return_histx[0],centerx[0] = np.histogram(particles[0]['vx'], bins=bins)
#c#         return_histx[1],centerx[1] = np.histogram(particles[1]['vx'], bins=bins)
#c#         centerx[0] = (centerx[0][1:]+centerx[0][:-1])/2.
#c#         centerx[1] = (centerx[1][1:]+centerx[1][:-1])/2.
#c#         return_hist_dict['vx'] = [return_histx,centerx]
#c# # the vy hist
#c#         return_histy[0],centery[0] = np.histogram(particles[0]['vy'], bins=bins)
#c#         return_histy[1],centery[1] = np.histogram(particles[1]['vy'], bins=bins)
#c#         centery[0] = (centery[0][1:]+centery[0][:-1])/2.
#c#         centery[1] = (centery[1][1:]+centery[1][:-1])/2.
#c#         return_hist_dict['vy'] = [return_histy,centery]
#c# # the vz hist
#c#         return_histz[0],centerz[0] = np.histogram(particles[0]['vz'], bins=bins)
#c#         return_histz[1],centerz[1] = np.histogram(particles[1]['vz'], bins=bins)
#c#         centerz[0] = (centerz[0][1:]+centerz[0][:-1])/2.
#c#         centerz[1] = (centerz[1][1:]+centerz[1][:-1])/2.
#c#         return_hist_dict['vz'] = [return_histz,centerz]
#c# 
#c#         return return_hist_dict
#c# 
#c# #---------------------------------------------------------------------------------------------------------------
#c# #   Method      : vdist_1d_par
#c# #
#c# #   Discription : This method accepts a point and a width in simulation units (c/wpi) to define a box.
#c# #               : In that box we bin all of the particles to form the effective distrobution function
#c# #
#c# #   Args        : location [x,y] ( where you want the center of you box to be located at)
#c# #               : width [x,y] (the width of the box to bin particles 
#c# #               : dump_num (This spesifies the particular runs dump file 
#c# #
#c# #   Comments    : It would be pretty easy and potential usefull to allow this to wrap around the edges
#c# #               : so we can pick zero as a boundry and the code will know what todo.
#c# #---------------------------------------------------------------------------------------------------------------
#c#     #def vdist_1d_par([location, width,dump_index]):
#c#     def vdist_1d_par(self,r_simUnit=[76.8,76.8],x_width_simUnit=1.0,dump_num='000',bins=51,vflag=''):
#c# # First we handel the input
#c#         x_simUnit = r_simUnit[0]
#c#         y_simUnit = r_simUnit[0]
#c#         if isinstance(x_width_simUnit,list): # Square box is assumed
#c#             y_width_simUnit = x_width_simUnit[1]
#c#             x_width_simUnit = x_width_simUnit[0]
#c#         else:
#c#             y_width_simUnit = x_width_simUnit
#c#         if len(vflag) > 0 and vflag[0].lower() == 'v':
#c#             verbose = 1
#c#             vflag = vflag[1:]
#c#         else: 
#c#             verbose = 0
#c#         if verbose:
#c#             def verboseprint(*args):
#c#                 for arg in args:
#c#                     print arg,
#c#                 print
#c#         else:   
#c#             verboseprint = lambda *a: None      # do-nothing
#c# 
#c# # Calling get particles in box to make the vdist
#c#         verboseprint('Reading Ions and Electrons from the Dump File')
#c#         particles = self.get_part_in_box([r_simUnit,[x_width_simUnit,y_width_simUnit],dump_num,vflag])
#c# # Reading in fields to calculate vpar
#c#         verboseprint('Reading in the Fields form the Dump File')
#c#         dump_field_dict = self.read_fields_from_dump([dump_num])
#c# # Now there are two ways to do this, a faster and a slower
#c# # Faster: just take the average B field and and use that
#c# #         for every particle
#c# # Slower: Use an interpolated Bfield for each point
#c# #
#c# # Right now im only coding the faster one
#c#         verboseprint('Calculating B field')
#c#         
#c#         bx_interp = interp_field(dump_field_dict['bx'],self.param_dict['lx'],self.param_dict['ly'],r_simUnit[0],r_simUnit[1])
#c#         by_interp = interp_field(dump_field_dict['by'],self.param_dict['lx'],self.param_dict['ly'],r_simUnit[0],r_simUnit[1])
#c#         bz_interp = interp_field(dump_field_dict['bz'],self.param_dict['lx'],self.param_dict['ly'],r_simUnit[0],r_simUnit[1])
#c#         bmag_interp = (bx_interp**2+by_interp**2+bz_interp**2)**(.5)
#c# 
#c#         b_perp1y = -1.*bz_interp/(by_interp**2+bz_interp**2)**(.5)
#c#         b_perp1z = by_interp/(by_interp**2+bz_interp**2)**(.5)
#c# 
#c#         b_perp2x = (by_interp*b_perp1z - bz_interp*b_perp1y)
#c#         b_perp2y = (-1.*bx_interp*b_perp1z)
#c#         b_perp2z = (bx_interp*b_perp1y)
#c#         b_perpmag = ((by_interp*b_perp1z - bz_interp*b_perp1y)**2+(bx_interp*b_perp1z)**2+(bx_interp*b_perp1y)**2)**(.5)
#c#         b_perp2x = b_perp2x/b_perpmag
#c#         b_perp2y = b_perp2y/b_perpmag
#c#         b_perp2z = b_perp2z/b_perpmag
#c#         
#c#         verboseprint('Rotating Velocties')
#c# # Also just doing the electons to start off with 
#c#         vpar = (bx_interp*particles[1]['vx']+by_interp*particles[1]['vy']+bz_interp*particles[1]['vz'])/bmag_interp
#c#         vperp1 = particles[1]['vy']*b_perp1y+particles[1]['vz']*b_perp1z
#c#         vperp2 = particles[1]['vx']*b_perp2x+particles[1]['vy']*b_perp2y+particles[1]['vz']*b_perp2z
#c# 
#c#         verboseprint('Generating Histograms')
#c# 
#c#         return_hist_dict = {}
#c#         return_hist = [0.,0.]
#c#         extent = [[0.,0.,0.,0.],[0.,0.,0.,0.]]
#c#         return_hist[0],xedge,yedge = np.histogram2d(particles[0]['vx'], particles[0]['vy'], bins=bins)
#c#         extent[0] = [xedge[0],xedge[-1],yedge[0],yedge[-1]]
#c#         return_hist[1],extent[1] = np.histogram(vpar,bins=bins)
#c#         #extent[1] = [xedge[0],xedge[-1],yedge[0],yedge[-1]]
#c#         return_hist_dict = [return_hist[1],(extent[1][1:]+extent[1][:-1])/2.]
#c# #@# # the vx vs vy hist
#c# #@#         return_hist[0],xedge,yedge = np.histogram2d(particles[0]['vx'], particles[0]['vy'], bins=bins)
#c# #@#         extent[0] = [xedge[0],xedge[-1],yedge[0],yedge[-1]]
#c# #@#         return_hist[1],xedge,yedge = np.histogram2d(particles[1]['vx'], particles[1]['vy'], bins=bins)
#c# #@#         extent[1] = [xedge[0],xedge[-1],yedge[0],yedge[-1]]
#c# #@#         return_hist_dict['vxy'] = [return_hist,extent]
#c# #@# # the vx vs vz hist
#c# #@#         return_hist[0],xedge,yedge = np.histogram2d(particles[0]['vx'], particles[0]['vz'], bins=bins)
#c# #@#         extent[0] = [xedge[0],xedge[-1],yedge[0],yedge[-1]]
#c# #@#         return_hist[1],xedge,yedge = np.histogram2d(particles[1]['vx'], particles[1]['vz'], bins=bins)
#c# #@#         extent[1] = [xedge[0],xedge[-1],yedge[0],yedge[-1]]
#c# #@#         return_hist_dict['vxz'] = [return_hist,extent]
#c# #@# # the vy vs vz hist
#c# #@#         return_hist[0],xedge,yedge = np.histogram2d(particles[0]['vy'], particles[0]['vz'], bins=bins)
#c# #@#         extent[0] = [xedge[0],xedge[-1],yedge[0],yedge[-1]]
#c# #@#         return_hist[1],xedge,yedge = np.histogram2d(particles[1]['vy'], particles[1]['vz'], bins=bins)
#c# #@#         extent[1] = [xedge[0],xedge[-1],yedge[0],yedge[-1]]
#c# #@#         return_hist_dict['vyz'] = [return_hist,extent]
#c# 
#c#         return return_hist_dict
#c# 
#c# #---------------------------------------------------------------------------------------------------------------
#c# #   Method      : vdist_2d
#c# #
#c# #   Discription : This method accepts a point and a width in simulation units (c/wpi) to define a box.
#c# #               : In that box we bin all of the particles to form the effective distrobution function
#c# #
#c# #   Args        : location [x,y] ( where you want the center of you box to be located at)
#c# #               : width [x,y] (the width of the box to bin particles 
#c# #               : dump_num (This spesifies the particular runs dump file 
#c# #
#c# #   Comments    : It would be pretty easy and potential usefull to allow this to wrap around the edges
#c# #               : so we can pick zero as a boundry and the code will know what todo.
#c# #---------------------------------------------------------------------------------------------------------------
#c#     #def vdist_2d([location, width,dump_index]):
#c#     def vdist_2d(self,r_simUnit=[76.8,76.8],x_width_simUnit=1.0,dump_num='000',bins=51,vflag=''):
#c# # First we handel the input
#c#         x_simUnit = r_simUnit[0]
#c#         y_simUnit = r_simUnit[0]
#c#         if isinstance(x_width_simUnit,list): # Square box is assumed
#c#             y_width_simUnit = x_width_simUnit[1]
#c#             x_width_simUnit = x_width_simUnit[0]
#c#         else:
#c#             y_width_simUnit = x_width_simUnit
#c#         if len(vflag) > 0 and vflag[0].lower() == 'v':
#c#             verbose = 1
#c#             vflag = vflag[1:]
#c#         else: 
#c#             verbose = 0
#c#         if verbose:
#c#             def verboseprint(*args):
#c#                 for arg in args:
#c#                     print arg,
#c#                 print
#c#         else:   
#c#             verboseprint = lambda *a: None      # do-nothing
#c# 
#c# # Calling get particles in box to make the vdist
#c#         verboseprint('Reading Ions and Electrons from the Dump File')
#c#         particles = self.get_part_in_box([r_simUnit,[x_width_simUnit,y_width_simUnit],dump_num,vflag])
#c# 
#c# # Here we are generating 3 histograms from the particles we read
#c#         verboseprint('Generating Histograms')
#c#         return_hist_dict = {}
#c#         return_histxy = [0.,0.]
#c#         extentxy = [[0.,0.,0.,0.],[0.,0.,0.,0.]]
#c#         return_histxz = [0.,0.]
#c#         extentxz = [[0.,0.,0.,0.],[0.,0.,0.,0.]]
#c#         return_histyz = [0.,0.]
#c#         extentyz = [[0.,0.,0.,0.],[0.,0.,0.,0.]]
#c# 
#c# # the vx vs vy hist
#c#         return_histxy[0],xedge,yedge = np.histogram2d(particles[0]['vx'], particles[0]['vy'], bins=bins)
#c#         extentxy[0] = [xedge[0],xedge[-1],yedge[0],yedge[-1]]
#c#         return_histxy[1],xedge,yedge = np.histogram2d(particles[1]['vx'], particles[1]['vy'], bins=bins)
#c#         extentxy[1] = [xedge[0],xedge[-1],yedge[0],yedge[-1]]
#c#         return_hist_dict['vxy'] = [return_histxy,extentxy]
#c# # the vx vs vz hist
#c#         return_histxz[0],xedge,yedge = np.histogram2d(particles[0]['vx'], particles[0]['vz'], bins=bins)
#c#         extentxz[0] = [xedge[0],xedge[-1],yedge[0],yedge[-1]]
#c#         return_histxz[1],xedge,yedge = np.histogram2d(particles[1]['vx'], particles[1]['vz'], bins=bins)
#c#         extentxz[1] = [xedge[0],xedge[-1],yedge[0],yedge[-1]]
#c#         return_hist_dict['vxz'] = [return_histxz,extentxz]
#c# # the vy vs vz hist
#c#         return_histyz[0],xedge,yedge = np.histogram2d(particles[0]['vy'], particles[0]['vz'], bins=bins)
#c#         extentyz[0] = [xedge[0],xedge[-1],yedge[0],yedge[-1]]
#c#         return_histyz[1],xedge,yedge = np.histogram2d(particles[1]['vy'], particles[1]['vz'], bins=bins)
#c#         extentyz[1] = [xedge[0],xedge[-1],yedge[0],yedge[-1]]
#c#         return_hist_dict['vyz'] = [return_histyz,extentyz]
#c# 
#c#         return return_hist_dict
#c# 
#c# #---------------------------------------------------------------------------------------------------------------
#c# #   Method      : vdist_2d_par
#c# #
#c# #   Discription : This method accepts a point and a width in simulation units (c/wpi) to define a box.
#c# #               : In that box we bin all of the particles to form the effective distrobution function
#c# #
#c# #   Args        : location [x,y] ( where you want the center of you box to be located at)
#c# #               : width [x,y] (the width of the box to bin particles 
#c# #               : dump_num (This spesifies the particular runs dump file 
#c# #
#c# #   Returns     : return_hist_dict['axis'][i/e][hist_array]
#c# #               :   'axis'  slecets the data for the 2d projection of the distabution function 
#c# #                   'axis' = 'parperp1', 'parperp1', 'perp1perp2'
#c# #                   i/e is an int that selects either the electron or ion data. i = 0, e = 1
#c# #                   The rest is just an object that is returned directly from the histogram2d method
#c# #
#c# #
#c# #   Comments    : It would be pretty easy and potential usefull to allow this to wrap around the edges
#c# #               : so we can pick zero as a boundry and the code will know what todo.
#c# #---------------------------------------------------------------------------------------------------------------
#c#     #def vdist_2d_par([location, width,dump_index]):
#c#     def vdist_2d_par(self,r_simUnit=[76.8,76.8],x_width_simUnit=1.0,dump_num='000',bins=51,vflag=''):
#c# # First we handel the input
#c#         x_simUnit = r_simUnit[0]
#c#         y_simUnit = r_simUnit[0]
#c#         if isinstance(x_width_simUnit,list): # Square box is assumed
#c#             y_width_simUnit = x_width_simUnit[1]
#c#             x_width_simUnit = x_width_simUnit[0]
#c#         else:
#c#             y_width_simUnit = x_width_simUnit
#c#         if len(vflag) > 0 and vflag[0].lower() == 'v':
#c#             verbose = 1
#c#             vflag = vflag[1:]
#c#         else: 
#c#             verbose = 0
#c#         if verbose:
#c#             def verboseprint(*args):
#c#                 for arg in args:
#c#                     print arg,
#c#                 print
#c#         else:   
#c#             verboseprint = lambda *a: None      # do-nothing
#c# 
#c# # Calling get particles in box to make the vdist
#c#         verboseprint('Reading Ions and Electrons from the Dump File')
#c#         particles = self.get_part_in_box([r_simUnit,[x_width_simUnit,y_width_simUnit],dump_num,vflag])
#c# # Reading in fields to calculate vpar
#c#         verboseprint('Reading in the Fields form the Dump File')
#c#         dump_field_dict = self.read_fields_from_dump([dump_num])
#c# # Now there are two ways to do this, a faster and a slower
#c# # Faster: just take the average B field and and use that
#c# #         for every particle
#c# # Slower: Use an interpolated Bfield for each point
#c# #
#c# # Right now im only coding the faster one
#c#         verboseprint('Calculating B field')
#c#         
#c#         bx_interp = interp_field(dump_field_dict['bx'],self.param_dict['lx'],self.param_dict['ly'],r_simUnit[0],r_simUnit[1])
#c#         by_interp = interp_field(dump_field_dict['by'],self.param_dict['lx'],self.param_dict['ly'],r_simUnit[0],r_simUnit[1])
#c#         bz_interp = interp_field(dump_field_dict['bz'],self.param_dict['lx'],self.param_dict['ly'],r_simUnit[0],r_simUnit[1])
#c#         bmag_interp = (bx_interp**2+by_interp**2+bz_interp**2)**.5
#c# 
#c#         if by_interp > 0.:
#c#             b_perp1x = 0.
#c#             b_perp1y = -1.*bz_interp/(bz_interp**2 + by_interp**2)**(.5)
#c#             b_perp1z = by_interp/(bx_interp**2 + by_interp**2)**(.5)
#c#         else:
#c#             b_perp1x = 0.
#c#             b_perp1y = bz_interp/(bz_interp**2 + by_interp**2)**(.5)
#c#             b_perp1z = -1.*by_interp/(bx_interp**2 + by_interp**2)**(.5)
#c# 
#c#         b_perp2x = (by_interp*b_perp1z - bz_interp*b_perp1y)
#c#         b_perp2y = (bz_interp*b_perp1x - bx_interp*b_perp1z)
#c#         b_perp2z = (bx_interp*b_perp1y - by_interp*b_perp1x)
#c#         b_perpmag = (b_perp2x**2+b_perp2y**2+b_perp2z**2)**.5
#c#         b_perp2x = b_perp2x/b_perpmag
#c#         b_perp2y = b_perp2y/b_perpmag
#c#         b_perp2z = b_perp2z/b_perpmag
#c# 
#c#         verboseprint('Rotating Velocties')
#c# # Also just doing the electons to start off with 
#c#         vpar_i = (bx_interp*particles[0]['vx']+by_interp*particles[0]['vy']+bz_interp*particles[0]['vz'])/bmag_interp
#c#         vperp1_i = particles[0]['vx']*b_perp1x+particles[0]['vy']*b_perp1y+particles[0]['vz']*b_perp1z
#c#         vperp2_i = particles[0]['vx']*b_perp2x+particles[0]['vy']*b_perp2y+particles[0]['vz']*b_perp2z
#c# 
#c#         vpar_e = (bx_interp*particles[1]['vx']+by_interp*particles[1]['vy']+bz_interp*particles[1]['vz'])/bmag_interp
#c#         vperp1_e = particles[1]['vx']*b_perp1x+particles[1]['vy']*b_perp1y+particles[1]['vz']*b_perp1z
#c#         vperp2_e = particles[1]['vx']*b_perp2x+particles[1]['vy']*b_perp2y+particles[1]['vz']*b_perp2z
#c# 
#c#         verboseprint('Generating Histograms')
#c# 
#c#         return_hist_dict = {}
#c#         
#c#         return_hist_dict['parperp1'] = [np.histogram2d(vperp1_i,vpar_i, bins=bins),np.histogram2d(vperp1_e,vpar_e, bins=bins)]
#c#         return_hist_dict['parperp2'] = [np.histogram2d(vperp2_i,vpar_i, bins=bins),np.histogram2d(vperp2_e,vpar_e, bins=bins)]
#c#         return_hist_dict['perp1perp2'] = [np.histogram2d(vperp2_i,vperp1_i, bins=bins),np.histogram2d(vperp2_e,vperp1_e, bins=bins)]
#c# 
#c#         return return_hist_dict

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


    def interp_field(self,field,lx='not_set',ly='not_set',r0='not_set'):
        """
#---------------------------------------------------------------------------------------------------------------
#   Method      : interp_field
#
#   Discription : This method takes a field and a floating point, and returns the linear fit value 
#               : between the grid points
#
#   Args        : field  The field you are interpolating
#               : r0[0] The xpoint to interpolate at
#               : r0[1] The ypoint to interpolate at
#
#   Comments    : I think this is working ok? It would be smart to make this an object method that just reads
#               : the internally saved field. so CODE IN THE FUTURE
#---------------------------------------------------------------------------------------------------------------
        """
        if lx=='not_set':lx=self.param_dict['lx']
        if ly=='not_set':ly=self.param_dict['ly']
        if r0=='not_set':r0=self._r0

        nx = len(field[0,:])
        ny = len(field[:,0])
        ip = int(np.floor(1.0*r0[0]/lx*nx))
        jp = int(np.floor(1.0*r0[1]/ly*ny))

        if ip + 1 > nx-1: ipp1 = 0
        else: ipp1 = ip+1

        if jp + 1 > ny-1: jpp1 = 0
        else: jpp1 = jp+1
        print ip,jp,nx,ny
        weight_x = 1.0*r0[0]/lx*nx - np.floor(1.0*r0[0]/lx*nx)
        weight_y = 1.0*r0[1]/ly*ny - np.floor(1.0*r0[1]/ly*ny)
        return np.array([(1.-weight_x)*(1.-weight_y)*field[jp,ip] + (weight_x)*(1.-weight_y)*field[jpp1,ip] \
               + (1.-weight_x)*(weight_y)*field[jp,ipp1]+ (weight_x)*(weight_y)*field[jpp1,ipp1]])
    


#Orphaned method that I use in load param but not quite sure how to fit it in
def convert(val):
    constructors = [int, float, str]
    for c in constructors:
        try:
            return c(val)
        except ValueError:
            pass
            





        

