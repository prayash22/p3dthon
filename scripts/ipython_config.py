# This is the ipython config file, this will load up good python stuff on start up
# Put this file in '~/.ipython/profile_default/startup/'
# and lastly change the path p3dthon_path to where ever you keep your own fork 
import sys
p3dthon_path = '/glade/u/home/colbyh/pythonprogs/2014.04.p3d_etal/p3dthon/'
sys.path.append(p3dthon_path+'objects/')
sys.path.append(p3dthon_path+'scripts/')
from p3d_runs import p3d_run
import sub
from  sub import *


