# This is the ipython config file, this will load up good python stuff on start up
# Put this file in '~/.ipython/profile_default/startup/'
# and lastly change the path p3dthon_path to where ever you keep your own fork 
import sys
import os
import time
p3dthon_path = '/glade/u/home/colbyh/p3dthon/'
sys.path.append(p3dthon_path+'objects/')
sys.path.append(p3dthon_path+'objects/testparticles/')
sys.path.append(p3dthon_path+'scripts/')
from testparticle import TPRun
from p3d_runs import p3d_run
from  sub import *


