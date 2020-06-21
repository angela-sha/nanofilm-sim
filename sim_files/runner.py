from simulator import Simulator
from simulatorRL import SimulatorRL 

import math
import sys
import os
#sys.path = sys.path + ['/home/alexfeistritzer/dash_attempt/md_engine/build/python/build/lib.linux-x86_64-2.7']
sys.path = sys.path + ['/home/cmliepold/Dans_GPU/md_engine-master7/build/python/build/lib.linux-x86_64-2.7']
from Sim import *
import numpy as np

#reduced temperature and density
temperature = .9
density = .93

#sigma, epsilon, alpha, r_a, a_g, sig_g, r_g
params = (1, 1, .11268, 1.1051, .96, 1.00547, .05153)	#reduced parameters
unred_params = (340, 58.6, 758789.8)					#epsilon, sigma, mass

#unreduced_param_list = (340, 5.298, 758789.8)			#epsilon, sigma, mass
#unrounded 110.2586974

box_params = (110.26, 173.27, 110.26)					#size of the box
command_list = [(100, 4000), (100, 100125)]				#100 ts, 4000 ts equilibriate; 400 ts, 53125 ts
#command_list = [(100, 4000), (100, 531250)] 			#parameters for simulations

sim = Simulator(box_params, density, temperature, 1, Vector(-.166,0,0), params, command_list)
fn = "attempt"
sim.update_fn(fn)
sim.run()
file = fn + ".xyz"