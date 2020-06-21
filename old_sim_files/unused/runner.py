'''The runner script is meant to be the only thing that needs edited from simulation to simulation. In this way,
we reduce the amount of changes that need to be made and also make it easier to run multiple simulations simultaneously.'''

# -----> imports

# Other objects that I've creaetd
from analyzer import Analyzer
from SimulatorRL import SimulatorRL
from simulator import Simulator
from PressureSensorRL import PressureSensorRL
from PressureSensor import PressureSensor
from Grapher import Grapher
from BinObject import BinObject

import sys

# location of the DASH simulator on the supercomputing cluster
sys.path = sys.path + ['/home/cmliepold/Dans_GPU/md_engine-master7/build/python/build/lib.linux-x86_64-2.7']
from Sim import *

#I don't even know if these need to be here, but these, in general, are part of all objects
import math
import os
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import Voronoi
from scipy import interpolate
import matplotlib.animation as manimation


#This definitely still needs fixed
temperature = .9 			#reduced; need to figure out the non-reduced another day
density = .93

# -----> Defining the simulation parameters
# so in order, these represent (sigma, epsilon, alpha, r_a, a_g, sig_g, r_g)
param_list = [(1, 1, .11268, 1.1051, .96, 1.00547, .05153)]	#parameters after being reduced, used for RL


unreduced_param_list = (340, 5.298, 758789.8)				#epsilon, sigma, mass

#110.2586974
box_params = (110.26, 173.27, 110.26)						#size of the box, (x,y,z)

command_list = [(100, 4000), (0, 100, 531250)] #1000000				#how to run the simulations themselves

g = Grapher("convolution")

#added in order to repeat a simulation multiple times
for i in range(1):
	for p in param_list:

		sim = SimulatorRL( box_params, density, temperature, 1, Vector(-.0166,0,0), p, command_list )
		fn = "attempt3D_1"		#filename
		sim.update_fn(fn)
		sim.run()

		
		file = fn + ".xyz"

		#add the pressure sensors
		ps = PressureSensorRL((50.1182, 60.1418, 76.6113, 96.6586), 2, p[2:], unreduced_param_list, g)
		pressure_sensor_list = [ps]

		ana = Analyzer( unreduced_param_list, pressure_sensor_list, fn )

		f = open(file, "r")
		number_of_particles = f.readline().strip()
		print(number_of_particles)

		ana.perform_analysis(f, number_of_particles)
		f.close()
		

#g.full_graph()




