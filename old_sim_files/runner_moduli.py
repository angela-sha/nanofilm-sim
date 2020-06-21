#imports
from analyzer import Analyzer
from SimulatorRL import SimulatorRL
from simulator import Simulator
from PressureSensorRL import PressureSensorRL
from PressureSensor import PressureSensor
from Grapher import Grapher
from BinObject import BinObject

import math
import sys
import os
sys.path = sys.path + ['/home/cmliepold/Dans_GPU/md_engine-master7/build/python/build/lib.linux-x86_64-2.7']
from Sim import *
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import Voronoi
from scipy import interpolate
import matplotlib.animation as manimation

#needs to be fixed
temperature = .9 	#reduced; need to figure out the non-reduced another day
density = .93

#sigma, epsilon, alpha, r_a, a_g, sig_g, r_g
param_list = [(1, 1, .11268, 1.1051, .96, 1.00547, .05153)]	#parameters after being reduced, used for RL
unreduced_param_list = (340, 5.298, 758789.8)				#epsilon, sigma, mass
#unrounded 110.2586974
box_params = (110.26, 173.27, 110.26)						#size of the box

command_list = [(100, 4000), (0, 5000, 531250)]
#command_list = [(100, 4000), (0, 100, 531250)] 				#parameters for simulations

#graphing options (Grapher object)
#g = Grapher("component")
#g = Grapher("binhua_method")
#g = Grapher("component")
g = Grapher("just_plot_the_goddamn_graph")

#run simulation  
sim = Simulator(box_params, density, temperature, 1, Vector(-.0166,0,0), param_list[0], command_list)
fn = "attempt"
sim.update_fn(fn)
sim.run()
file = fn + ".xyz"

#add pressure sensors
ps = PressureSensorRL((50.1182, 60.1418, 76.6113, 96.6586), 2, param_list[0][2:], unreduced_param_list, g)
pressure_sensor_list = [ps]

ana = Analyzer(unreduced_param_list, pressure_sensor_list, fn)

f = open(file, "r")
number_of_particles = f.readline().strip()
print(number_of_particles)

ana.perform_analysis(f, number_of_particles)

f.close()

#g.full_graph()
g.graph("Graph Title", "graph")

