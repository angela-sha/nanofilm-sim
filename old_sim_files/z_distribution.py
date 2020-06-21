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

fn = "attempt3D_1"

temperature = .9 			#reduced; need to figure out the non-reduced another day
density = .93

g = Grapher("convolution")
box_params = (110.26, 173.27, 110.26)
unreduced_param_list = (340, 5.298, 758789.8)
param_list = [(1, 1, .11268, 1.1051, .96, 1.00547, .05153)]

file = fn + ".xyz"

#add the pressure sensors
ps = PressureSensorRL((50.1182, 60.1418, 76.6113, 96.6586), 2, param_list[0][2:], unreduced_param_list, g)
pressure_sensor_list = [ps]

ana = Analyzer( unreduced_param_list, pressure_sensor_list, fn )

f = open(file, "r")
number_of_particles = f.readline().strip()
print(number_of_particles)

ana.perform_analysis(f, number_of_particles)
f.close()