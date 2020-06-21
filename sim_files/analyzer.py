'''
This class runs on an xyz file (taken as the second argument).
'''
from parse_xyz import Trajectory
from process_coord import Coordinate
#from lattice import Lattice

import numpy as np
import sys

params = (1, 1, .11268, 1.1051, .96, 1.00547, .05153)		#sigma, epsilon, alpha, r_a, a_g, sig_g, r_g
unred_params = (340, 58.6, 758789.8)						#epsilon, sigma, mass

#parse file, store trajectory data in numpy array
in_file = sys.argv[1]
pos = Trajectory(in_file, unred_params)
num_ts = pos.parse_file()
num_particles = pos.get_num_particles()
area_data = pos.get_area_data()				#area is unreduced, in nm^2
pos_data = pos.get_pos_data()				#positions are reduced

print("Timesteps: " + str(num_ts))
print("Particles: " + str(num_particles))

#start analyzing coordinate date (z-positions, position graphs)
coord = Coordinate(num_particles, num_ts, area_data, pos_data, unred_params)
#coord.cluster_kde(1041, True, 'kde_graph_1041.png')
'''
coord.graph_particle_path(300, 40, 780, 'particle_300_path.png')
coord.graph_particle_path(400, 40, 780, 'particle_400_path.png')
coord.graph_particle_path(600, 40, 780, 'particle_600_path.png')
coord.graph_particle_path(700, 40, 780, 'particle_700_path.png')
'''
#test functions here: get_threshold_value, check_layer, graph_pos_dist

#graph coordinate data here
#coord.graph_pos_dist(441, 'plot_441.png')
coord.graph_2D(300, 200, 250, 600, False, 'graph_test_300.png')
coord.graph_2D(400, 200, 250, 600, False, 'graph_test_400.png')
coord.graph_2D(500, 200, 250, 600, False, 'graph_test_500.png')
coord.graph_2D(600, 200, 250, 600, False, 'graph_test_600.png')
coord.graph_2D(700, 200, 250, 600, False, 'graph_test_700.png')
coord.graph_2D(800, 200, 250, 600, False, 'graph_test_800.png')
coord.graph_2D(900, 200, 250, 600, False, 'graph_test_900.png')
coord.graph_2D(950, 200, 250, 600, False, 'graph_test_950.png')
coord.graph_2D(1000, 200, 250, 600, False, 'graph_test_1000.png')
coord.graph_2D(1040, 200, 250, 600, False, 'graph_test_1040.png')


