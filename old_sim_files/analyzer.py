'''
An instance of this class is run after a simulation has been completed. It is in charge of giving the pressure defined by any and all pressure sensors that are instantiated, and essentially acts as a universal holder for pressure sensors.
'''
import matplotlib
matplotlib.use('Agg')
from scipy.spatial import Voronoi
import math
#import matplotlib.animation as manimation
import matplotlib.pyplot as plt
import numpy as np

class Analyzer:

	'''Constructor for the Analyzer object.
	Parameters:
		unreduced_params - the reducing factors, i.e. epsilon, sigma, mass before they are converted to one. Sent as a tuple of (epsilon, sigma, mass)
		pressure_sensor_list - a list of all of the pressure sensor instances to be used in this analysis
		fn - filename, for the graph generated
	Returns:
		null
	Note:
		number_of_particles is updated later on, when perform analysis is called.
	'''
	def __init__(self, unreduced_params, pressure_sensor_list, fn):
		self.number_of_particles = 0
		self.unred_eps = unreduced_params[0]
		self.unred_sig = unreduced_params[1]
		self.unred_mass = unreduced_params[2]
		self.pressure_sensor_list = pressure_sensor_list
		self.fn = fn

	'''Helper function for parse_file(). Gives the bounds of the area of the simulation, given the raw file string from the .xyz file that lists the bounds of the simulation box
	Parameters: 
		str - the raw string from the .xyz file for simulation box bounds
	Returns: 
		a float of the area of the simulation
	'''
	def find_bounds(self, str):

		#get the numerical data
		bounds_raw    = str.strip().split()
		bounds_lo_raw = bounds_raw[2] + bounds_raw[3] + bounds_raw[4]		
		bounds_hi_raw = bounds_raw[6] + bounds_raw[7] + bounds_raw[8]

		#remove non-numeric text
		bounds_lo_str = bounds_lo_raw[1:-1].split(",")
		bounds_hi_str = bounds_hi_raw[1:-1].split(",")

		#parse to float
		bounds_lo = [float(i) for i in bounds_lo_str]
		bounds_hi = [float(j) for j in bounds_hi_str]

		#calculate and return area
		return abs(bounds_hi[0] - bounds_lo[0]) * abs(bounds_hi[1] - bounds_lo[1]) * ((self.unred_sig) ** 2) / float(self.number_of_particles) 



	'''Parses the data from the .xyz file out of the file into the format [[simulation_area, (xyz), ...], [...], ...].
	Parameters: 
		file - an open .xyz file, in read mode
	Returns: 
		list in the form described above, where simulation_area is a float and (xyz) is a tuple of the x, y, z position of a particle
	'''
	def parse_file(self, file):
		simulation_holder = []							#stores all of the timesteps
		time_step_holder = []							#stores an individual timestep

		#get and save area of simulation box
		bounds = self.find_bounds(file.readline())
		time_step_holder.append(bounds)
		print(bounds)

		#determines a new timestep
		reset = False		

		for line in file:
			if line.strip() == str(self.number_of_particles):		#check to see if the next timestep has been reached
				reset = True
				continue

			if reset:					#resets time_step_holder and appends its data to simulation_holder
				#simulation_holder.append(time_step_holder)
				yield time_step_holder
				time_step_holder = []
				reset = False

				#get and save area of simulation box
				bounds = self.find_bounds(line)
				time_step_holder.append(bounds)

				continue

			#parses xyz data into a tuple and appends it to time_step_holder
			line = line.strip().split()
			line = line[1:]
			line = [float(i) for i in line]
			pos = tuple(line)
			time_step_holder.append(pos)

		#return simulation_holder
		yield time_step_holder




	'''Used to find all the particles in a given region.
	Parameters: 
		area_to_search - list of parsed particle positions from which the points are pulled from
		bounds - what area is searched, given as a tuple of (x1, x2, y1, y2)
	Returns: 
		list with all of the particles within the area, each as a tuple of (x, y, z)
	'''
	def find_particles_in(self, area_to_search, bounds):
		particles_in_region = []				#where to store the particles to return

		x1 = bounds[0]
		x2 = bounds[1]
		y1 = bounds[2]
		y2 = bounds[3]
		for particle in area_to_search:
			x = particle[0]
			y = particle[1]

			#check if in bounds
			if (x >= x1 and x <= x2 and y >= y1 and y <= y2):
				particles_in_region.append(particle)


		return particles_in_region




	'''Finds the next-nearest neighbors, given a dictionary of nearest neighbors
	Parameters: 
		dict - dictionary of nearest neighbors
	Returns: 
		dictionary of nearest and next-nearest neighbors
	'''
	def find_next_nearest(self, dict):
		new_dict = {}

		for point, neighbor_list in dict.iteritems():			#finds all neighbors
			search = []
			for n in neighbor_list:
				search.append(n)

			new_neighbor_list = []
			for n in search:									#finds the next-nearest neighbors
				nnn_list = dict[n]
				for nnn in nnn_list:
					new_neighbor_list.append(nnn)

			new_neighbor_list = list(set(new_neighbor_list))

			new_dict[point] = new_neighbor_list

		return new_dict




	'''Given a list of particles, it will calculate all of the nearest neighbors, and optionally, next-nearest neighbors, returning a dictionary where the particle position tuple is the key and a list of particle tuples that consist of the nearest (and next-nearest) neighbors
	Parameters: 
		positions: list of particle position in tuple form
		next_nearest: boolean if next-nearest neighbors are to be returned as well
	Returns: 
		dictionary of { particle: [neighbors] }
	'''
	def find_nearest_neighbors(self, positions, next_nearest):

		#xy_position = [(positions[i][0], positions[i][1]) for i in range(len(positions)) ]	#xy pos only
		xyz_position = [(positions[i][0], positions[i][1], positions[i][2]) for i in range(len(positions))]

		nn_dict = {}

		#vor = Voronoi(xy_position)			#generate a Voronoi diagram of the particles
		vor = Voronoi(xyz_position)

		#finds half of the region (only one-way)
		for points in vor.ridge_points:		
			key = positions[points[0]]		

			#checks to see if the particle is already in the dictionary, and adds it
			if key in nn_dict:
				old = nn_dict[key]
				old.append(positions[points[1]])
				val = old
				nn_dict[key] = val
			else:
				nn_dict[key] = [ positions[points[1]] ]

		#repeats after inerting the Voronoi list
		new_vor = [ [i[1], i[0]] for i in vor.ridge_points ]
		for points in new_vor:
			key = positions[points[0]]

			if key in nn_dict:
				old = nn_dict[key]
				old.append(positions[points[1]])
				val = old
				nn_dict[key] = val
			else:
				nn_dict[key] = [ positions[points[1]] ]

		#include next-nearest, if applicable
		if next_nearest:
			full_nn_dict = self.find_next_nearest(nn_dict)
			return full_nn_dict

		return nn_dict



	'''Runs the actual analysis for a simulation.
	Parameters:
		file_name - the (open, in read mode) file where the .xyz data is stored
		number_of_molecules - the number of molecules in th simulation, parsed from the file, and used for post-analysis calculations
	Returns:
		null
	'''
	def perform_analysis(self, file_name, number_of_molecules):
		self.number_of_particles = number_of_molecules

		print("Running an analysis")

		gen_fun = self.parse_file(file_name)
		counter = 0

		#FFMpegWriter = manimation.writers['ffmpeg']
		#metadata = dict(title='Movie Test', artist='Matplotlib', comment='Movie support!')
		#writer = FFMpegWriter(fps=15, metadata=metadata)

		#fig = plt.figure()
		#im, = plt.imshow((x, y, e), animated=True)

		#with writer.saving(fig, "entropy_test1.mp4", 100):
		for ts in gen_fun:

			#Removed for now, but kinda essential
			
			#for ts in timestep_data:
			area = ts[0]

			for sensor in self.pressure_sensor_list:
				particle_list = self.find_particles_in(ts[1:], sensor.frame_bounds)
				particle_dict = self.find_nearest_neighbors(particle_list, False)

				#sensor.calc_pressure2(particle_dict, area)
				sensor.calc_pot(particle_dict, area)

			#think about adding an "out of" counter
			print("Finished timestep " + str(counter))
			counter+=1
			
			'''
			#here's the entropy stuff
			if counter % 100 == 0:

				plt.figure(1)
				z_pos = [ts[1:][i][2] for i in range(len(ts[1:]))]
				print(np.mean(z_pos))
				print(np.std(z_pos))
				plt.plot(z_pos, marker='o', linestyle='None') #marker='.'
				#plt.imshow(np.array((entropy_x_norm, entropy_y_norm), dtype='uint8'))
				plt.title("Z Position at TS " + str(counter))
				plt.savefig("zdata/z_pos_" + str(counter) + ".png")
				plt.close()

				#print(counter)
				counter += 1
				continue

				particle_dict_entropy = self.find_nearest_neighbors(ts[1:], True)
				entropy_x = []
				entropy_y = []
				entropy_z = []
				entropy_e = []

				#go through and calculate entropy for each of these particles
				for particle, neighbors in particle_dict_entropy.iteritems():
					distance_sum = 0
					for n in neighbors:
						d = math.sqrt( (particle[0]-n[0])**2 + (particle[1]-n[1])**2 + (particle[2]-n[2])**2 )
						distance_sum += d
					distance = distance_sum / len(n)
					entropy = math.log(distance / self.unred_sig) / distance**2
					entropy_x.append(particle[0])
					entropy_y.append(particle[1])
					entropy_z.append(particle[2])
					entropy_e.append(entropy)


				#normalize the data
				#entropy_x_norm = (entropy_x - np.max(entropy_x))/-np.ptp(entropy_x)
				#entropy_y_norm = (entropy_y - np.max(entropy_y))/-np.ptp(entropy_y)
				#entropy_e_norm = (entropy_e - np.max(entropy_e))/-np.ptp(entropy_e)

				#save the image
				plt.figure(1, figsize=(13,13))
				plt.scatter(entropy_x, entropy_y, c=entropy_e, cmap='copper') #marker='.'
				#plt.imshow(np.array((entropy_x_norm, entropy_y_norm), dtype='uint8'))
				plt.title("Entropy at timestep " + str(counter))
				plt.colorbar()
				plt.savefig("entropytest5/entropy_" + str(counter) + ".png")
				plt.close()

			#print(counter)
			counter += 1

			#im.set_array(entropy_x_norm, entropy_y_norm, entropy_e_norm)
			#writer.grab_frame()


			'''
			#also kinda essential
		for sensor in self.pressure_sensor_list:
			sensor.generate_graph(self.fn)








			
			