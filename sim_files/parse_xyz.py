import numpy as np
import sys

class Trajectory:
	'''
	Basic constructor for Trajectory class. 
	Parameters: file name.
	Returns: null
	'''
	def __init__(self, file_name, unred_params):
		self.fn = file_name
		self.num_particles = 0
		self.area_data = []
		self.pos_data = []

		self.u_eps = unred_params[0]
		self.u_sig = unred_params[1]
		self.u_mass = unred_params[2]

	'''
	Finds area of box given string of format:
	bounds lo (0, 0, 0) hi (x, y, z)
	'''
	def get_area(self, str):
		#get data
		bounds = str.strip().replace('(', '*').replace(')', '*').split('*')
		bounds_lo = bounds[1].split(',')
		bounds_hi = bounds[3].split(',')

		#change to float
		bounds_lo = [float(i) for i in bounds_lo]
		bounds_hi = [float(j) for j in bounds_hi]

		#calculate area
		return abs((bounds_hi[0] - bounds_lo[0]) * (bounds_hi[1]-bounds_lo[1]) * (self.u_sig ** 2) / (100* float(self.num_particles)))

	'''
	Converts .xyz file into simulation_data list.
	Updates: area_data and pos_data (indexes are same)
	Returns: number of time steps
	timestep_data is in format [(bounds)]
	'''
	def parse_file(self):
		input_file = self.fn
		file = open(input_file, 'r')

		#store number of particles
		self.num_particles = int(file.readline().strip())

		timestep_data = []
		num_timesteps = 0

		new_timestep = False
		first = True
		
		for line in file:
			if first:
				first = False
				continue
			if line.strip() == str(self.num_particles):
				new_timestep = True
				self.pos_data.append(timestep_data)
				continue
			if new_timestep:
				timestep_data = []
				new_timestep = False
				num_timesteps+=1
				print("Finished timestep " + str(num_timesteps))

				#add new area of box to new timestep array
				self.area_data.append(self.get_area(line))
			else:		#otherwise, xyz stored in tuple and appended
				line = line.strip().split()
				line = line[1:]
				line = [float(i) for i in line]
				timestep_data.append(line)

		#append final timestep
		self.pos_data.append(timestep_data) #append last timestep
		num_timesteps += 1
		print("Finished timestep " + str(num_timesteps))
		
		self.area_data = np.array(self.area_data) 
		self.pos_data = np.array(self.pos_data)
		return num_timesteps

	'''
	Returns area data.
	'''
	def get_area_data(self):
		return self.area_data

	'''
	Returns position data.
	'''
	def get_pos_data(self):
		return self.pos_data

	'''
	Returns number of particles.
	'''
	def get_num_particles(self):
		return self.num_particles

	'''
	Prints area.txt file in format area\narea\n etc. and writes to file.
	'''
	def write_area(self): 
		area = open("area.txt", "w")
		for i in self.area_data:
			area.write(str(i)+"\n")
		area.close()


