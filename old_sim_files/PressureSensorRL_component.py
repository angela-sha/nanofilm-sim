'''
This pressure sensor is meant to generate things component by component, but is fundamentally the same as PressureSensorRL.py
'''

from PressureSensor import PressureSensor
from Grapher import Grapher
import math


class PressureSensorRL_Component(PressureSensor):

	'''Creates the sensor itself, only defining its location and bounds.
	Parameters: 
		bounds - tuple of (x1, x2, y1, y2)
		frame - float of how large around the frame you can calculate potentials
		params - the parameters that are not sigma, epsilon, as a tuple
		unreduced_params - (epsilon, sigma, mass) before they are converted to be one, used to calculate real values for the isotherms generated
	Returns:
		null
	'''
	def __init__(self, bounds, frame, params, unreduced_params):
		PressureSensor.__init__(self, bounds, frame, unreduced_params)
		self.alpha = params[0]
		self.r_a = params[1]
		self.a_g = params[2]
		self.sig_g = params[3]
		self.r_g = params[4]


	'''Calculates the force, using an RL potential (NOT LJ).
	Parameters:
		r - distance
	Returns:
		force, in real units
	'''
	def force_calc(self, r):

		if(r != float(0) ):
			p1 = 64 * math.pow(r, -65)
			p2 = 1 / self.alpha * math.exp(- (r - self.r_a) / self.alpha)
			p3 = self.a_g / (self.sig_g * self.sig_g) * (r - self.r_g) * math.exp(-math.pow(r - self.r_g, 2) / (2 * self.sig_g * self.sig_g))

			force = self.eps * (p1 - p2 + p3)
			return force

		return 0.0

	'''Redoes the pressure calculation so that it is done component by component, rather than all at once
	'''
	def calc_pressure(self, nn_dict, area_of_simulation):
		list = [i for i in nn_dict.iterkeys()]
		plate_particles = self.find_particles_in(list, self.bounds)

		pressure_sum = [0, 0, 0]

		for particle in plate_particles:
			neighbors = nn_dict[particle]

			for n in neighbors:
				#r = math.sqrt( ((n[0] - px) * (n[0] - px)) + ((n[1] - py) * (n[1] - py))  + ((n[2] - pz) * (n[2] - pz)) )
				forcex = self.force_calc(particle[0])
				forcey = self.force_calc(particle[1])
				forcez = self.force_calc(particle[2])
				pressurex = forcex * particle[0] #* self.ur_eps / ((self.ur_sig * 1e-9) ** 3)		#convert reduced units -> Pascals
				pressurey = forcey * particle[1]
				pressurez = forcez * particle[2]
				pressure_sum[0] += pressurex
				pressure_sum[1] += pressurey
				pressure_sum[2] += pressurez

		pressure_sum = [i / 2. for i in pressure_sum]
		self.pressures.append([i / self.plate_area for i in pressure_sum])
		self.areas.append(area_of_simulation)


	def generate_graph(self, fn):
		g = Grapher(self.areas, self.pressures, "component")
		g.graph("Attempt", fn)


