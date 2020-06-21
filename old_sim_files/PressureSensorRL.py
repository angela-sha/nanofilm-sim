'''
Subclass of PressureSensor, tailored to use the Rice-Lin potential function instead.
'''

from PressureSensor import PressureSensor
import math


class PressureSensorRL(PressureSensor):

	'''Creates the sensor itself, only defining its location and bounds.
	Parameters: 
		bounds - tuple of (x1, x2, y1, y2)
		frame - float of how large around the frame you can calculate potentials
		params - the parameters that are not sigma, epsilon, as a tuple
		unreduced_params - (epsilon, sigma, mass) before they are converted to be one, used to calculate real values for the isotherms generated
	Returns:
		null
	'''
	def __init__(self, bounds, frame, params, unreduced_params, g):
		PressureSensor.__init__(self, bounds, frame, unreduced_params, g)
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

	def potential_calc(self, r):
		if(r != float(0) ):
			p1 = math.pow((self.sig/r),64)
			p2 = math.exp(- (r - self.r_a) / self.alpha)
			p3 = self.a_g * math.exp( -math.pow(r - self.r_g, 2) / (2 * self.sig_g * self.sig_g) )

			pot = self.eps * (p1 - p2 + p3)
			return pot

		return 0.0


