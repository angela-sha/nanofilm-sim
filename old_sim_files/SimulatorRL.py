'''
This Simulator subclass is tailored for use for the Rice-Line potential, created by Chris Liepold. It only changes the pairwise interaction and the parameters used for this interaction.
'''

import math
import sys
sys.path = sys.path + ['/home/cmliepold/Dans_GPU/md_engine-master7/build/python/build/lib.linux-x86_64-2.7']
from Sim import *
from simulator import Simulator


class SimulatorRL(Simulator):

	'''Adding extra parameters to the initialization.
	Parameters:
		boxparam - box parameters (3 tuple, (length, width, height) ), in reduced units
		density - packing denisty (not in reduced units)
		temp - temperature (in reduced units) 
		dr - deformation rate
		params - tuple:
			alpha; r_a; A_g; sigma_g; r_g
	'''
	def __init__(self, boxparam, density, temp, dr, multiplier, params, command_list):
		Simulator.__init__(self, boxparam, density, temp, dr, multiplier, params[:2], command_list)
		self.alpha = params[2]
		self.r_a = params[3]
		self.a_g = params[4]
		self.sig_g = params[5]
		self.r_g = params[6]

	'''Set the pairwise interaction fix to be RLCut
	'''
	def set_pair_interaction_fix(self, state, handle):
		return FixRLCut(state, handle)


	'''Add the extra parameters to the fix.
	'''
	def set_parameters(self, nonbond):
		nonbond.setParameter('sig', 'spc1', 'spc1', self.sig)
		nonbond.setParameter('eps', 'spc1', 'spc1', self.eps)
		nonbond.setParameter('alpha', 'spc1', 'spc1', self.alpha)
		nonbond.setParameter('r_a', 'spc1', 'spc1', self.r_a)
		nonbond.setParameter('a_g', 'spc1', 'spc1', self.a_g)
		nonbond.setParameter('sig_g', 'spc1', 'spc1', self.sig_g)
		nonbond.setParameter('r_g', 'spc1', 'spc1', self.r_g)


