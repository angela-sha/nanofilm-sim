import math
import random
import sys
#sys.path = sys.path + ['/home/alexfeistritzer/dash_attempt/md_engine/build/python/build/lib.linux-x86_64-2.7']
sys.path = sys.path + ['/home/cmliepold/Dans_GPU/md_engine-master7/build_tmp/python/build/lib.linux-x86_64-2.7']
from Sim import *
import numpy as np

class Simulator: 

	''' Basic constructor for simulator (LJ simulation)
		Parameters:
			boxparam: (length, width, height) in reduced units
			density: packing density (nonreduced units)
			temp: temperature (reduced units)
			dr: deformation rate
			multiplier: mulitplier for deformation (Vector(x,y,z))
			params: (epsilon, sigma) tuple
		Returns:
			null
	'''
	def __init__(self, boxparam, density, temp, dr, multiplier, params, command_list):
		self.length = boxparam[0]
		self.width = boxparam[1]
		self.height = boxparam[2]

		self.pack_len = self.calc_pack_len(density)
		self.temp = temp
		self.fn = ""				#the name of the xyz file
		self.deform_multiplier = multiplier
		self.deform_rate = dr
		self.command_list = command_list 

		self.eps = params[0]
		self.sig = params[1]

	'''Calculates the packing density in reduced units from the reduced density.
	Parameters: 
		d - reduced density
	Returns:
		particle spacing in reduced units
	'''
	def calc_pack_len(self, d):

		#comes from l = sqrt(2/ [rho~ * sqrt(3)] )
		l = math.sqrt(2. / (d* math.sqrt(3)))
		return l

	'''Updates the filename for the iteration
	Parameters: 
		str - the filename
	Returns: 
		null
	'''
	def update_fn(self, str):
		self.fn = str

	'''Sets the pairwise interaction potential fix; use this to update in subclasses.
	Parameters: 
		state - state
		handle - name of the fix used in the simulation
	Returns: 
		the fix as an object, passed into the simulation
	'''
	def set_pair_interaction_fix(self, state, handle):
		return FixLJCut(state, handle)


	'''Sets all of the parameters necessary to run the simulation, based entirely on the fixpair from set_pair_interaction_fix. Update in the subclasses.
	Parameters: 
		nonbond - the name of the fix that set_pair_interaction_fix returns
	Returns:
		null
	Note: 
		sets all of the parameters individually
	'''
	def set_parameters(self, nonbond):
		nonbond.setParameter('sig', 'spc1', 'spc1', self.sig)
		nonbond.setParameter('eps', 'spc1', 'spc1', self.eps)

	'''Runs the simulation.
	Parameters:
		none
	Returns:
		null
	'''
	def run(self): 

		#defines simulation box
		state = State()
		state.deviceManager.setDevice(0)
		state.bounds = Bounds(state, lo = Vector(0, 0, 0), hi = Vector(self.length, self.width, self.height))
		state.rCut = 2.5
		state.padding = .6
		state.periodicInterval = 7
		state.shoutEvery = 100

		#add particles into the simulation
		state.atomParams.addSpecies(handle='spc1', mass=1, atomicNum=1)
		nonbond = self.set_pair_interaction_fix(state, 'cut')
		self.set_parameters(nonbond)

		state.activateFix(nonbond)

		#place atoms inside simulation box
		num_particles_length = int(.8 * self.length / self.pack_len)
		num_particles_width = int(.8 * self.width / self.pack_len)
		x = 0.0	
		y = 0.0
		z = self.height/2.

		hang_length = self.length - (num_particles_length * self.pack_len)	#space left over in row
		hang_width = self.width - (num_particles_width * self.pack_len)		#space left over in column
		offrow_length = 0
		offrow_width = 0

		if(hang_length >= self.pack_len/2):			#checks to see if a half-length can fit in overhang
			offrow_length = num_particles_length		
		else:
			offrow_length = num_particles_length - 1	#if hang is too small, reduces the amount of particles in the off rows

		if(hang_width >= self.pack_len/2):			#repeat with columns
			offrow_width = num_particles_width		
		else:
			offrow_width = num_particles_width - 1

		#fill box with particles
		for i in range(num_particles_width):
			if i%2 == 0:
				for j in range(num_particles_length):
					state.addAtom('spc1', Vector(float(x),float(y),float(z)), 0.0)
					x += self.pack_len
			else:
				x += self.pack_len / 2.
				for j in range(offrow_length):
					state.addAtom('spc1', Vector(x,y,z), 0.0)
					x += self.pack_len
			y+= self.pack_len
			x = 0.0

		
		#initializing fixes
		InitializeAtoms.initTemp(state, 'all', self.temp)
		fixNVT = FixNoseHoover(state, 'temp', 'all', self.temp, timeConstant=1)
		state.activateFix(fixNVT)

		deform = FixDeform(state, handle='def', groupHandle='all', deformRate=self.deform_rate, multiplier=self.deform_multiplier, applyEvery=1) 
		wallDist = self.height/2.
		#harmonic potential to keep particles in place
		topWall = FixWallHarmonic(state, handle="topwall", groupHandle='all', origin=Vector(0, 0, state.bounds.hi[2]), forceDir=Vector(0, 0, -1),dist=wallDist, k=100)
		bottomWall = FixWallHarmonic(state, handle="bottomwall", groupHandle='all',origin=Vector(0, 0,state.bounds.lo[2]),forceDir=Vector(0, 0, 1), dist=wallDist, k=1000)
		state.activateFix(topWall)
		state.activateFix(bottomWall)

		pressure = FixPressureBerendsen(state, 'constP', .2, 10, 1)
		integVerlet = IntegratorVerlet(state)

		#print statements
		writeconfig = WriteConfig(state, fn=str(self.fn), writeEvery=100, format='xyz', handle='writer')
		state.activateWriteConfig(writeconfig)

		writeconfig.writeEvery = self.command_list[0][0]
		integVerlet.run(self.command_list[0][1])

		state.activateFix(deform)
		writeconfig.writeEvery = self.command_list[1][0]
		integVerlet.run(self.command_list[1][1])

