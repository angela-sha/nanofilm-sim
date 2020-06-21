'''
Creates a Bin object, to make data analysis and storage easier when generating some
statistical information about the pressure data. It is meant to pull out the averages 
and quartiles, and then correct and re-interpret the data as such.
'''

import numpy as np
import matplotlib.pyplot as plt

class BinObject():

	'''Initializes the object.
	Parameters:
		min - the minimum area-value for the bin
		max - the maximum area-value for the bin
	Returns:
		null
	'''
	def __init__(self, min, max):
		self.min = min
		self.max = max
		self.tuples = []		#stores the tuples

		#value storage
		self.areas = []
		self.pressurex = []
		self.pressurey = []

		#parameter storage, separated beetween x and y
		self.avgx = 0			#average (mean)
		self.stdx = 0			#standard deviation
		self.q1x = 0			#quartile 1
		self.q3x = 0			#quartile 3

		self.avgy = 0
		self.stdy = 0
		self.q1y = 0
		self.q3y = 0

		#cutoffs for outliers
		self.high_range_x = 0
		self.high_range_y = 0
		self.low_range_x = 0
		self.low_range_y = 0



	'''Adds tuples of (area, pressure x, pressure y) to the object
	Parameters:
		t - a tuple in (area, px, py) form
	Returns:
		null
	'''
	def add_tuple(self, t):
		self.tuples.append(t)



	'''Computes mean, std, quartiles, outlier cutoffs.
	Parameters:
		none
	Returns:
		null, but updates object variables
	'''
	def do_stats(self):
		#first separate the variables
		for t in self.tuples:
			self.areas.append(t[0])
			self.pressurex.append(t[1])
			self.pressurey.append(t[2])

		px = np.array(self.pressurex)
		py = np.array(self.pressurey)

		#calculate the statistical values
		self.avgx = np.mean(px)
		self.stdx = np.std(px)
		self.avgy = np.mean(py)
		self.stdy = np.std(py)

		self.q1x = np.percentile(px, 25)
		self.q3x = np.percentile(px, 75)
		self.q1y = np.percentile(py, 25)
		self.q3y = np.percentile(py, 75)

		iqrx = self.q3x - self.q1x
		iqry = self.q3y - self.q1y

		self.high_range_x = self.q3x + 1.5 * iqrx
		self.high_range_y = self.q3y + 1.5 * iqry
		self.low_range_x = self.q1x - 1.5 * iqrx
		self.low_range_y = self.q1y - 1.5 * iqry



	'''Returns the averages and standard deviations.
	Parameters:
		none
	Returns:
		Averages and standard devations, as in the form below.
	'''
	def return_avgs(self):
		return( [(self.avgx, self.stdx), (self.avgy, self.stdy)] )



	'''Removes the outliers from the data set, and returns the cleansed data.
	Parameters:
		none
	Returns:
		Data in the form [ [(area), (pressure x)], [(area), (pressure y)] ]
	'''
	def remove_outliers(self):
		pressurex_graph = []
		pressurey_graph = []

		for val in range(len(self.areas)):
			if self.pressurex[val] < self.high_range_x and self.pressurex[val] > self.low_range_x and self.pressurex[val] > .1:
				pressurex_graph.append( (self.areas[val], self.pressurex[val]) )
			if self.pressurey[val] < self.high_range_y and self.pressurey[val] > self.low_range_y and self.pressurey[val] > .1:
				pressurey_graph.append( (self.areas[val], self.pressurey[val]) )

		return([pressurex_graph, pressurey_graph])







