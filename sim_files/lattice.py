#from sklearn.preprocessing import MinMaxScaler
	#use scaler when variables have different scales
#from sklearn.cluster import DBSCAN
	#clustering package
import numpy as np
import sys
#plotting inputs
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import griddata
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('seaborn-white')

from scipy.spatial import Voronoi
from scipy.spatial import cKDTree

class Lattice:
	'''Constructor for the Lattice object.
	Parameters:
	'''
	def __init__(self, pos_list):
		self.num_particles = 0
		self.tree = self.kd(pos_list)

	'''
	Given list of particle positions, creates kd_tree for lookup.
	Returns: kd tree (cKDTree)
	'''
	def kd(self, pos_list):
		return cKDTree(pos_list)

	'''
	Given list of particle positions, calculates nearest neighbors.
	Returns dictionary: particle position tuuple is key, list of tuples of nearest neighbors.
	Parameters: 
		pos_list: list of particle positions (list of tuples)
		n: number of neighbors to find
	'''
	def find_nearest_neighbors(self, pos_list, n):
		xyz_position = np.array([pos_list[i][0], pos_list[i][1], pos_list[i][2]] for i in range (len(positions)))

		n_dict = {}
		self.tree.query(pos_list, n)




