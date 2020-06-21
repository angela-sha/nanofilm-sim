#from sklearn.preprocessing import MinMaxScaler
	#use scaler when variables have different scales
from sklearn.cluster import DBSCAN
from sklearn.neighbors import KernelDensity
from sklearn import metrics
	#clustering package
import numpy as np
import sys
#plotting inputs
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import griddata
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm

#style customizations for matplotlib
import seaborn as sns
sns.set()
#set macros for text size
SMALL_SIZE = 8
BIG_SIZE = 10
plt.rc('font', size=SMALL_SIZE)         #default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)    #fontsize of the axes title
plt.rc('axes', labelsize=SMALL_SIZE)    #fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)   #fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)   #fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)   #legend fontsize
plt.rc('figure', titlesize=BIG_SIZE)  	#fontsize of the figure title

class Coordinate:
	''' Constructor for coordinate object.
	Parameters: area and position data lists.
	Returns: null.
	'''
	def __init__(self, num_particles, num_ts, area_data, pos_data, unred_params):
		self.num_particles = num_particles
		self.num_ts = num_ts
		self.area_data = area_data
		self.pos_data = pos_data

		self.u_eps = unred_params[0]
		self.u_sig = unred_params[1]
		self.u_mass = unred_params[2]

	'''
	Returns the coord_data of one timestep, given index of coordinate (0-2).
	Prerequisite: timestep contains at least one particle.
	'''
	def get_coord_data(self, ts, idx):
		ts_data = self.pos_data[ts]
		return ts_data[:,idx]

	'''
	Calculate IQR and max value not counted as outlier.
	'''
	def outlier_range(self, z_data):
		q1 = np.quantile(z_data, .25)
		q3 = np.quantile(z_data, .75)
		iqr = q3-q1
		return q3+1.5*iqr

	'''
	Calculate threshold value, which is defined by a percentage between the given 
	"layers" (0 being the lower layer, 100 being the upper layer) 
	'''
	def get_threshold_value(self, ts, percentile):
		#first separate the z-values into two layers
		vals = self.pos_data[ts]
		avg = np.mean(vals, axis=0)
		avg_z = avg[2]
		return avg_z

	'''
	Calculate monolayer threshold value. Extracted from KDE graph, given a timestep.
	
	def get_monolayer_threshold(self, ts, percentile):
		max = find_KDE_max(self, ts);
		#first separate the z-values into two layers
		vals = self.pos_data[ts]
		avg = np.mean(vals, axis=0)
		avg_z = avg[2]
		return avg_z
	'''
	#find maximum first

	'''
	Given a particle index, checks if contained within the given layer. 
	Returns a boolean value. 
	Note: inputted value in nm.
	'''
	def check_layer(self, z, num):
		#min_val = self.get_threshold_value(818, 50)
		if num == 1:
			min_val = 322.5
			max_val = 327.5
		if num == 1.5:
			min_val = 327.5
			max_val = 330
		if num == 2:
			min_val = 330
			max_val = 335
		#if self.nm_to_red(z) >= min_val and self.nm_red(z) <= max_val:
		if z >= min_val and z <= max_val:
			return True
		else:
			return False

	'''
	Returns reduced value from nm.
	'''
	def nm_to_red(self, val):
		return val/self.u_sig*10.

	'''
	Returns nm value from reduced value.
	'''
	def red_to_nm(self, val):
		return val*self.u_sig/10.

	'''
	Returns positions of a given particle from a given timestep to end of sim.
	In form [[x,y,z], [...], ...]
	Note: if range is set to 0, this indicates to end.
	'''
	def track_particle(self, idx, ts, ts_range):
		if range == 0:
			post_ts_data = self.pos_data[ts:]
		else:
			post_ts_data = self.pos_data[ts:ts+ts_range]
		return post_ts_data[:,idx]

	'''
	Graphs 2-dimensional graph of one particle's positions. 
	Input parameters: index of particle, timestep to start at.
	'''
	def graph_particle_path(self, idx, ts, ts_range, fn):
		path = self.track_particle(idx, ts, ts_range)
		path = path*self.u_sig/10.
		x_data = path[:,0]
		y_data = path[:,1]
		z_data = path[:,2]
		fig = plt.figure()
		#ax = fig.add_subplot(1, 1, 1)
		plt.xlabel('x position (nm)', labelpad=5)
		plt.ylabel('y position (nm)', labelpad=10) 
		#plt.plot(x_data, y_data, '-', linewidth=.5)
		
		c = cm.jet((z_data-np.min(z_data))/(np.max(z_data)-np.min(z_data)))
		ax = plt.gca()
		for i in np.arange(len(x_data)-1):
			ax.plot([x_data[i],x_data[i+1]], [y_data[i], y_data[i+1]], c=c[i], lw=.5)
		#plt.colorbar()
		plt.savefig(fn, dpi=600)

		plt.close()

	'''
	Graphs 2-dimensional graph of positions given an x-y boundary.
	Note: particles are to-scale in this representation.
	Note: grid is a boolean value for whether or not there is a grid behind.
	Note: input values are unreduced.
	'''
	def graph_2D(self, ts, x_lo, y_lo, size, grid, fn):
		vals = self.pos_data[ts]
		vals = vals*self.u_sig/10.
		plot = []
		plot_2 = []
		plot_3 = []
		#height = []

		x_hi = x_lo+size
		y_hi = y_lo+size
		
		for val in vals:
			if val[0] > x_lo and val[0] < x_hi and val[1] > y_lo and val[1] < y_hi:
				if self.check_layer(val[2], 1):
					plot.append((val[0], val[1]))
				elif self.check_layer(val[2], 2):
					plot_2.append((val[0], val[1]))
				elif self.check_layer(val[2], 1.5):
					plot_3.append((val[0], val[1]))
				#height.append(val[2])

		#calculate the threshold value used for the plots
		#print self.get_threshold_value(818, 50)
		fig = plt.figure()
		ax = fig.add_subplot(1, 1, 1)
		ax.set_xlim(x_lo, x_hi)
		ax.set_ylim(y_lo, y_hi)
		#scale axes (convert to nm)
		#ticks = ticker.FuncFormatter(lambda a, pos: '{0:g}'.format(self.red_to_nm(a)))
		#ax.xaxis.set_major_formatter(ticks)
		#ax.yaxis.set_major_formatter(ticks)
		ax.set_facecolor('w')
		if grid:
			ax.grid(color='k', lw=.5)
		plt.xlabel('x position (nm)', labelpad=5)
		plt.ylabel('y position (nm)', labelpad=10) 
		ax.set_aspect(1)			#same scaling on x-y axes

		r = 2.5#self.nm_to_red(2.5) #convert to reduced radius of the nanoparticle
		circles = [plt.Circle(pi, radius=r, lw=None) for pi in plot]
		circles_2 = [plt.Circle(pi, radius=r, lw=None) for pi in plot_2]
		circles_3 = [plt.Circle(pi, radius=r, lw=None) for pi in plot_3]
		c = matplotlib.collections.PatchCollection(circles, edgecolor='b')
		c2 = matplotlib.collections.PatchCollection(circles_2, color='r', edgecolor='r')
		c3 = matplotlib.collections.PatchCollection(circles_3, color='m', edgecolor='m')
		ax.add_collection(c)
		ax.add_collection(c3)
		ax.add_collection(c2)

		plt.savefig(fn, dpi=600)

		plt.close()

	'''
	Graphs the position distribution of the particles given the timestep index. 
	'''
	def graph_pos_dist(self, ts, fn):
		#run the file
		#eqlm_ts = 0 #subtract one for index

		fig = plt.figure()
		ax = fig.add_subplot(1, 1, 1, projection='3d')

		temp_x = self.get_coord_data(ts, 0)
		temp_y = self.get_coord_data(ts, 1)
		temp_z = self.get_coord_data(ts, 2)

		temp_x = temp_x*self.u_sig/10.
		temp_y = temp_y*self.u_sig/10.
		temp_z = temp_z*self.u_sig/10.

		ax.scatter(temp_x, temp_y, temp_z, s=.3,c='r', marker='o', alpha =.5)
		ax.set_xlabel('x position (nm)')
		ax.set_ylabel('y position (nm)')
		ax.set_zlabel('z position (nm)')

		#note:  limits are manually set from visual data. comment out for autoscaling
		ax.set_xlim([0,650])
		ax.set_ylim([0,1000])
		ax.set_zlim([320,334])

		plt.savefig(fn, dpi=600)

		plt.close()

	'''
	Clusters univariate (z-position) data given timestep, using kernel density estimation.
	Takes a boolean that defines whether or not a graph output will be produced.
	If graph is not outputted, returns list of minimums of 
	'''
	def cluster_kde(self, ts, graph, fn):
		#define the z-values and bins
		z = self.get_coord_data(ts, 2)
		z_d = np.linspace(54.60, 57.00, 240)

		#instantiate/fit the kde model
		kde = KernelDensity(bandwidth=0.1, kernel='gaussian')
		kde.fit(z[:, None])
		# score_samples returns the log of the probability density
		logprob = kde.score_samples(z_d[:, None])

		if graph:
			fig = plt.figure()
			plt.plot(z_d[:, None], np.exp(logprob), '-', linewidth=1)
			plt.plot(z, np.full_like(z, -0.01), '|k', markeredgewidth=.75)
			plt.savefig(fn, dpi=600)
			plt.close()

	'''
	Clusters given timestep using DBSCAN. 
	'''
	def cluster_DBSCAN(self, ts):
		temp_z.reshape(-1, 1)
		db = DBSCAN(eps = 1.5, min_samples = 25, metric='euclidean',
			algorithm='kd_tree')
		db.fit(temp_z)
		db.fit(self.pos_data[ts])
		core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
		core_samples_mask[db.core_sample_indices_] = True
		labels = db.labels_

		# Number of clusters in labels, ignoring noise if present.
		n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
		n_noise_ = list(labels).count(-1)

		print('Estimated number of clusters: %d' % n_clusters_)
		print('Estimated number of noise points: %d' % n_noise_)
		#plot the corresponding graph

	#-----------------------------FUNCTIONS USING KDE-----------------------------

	'''
	Find number of minimums, and minimum values (x-value of graph, z-value of 
	simulation) of KDE, given timestep.
	
	def find_KDE_min(self, ts):
		self.cluster_kde
	
	Find number of maximums, and maximum values (x-value of graph, z-value of 
	simulation) of KDE, given timestep
	def 
	Returns the percentile of particles in the "first layer" using KDE, given timestep.
	
	def calculate_percentiles(self, ts):
	'''
