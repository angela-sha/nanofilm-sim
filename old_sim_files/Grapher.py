'''
Creates a graphing object used for graphing the results of the pressure sensor.
Only meant to store a single simulation run before returning a graph, and not meant to hold the results of multiple simulations.

Can be inherited or modified to output different types of graphs describing different data; in this case, it is used to print out the pressure - area isotherms later used to calculate the strain/stress moduli.
'''

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
from BinObject import BinObject

class Grapher:
	'''Constructor for the grapher object.
	Parameters: 
		tag - string that describes what kind of graphing that will be used; options:
			"best_fit" - uses a quintic interpolation from numpy to generate a line that passes through the data
			"connect" - basic line graph, no fitting function
	Returns:
		null
	'''
	def __init__(self, tag):

		#for a single timestep
		self.area = []
		self.pressure = []
		self.tag = tag

		#for all of the time steps
		self.area_holder = []
		self.pressurex_holder = []
		self.pressurey_holder = []


	'''Used to update the area for a single simulation, so as to reuse the same grapher for each sensor.
	Parameters:
		a - an array of areas, already unreduced and in units of nm^2/particle
	Returns:
		null
	'''
	def update_area(self, a):
		self.area = a

	'''Used to update the area for a single simulation, so as to reuse the same grapher for each sensor.
	Parameters:
		p - an array of pressures, already unreduced, and in the form [ [px, py], [px, py], ... ] where each pressure comes from a different timestep
	Returns:
		null
	'''
	def update_pressure(self, p):
		self.pressure = p

	'''Returns the area_holder array; used when making the new graphs from the sensor object
	Paramers:
		none
	Returns:
		null
	'''
	def return_area(self):
		return self.area_holder

	'''Returns the pressurex_holder array; used when making the new graphs from the sensor object
	Paramers:
		none
	Returns:
		null
	'''
	def return_pressurex(self):
		return self.pressurex_holder

	'''Returns the pressurey_holder array; used when making the new graphs from the sensor object
	Paramers:
		none
	Returns:
		null
	'''
	def return_pressurey(self):
		return self.pressurey_holder

	'''Generates a graph to save to the current working directory.
	Parameters:
		title - string of the title of the graph
		fn - the name of the graph when it is saved as a .png
	Returns:
		null, but generates a graph
	'''
	def graph(self, title, fn):
		#generates a best fit (actually quintic spline) line from the data
		if self.tag == "just_plot_the_goddamn_graph":
			fig = plt.figure()
			ax = fig.add_subplot(1,1,1)

			x = []
			y = []

			for p in self.pressure:
				x.append(p[0])
				self.pressurex_holder.append(p[0])
				y.append(p[1])
				self.pressurey_holder.append(p[1])

			for a in self.area:
				self.area_holder.append(a)

			ax.plot(self.area, x, 'r.', self.area, y, 'b.')
			#ax.plot(self.area, x, 'k', linewidth=1)
			#ax.plot(self.area, y, 'g', linewidth=1)
			ax.set_title("Potential Values without Lines")
			plt.savefig(fn)
			plt.close()

		elif self.tag == "best_fit":
			z = np.polyfit(self.area, self.pressure, 5)
			f = np.poly1d(z)
			x = np.linspace(self.area[0], self.area[-1], 1000)
			y = f(x)

			plt.figure(1)
			plt.plot(self.area, self.pressure, 'o', x, y, 'b')

			plt.xlabel("Area / nanoparticle")
			plt.ylabel("Pressure")

			plt.title(title)

			plt.savefig(fn)
			plt.close()

		#line graph, point by point
		elif self.tag == "connect":
			plt.figure(1)
			plt.plot(self.area, self.pressure, 'k')
			plt.xlabel("Area / nanoparticle")
			plt.ylabel("Pressure")
			plt.title(title)

			plt.savefig(fn)
			plt.close()

		#makes the graph component by component
		elif self.tag == "component":
			print('Beginning graphing!')
			#parse out the data
			x = []
			y = []

			for p in self.pressure:
				x.append(p[0])
				self.pressurex_holder.append(p[0])
				y.append(p[1])
				self.pressurey_holder.append(p[1])

			for a in self.area:
				self.area_holder.append(a)

			#print(len(self.area))
			#print(len(x))
			#print(len(y))

			#things up to here are functioning fine

			#generate the plots
			f, (gx, gy) = plt.subplots(2)
			gx.plot(self.area, x, 'ro')
			gx.legend(["x"], fontsize='xx-small')
			gy.plot(self.area, y, 'bo')
			gy.legend(["y"], fontsize='xx-small')

			gx.set_title(title)

			plt.savefig(fn)
			plt.close()

			print("Completed part 1")

			f = open("area.txt", "w")
			g = open("x.txt", "w")
			h = open("y.txt", "w")
			for i in range(len(self.area)):
				f.write(str(self.area[i]) + "\n")
				g.write(str(x[i]) + "\n")
				h.write(str(y[i]) + "\n")
			f.close()
			g.close()
			h.close()

			'''
			#hopefully this works to graph stuffs
			plt.figure(1)
			self.area.reverse()
			x.reverse()
			y.reverse()

			#print(len(self.area)) ->1000

			#5E5 -> way too much
			#5E4 <- not enough
			tck_x = interpolate.splrep(self.area, x, s=1E4)
			tck_y = interpolate.splrep(self.area, y, s=1E4)
			#print(self.area[0])
			#print(self.area[-1])
			area_new = np.linspace(self.area[0], self.area[-1], len(self.area)*10)

			#print(area_new)

			x_new = interpolate.splev(area_new, tck_x, der=0)
			y_new = interpolate.splev(area_new, tck_y, der=0)

			#print(x_new)
			#print(y_new)

			plt.plot(area_new, x_new, area_new, y_new)
			plt.legend(['Perpendicular', 'Parallel'])
			plt.title("Isotherm")
			plt.savefig(fn + "_isotherm")
			plt.close()

			x_der = interpolate.splev(area_new, tck_x, der=1)
			y_der = interpolate.splev(area_new, tck_y, der=1)
			plt.figure(2)
			k = -area_new / 2 * (y_der + x_der)
			g = -area_new / 2 * (y_der - x_der)
			plt.plot(area_new, k, area_new, g)
			plt.legend(['Compressive Modulus, K', "Shear Modulus, G"])
			plt.title("Moduli Plot")
			plt.savefig(fn + "_moduli_plot")
			plt.close()
			'''

		#I guess like a moving average? Averages some number of data points and then uses that number
		elif self.tag == "binhua_method":
			print('Beginning graphing!')
			#parse out the data
			x = []
			y = []

			for p in self.pressure:
				x.append(p[0])
				self.pressurex_holder.append(p[0])
				y.append(p[1])
				self.pressurey_holder.append(p[1])

			for a in self.area:
				self.area_holder.append(a)

			#This is where we're gonna try and do the averaging
			x_avgd = []
			y_avgd = []
			a_avgd = []

			sumx = 0
			sumy = 0
			suma = 0
			counter = 1
			average_bin = 10
			for i in range(len(x)):
				if counter % average_bin == 0:
					ax = sumx/average_bin
					ay = sumy/average_bin
					aa = suma/average_bin

					if( (ax not in x_avgd) and (ay not in y_avgd) and (aa not in a_avgd)):
						x_avgd.append(ax)
						y_avgd.append(ay)
						a_avgd.append(aa)
					
					sumx = x[i]
					sumy = y[i]
					suma = self.area[i]

				else:
					sumx += x[i]
					sumy += y[i]
					suma += self.area[i]

				counter+=1
				

			x = x_avgd[1:]
			y = y_avgd[1:]
			self.area = a_avgd[1:]


			print(self.area)
			#things up to here are functioning fine

			#generate the plots
			f, (gx, gy) = plt.subplots(2)
			gx.plot(self.area, x, 'ro')
			gx.legend(["x"], fontsize='xx-small')
			gy.plot(self.area, y, 'bo')
			gy.legend(["y"], fontsize='xx-small')

			gx.set_title(title)

			plt.savefig(fn)
			plt.close()

			print("Completed part 1")

			#hopefully this works to graph stuffs
			plt.figure(1)
			self.area.reverse()
			x.reverse()
			y.reverse()

			#print(len(self.area)) ->1000

			#5E5 -> way too much
			#5E4 <- not enough
			tck_x = interpolate.splrep(x=self.area, y=x, s=666)
			tck_y = interpolate.splrep(x=self.area, y=y, s=666) #1E4
			#print(self.area[0])
			#print(self.area[-1])
			area_new = np.linspace(self.area[0], self.area[-1], len(self.area)*10)

			#print(area_new)

			x_new = interpolate.splev(area_new, tck_x, der=0)
			y_new = interpolate.splev(area_new, tck_y, der=0)

			#print(x_new)
			#print(y_new)

			plt.plot(area_new, x_new, area_new, y_new)
			plt.legend(['Perpendicular', 'Parallel'])
			plt.title("Isotherm")
			plt.savefig(fn + "_isotherm")
			plt.close()

			x_der = interpolate.splev(area_new, tck_x, der=1)
			y_der = interpolate.splev(area_new, tck_y, der=1)
			plt.figure(2)
			k = -area_new / 2 * (y_der + x_der)
			g = -area_new / 2 * (y_der - x_der)
			plt.plot(area_new, k, area_new, g)
			plt.legend(['Compressive Modulus, K', "Shear Modulus, G"])
			plt.title("Moduli Plot")
			plt.savefig(fn + "_moduli_plot")
			plt.close()
		
		#attempting to use a convolution (moving average) to best graph the data
		elif self.tag == "convolution":	
			print('Beginning graphing!')
			#parse out the data
			x = []
			y = []

			for p in self.pressure:
				x.append(p[0])
				self.pressurex_holder.append(p[0])
				y.append(p[1])
				self.pressurey_holder.append(p[1])

			for a in self.area:
				self.area_holder.append(a)


			#this is the bit that was added
			#####################################################
			f = open("area.txt", "w")
			g = open("pres_x.txt", "w")
			h = open("pres_y.txt", "w")

			for i in range(len(x)):
				f.write(str(self.area[i]) + "\n")
				g.write(str(x[i])+"\n")
				h.write(str(y[i])+"\n")

			f.close()
			g.close()
			h.close()

			exit()

			#######################################################

			#print(len(self.area))
			#print(len(x))
			#print(len(y))

			#things up to here are functioning fine

			#generate the plots
			f, (gx, gy) = plt.subplots(2)
			gx.plot(self.area, x, 'r.', self.area, x, 'k')
			gx.legend(["x"], fontsize='xx-small')
			gy.plot(self.area, y, 'b.', self.area, y, 'y')
			gy.legend(["y"], fontsize='xx-small')

			gx.set_title(title)

			plt.savefig(fn)
			plt.close()

			print("Completed part 1")

			#hopefully this works to graph stuffs
			plt.figure(1)
			#self.area.reverse()
			x.reverse()
			y.reverse()

			plt.plot(np.convolve(self.area, x, mode='same'))
			plt.plot(np.convolve(self.area, y, mode='same'))

			plt.legend(['Perpendicular', 'Parallel'])
			plt.title("Isotherm")
			plt.savefig(fn + "_isotherm")
			plt.close()

			plt.figure(1)
			plt.plot(np.convolve(x, np.ones(len(x),)/len(x), mode='same'))
			plt.plot(np.convolve(y, np.ones(len(y),)/len(y), mode='same'))
			plt.legend(['Perpendicular', 'Parallel'])
			plt.title("Isotherm #2")
			plt.savefig(fn + "_isotherm2")
			plt.close()

	'''Generates a graph of all of the data, from all simulations that used this sensor/grapher combo.
	Parameters:
		none
	Returns:
		null, but generates a graph
	'''
	def full_graph(self):
		#generate the best fit lines
		z1 = np.polyfit(self.area_holder, self.pressurex_holder, 5)
		z2 = np.polyfit(self.area_holder, self.pressurey_holder, 5)
		f1 = np.poly1d(z1)
		f2 = np.poly1d(z2)
		x = np.linspace(min(self.area_holder), max(self.area_holder), 1000)
		y1 = f1(x)
		y2 = f2(x)

		#generate the plot
		plt.figure(1)
		plt.plot(self.area_holder, self.pressurex_holder, 'ro', self.area_holder, self.pressurey_holder, 'bo', x, y1, 'k', x, y2, 'y')

		#labels and saving
		plt.xlabel("Area / nanoparticle")
		plt.ylabel("Pressure")

		plt.title("Combined Graph")

		plt.savefig("combined_graph")
		plt.close()

		#do the statistical analysis
		self.stat_analysis()

	'''Performs the statistical analysis, as follows:
		-> divides the data into 100 bins
			-> finds min, max, average, std, and quartiles
		-> removes outliers from the data, using the standard IQR * 1.5 method
		-> generates a new graph of just averages and standard devations
		-> generates a graph similar to full_graph(), but with the removal of outliers
	Parameters:
		none
	Returns:
		null, but generates two graphs
	'''
	def stat_analysis(self):
		#store data in tuples so it doesn't get messed with too much
		tuple_holder = []
		for val in range(len(self.area_holder)):
			t = (self.area_holder[val], self.pressurex_holder[val], self.pressurey_holder[val])
			tuple_holder.append(t)

		#determining bin sizes
		min_val = min(self.area_holder)
		max_val = max(self.area_holder)

		bin_size = float(max_val - min_val) / 100.

		#stores the bin objects created
		bin_holder = []

		#generates bin objects using bin size, and adds all of the tuples made above to the appropriate bin
		for i in range(100):
			b = BinObject(min_val + i * bin_size, min_val + (i + 1) * bin_size)
			new_tuple_holder = []
			for t in tuple_holder:
				if t[0] >= b.min and t[0] < b.max:
					b.add_tuple(t)
				else:
					new_tuple_holder.append(t)

			tuple_holder = new_tuple_holder
			bin_holder.append(b)

		#lotta lists
		bin_mins = []			#holds the bin minimums, used for graphing averages
		avgsx = []				#holds the average of the x-direction pressure
		errx = []				#holds the standard devation of the x-direction pressure
		avgsy = []				#  ''
		erry = []				#  ''
		x_no_outliers = []		#holds the data without outliers for x-direction pressure
		y_no_outliers = []		#  ''
		x_areas = []			#holds the areas that have a non-outlier pressure for x-direction
		y_areas = []			#  ''

		#fill all of the lists above
		for b in bin_holder:
			b.do_stats()
			a = b.return_avgs()
			bin_mins.append(b.min)
			avgsx.append(a[0][0])
			errx.append(a[0][1])
			avgsy.append(a[1][0])
			erry.append(a[1][1])

			no_o = b.remove_outliers()
			for i in range(len(no_o[0])):
				x_areas.append(no_o[0][i][0])
				x_no_outliers.append(no_o[0][i][1])
			for j in range(len(no_o[1])):
				y_no_outliers.append(no_o[1][j][1])
				y_areas.append(no_o[1][j][0])


		#generate a graph of the averages/standard deviations 
		plt.figure(1)
		plt.errorbar(bin_mins, avgsx, yerr=errx, fmt="-o")
		plt.errorbar(bin_mins, avgsy, yerr=erry, fmt="-^")
		plt.title("Average and Standard Devation of 100 bins")
		plt.savefig("avgstd")
		plt.close()

		#generate a graph without the outliers, and add in best-fit lines
		z1 = np.polyfit(x_areas, x_no_outliers, 5)
		z2 = np.polyfit(y_areas, y_no_outliers, 5)
		f1 = np.poly1d(z1)
		f2 = np.poly1d(z2)
		x = np.linspace(min(x_areas), max(x_areas), 1000)
		y1 = f1(x)
		y2 = f2(x)

		plt.figure(1)
		plt.plot(x_areas, x_no_outliers, 'ro', y_areas, y_no_outliers, 'bo', x, y1, 'k', x, y2, 'y')

		plt.xlabel("Area / nanoparticle")
		plt.ylabel("Pressure")

		plt.title("Without Outliers")

		plt.savefig("withoutoutliers")
		plt.close()

















