import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from matplotlib.colors import LogNorm

#parse some files
f = open("area.txt","r")
X = []
for line in f:
	l = line.strip().split()
	X.append(float(l[0]))
f.close()
X = np.array(X)

#set up the histogram range again
Y = np.arange(54.60,57.00,.01)

#scaling code for z position units
Y = np.multiply(Y, 58.6)
Y = np.divide(Y, 10)
Y = np.subtract(Y, 322.5)
Y = np.divide(Y, 7)
#Y = np.multiply(Y, 5.298)

Z = np.loadtxt("z_data.txt",unpack=False)		#reads the txt file as a matrix

#so these are matrix manipulations that were necessary to get what i wanted but i'm not entirely sure why they were necessary
Z2 = np.flip(Z,1)
Z2 = np.flip(Z2,0)
Z2 = np.ndarray.transpose(Z2)
Z2 = np.multiply(Z2,58.6)
#Z2 = np.multiply(Z2,5.298)

#remove the first 39 data points -- cleaning data
X = X[40:]
Z = Z[40:]
Z2 = Z2[40:]


#so this chunk of code is useful if you want to make a histogram at a particular timestep
#i've attached an example one
#'i' is used to locate the particular timestep you want, as well as the 'while' loop
'''
i = 3500

while(i < 4000):
	plt.bar(Y[:-1], Z[i])
	plt.title("nm^2 / nanoparticle : " + str(X[i]))
	plt.savefig("narrow_hist_area_" + str(X[i]) + ".png")
	plt.close()
	i += 50

exit()
'''

#this makes the funky picture with all the colors :)
#i stole it from stackoverflow so i have no clue what most of this does :)))
plt.imshow(Z2, extent=(np.amin(X), np.amax(X), np.amin(Y), np.amax(Y)), cmap=cm.hot, aspect='auto', norm=LogNorm())
#plt.imshow(Z2, extent=(np.amin(X), np.amax(X), np.amin(Y), np.amax(Y)), cmap=cm.hot, aspect='auto', norm=LogNorm())
plt.xlim(np.amax(X), np.amin(X))
plt.ylim(np.amin(Y),np.amax(Y))
plt.xlabel("Area/nanoparticle, nm^2")
plt.ylabel("Z Position (# diameters)")
plt.colorbar()
plt.savefig("height_dist",dpi=500)
