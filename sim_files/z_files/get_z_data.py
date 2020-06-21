import numpy as np

file = open("attempt.xyz", "r")

simulation_holder = []							#stores all of the timesteps
time_step_holder = []							#stores an individual timestep

#get rid of the bounds
nop = file.readline()
l = nop.strip().split()
number_of_particles = l[0]
file.readline()

#determines a new timestep
reset = False		

final_z_holder = []		#list containing the data from each histogram
counter = 1

#location of the bins, based on the heights of the particles
#this was done after the analysis to figure out appropriate range for the histograms, just as an fyi
a = np.arange(54.60, 57.00, .01)
#a = np.arange(54.86,56.46,.01)

z_holder = []			#list holding data for a single histogram
for line in file:
	if line.strip() == str(number_of_particles):		#check to see if the next timestep has been reached
		reset = True
		continue

	if reset:					#resets time_step_holder and appends its data to simulation_holder
		#simulation_holder.append(time_step_holder)
		d,e = np.histogram(z_holder,a)
		final_z_holder.append(d)
		z_holder = []
		reset = False
		print("Finished timestep " + str(counter))
		counter += 1

		#we've already read the line with the bounds, so continue
		continue

	#parses xyz data into a tuple and appends it to time_step_holder
	line = line.strip().split()
	line = line[1:]
	line = [float(i) for i in line]
	z_holder.append(line[2])


#write to file
q = open("z_data.txt", "w")
for l in final_z_holder:
	for i in l:
		q.write(str(i) + "\t")
	q.write("\n")
q.close()

#now call graphing_z_stuff.py

