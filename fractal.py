import numpy as np 
from Bio import SeqIO
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt 
import sys
import re

def levy_walk(input_sequence):
	
	origin = np.zeros(3, dtype=int)
	present_start_point = origin 
	points_of_the_levy_walk = []
	sequence_length =  len(input_sequence)

	for i in range(0,sequence_length,1):
		symbol = input_sequence[i]

		if symbol == 'A':
			temp = np.array(present_start_point) 
			temp[0] -=1
			temp[2] = i+1
			points_of_the_levy_walk.append(temp)
			present_start_point = temp
		if symbol == 'T':
			temp = np.array(present_start_point) 
			temp[0] +=1
			temp[2] = i+1
			points_of_the_levy_walk.append(temp)
			present_start_point = temp
		if symbol == 'G':
			temp = np.array(present_start_point) 
			temp[1] -=1
			temp[2] = i+1
			points_of_the_levy_walk.append(temp)
			present_start_point = temp
		if symbol == 'C':
			temp = np.array(present_start_point) 
			temp[1] +=1
			temp[2] = i+1
			points_of_the_levy_walk.append(temp)
			present_start_point = temp

	return points_of_the_levy_walk

n_figure = plt.figure(2)
#n_ax = n_figure.add_subplot(111, projection='3d')
n_ax = n_figure.add_subplot(111)
n_ax.set_xlabel('A-T plane')
n_ax.set_ylabel('G-C plane')
#n_ax.set_zlabel('N')

figure = plt.figure(1)
ax = figure.add_subplot(111)
#ax = figure.add_subplot(111, projection='3d')
ax.set_xlabel('A-T plane')
ax.set_ylabel('G-C plane')
#ax.set_zlabel('N')

s_figure = plt.figure(3)
s_ax = s_figure.add_subplot(111)
#s_ax = s_figure.add_subplot(111, projection='3d')
s_ax.set_xlabel('A-T plane')
s_ax.set_ylabel('G-C plane')
#s_ax.set_zlabel('N')



def plot_points(points,color):
	global figure
	global ax	
	n = len(points)
	x = np.zeros(n,dtype=int)
	y = np.zeros(n,dtype=int)
	z = np.zeros(n,dtype=int)
	for i in range(0 , n, 1):
		point = points[i]
		x[i] = point[0]
		y[i] = point[1]
		z[i] = point[2]
		#print(str(x)+" * "+str(y)+" * "+str(z))
	if color == 'r':
		#ax.plot(x,y,z,label='NON ENHANCER SAMPLE',c=color)
		ax.plot(x,y,color)
	elif color == 'g':
		#n_ax.plot(x,y,z,label='ENHANCER SAMPLE',c=color)
		n_ax.plot(x,y,color)

	s_ax.plot(x,y,color)
	
	

def get_label(file_name, mode ):
	if mode == 0:
		if file_name.startswith('non'):
			return 0
		else:
			return 1

	else:
		if file_name.startswith('non'):
			return 0
		if file_name.startswith('strong'):
			return 1
		if file_name.startswith('weak'):
			return 2

def read_files(input_files):

	for i in range(len(input_files)):
		index = 0
		file_name = input_files[i]
		label = get_label(file_name,0)
		color = 'r'
		if label == 0:
			color = 'r'
		else:
			color = 'g'

		for record in SeqIO.parse(file_name, "fasta"):
			# if index == 100 and label ==1:
			# 	break
			# elif index == 200 and label ==0:
			# 	break
			sequence = re.sub('[0123456789]','',str(record.seq))
			plot_points(levy_walk(sequence.upper()),color)
			index +=1


files = [sys.argv[1],sys.argv[2],sys.argv[3]]
read_files(files)
plt.show()