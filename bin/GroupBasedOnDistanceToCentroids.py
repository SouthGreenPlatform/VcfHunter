#!/usr/bin/env python
#
#  Copyright 2017 CIRAD
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, see <http://www.gnu.org/licenses/> or
#  write to the Free Software Foundation, Inc.,
#  51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
#
# -*- coding: utf-8 -*-

import sys
sys.stdout.write('\nLoading modules\n\n')
import optparse
import numpy

def calc_euclid_dist(INDI, CENTROID):
	
	"""
		Calculate euclidean distance between individuals and centroids
		
		:param INDI: An array containing by line dot coordinate.
		:type INDI: array
		:param CENTROID: An array containing by line centroid coordinate.
		:type CENTROID: array
		:return: return an array with distance between individuals (raw) and centroids (column).
		:rtype: array
	"""
	
	list_dist_2_centro = []
	for n in INDI:
		for k in CENTROID:
			list_dist_2_centro.append(numpy.sqrt(numpy.sum((n-k)**2)))
	list_dist_2_centro = numpy.array(list_dist_2_centro).reshape(len(INDI),len(CENTROID))
	
	return list_dist_2_centro

def CreateNumpyArray(FILE, AXES, K, AP):
	
	"""
		Create a numpy array that will be used for euclidean distance calculation
		
		:param FILE: File name of coordinates
		:type FILE: str
		:param AXES: string of axis to work with (each axis should be separated by ":")
		:type AXES: str
		:param K: Group number
		:type K: int
		:param AP: Work on present/absent alleles 'y' or only on present allele 'n'
		:type AP: str
		:return: A tuple with 2 object: a numpy array and correponding marker names
		:rtype: void
	"""
	
	# recording axis that will be used
	axes = list(map(int, AXES.split(':')))
	# if len(axes) < 2:
		# sys.exit('This is embarrassing... The program exited without finishing because at least 2 axes should be filled in --dAxes argument and it is not the case')
	
	#-----------------------------------------------------------
	# Loading coordinates in an array (Absent lines are removed)
	#-----------------------------------------------------------
	my_array = []
	# Opening coordinate file
	file = open(FILE)
	# recording header
	header = file.readline().replace('"','').split()
	# recording coordinates
	col_header = []
	
	if AP == 'n':
		for line in file:
			data = line.split()
			if data:
				if data[0].replace('"','').replace('.',':').split(':')[3] == 'P':
					col_header.append(data[0].replace('"','').replace('.',':'))
					my_list = []
					for ax in axes:
						my_list.append(data[ax])
					my_array.append(list(map(float, my_list)))
		file.close()
	elif AP == 'y':
		for line in file:
			data = line.split()
			if data:
				col_header.append(data[0].replace('"','').replace('.',':'))
				my_list = []
				for ax in axes:
					my_list.append(data[ax])
				my_array.append(list(map(float, my_list)))
		file.close()
	else:
		sys.exit('Unrecognized argument passed to --AP '+AP+' Only possible values are "y" or "n".')
	
	DicoCol = {}
	for i in range(len(col_header)):
		DicoCol[col_header[i]] = i
	
	return (numpy.array(my_array), col_header, DicoCol)


def __main__():
	#Parse Command Line
	parser = optparse.OptionParser(usage="python %prog [options]\n\nProgram designed by Guillaume MARTIN : guillaume.martin@cirad.fr")
	# Wrapper options. 
	parser.add_option( '',	'--centroid',		dest='centroid',		default=None,			help='The centroid coordinate file')
	parser.add_option( '',	'--centcores',		dest='centcores',		default=None,			help='The centroid correspondence with group')
	parser.add_option( '',	'--VarCoord',		dest='VarCoord',		default=None, 			help='The tabulated file of variables coordinates in new axis. (The --prefix + _variables_coordinates.tab file generated when running this script with PCA type)')
	parser.add_option( '',	'--MaxDist',		dest='MaxDist',			default=None,			help='Maximal distance between allele and centroid')
	parser.add_option( '',	'--Axes',			dest='Axes',			default=None,			help='Axes to use. Axis should be separated by ":".')
	parser.add_option( '',	'--groups',			dest='groups',			default=None,			help='Groups to keep. Groups should be separated ":".')
	parser.add_option( '',	'--mat',			dest='mat',				default=None,			help='Allele file in which group will be re-attributed.')
	parser.add_option( '',	'--eval',			dest='eval',			default='n',			help='If yes this program only calculate mean distances between centroids. Possible values: "y", or "n". [Default: %default]')
	parser.add_option( '',	'--prefix',			dest='prefix',			default='Clustered',	help='Prefix for output files. [Default: %default]')
	
	
	
	(options, args) = parser.parse_args()
	
	
	# Loading centroid information
	sys.stdout.write('Recording centroid matrix\n')
	grpToKeep = options.groups.split(':')
	CentroidGroup = {}
	file = open(options.centcores, 'r')
	line = file.readline()
	for line in file:
		data = line.split()
		if data:
			if data[1] in grpToKeep:
				CentroidGroup[data[0]] = data[1]
	file.close()
	
	axes = list(map(int, options.Axes.split(':')))
	col_header = []
	CentroArray = []
	file = open(options.centroid)
	header = ['centro'] + file.readline().split()
	for line in file:
		data = line.split()
		if data:
			if data[0] in CentroidGroup:
				col_header.append(CentroidGroup[data[0]])
				my_list = []
				for ax in axes:
					my_list.append(data[ax])
				CentroArray.append(list(map(float, my_list)))
	file.close()
	
	CentroArray = numpy.array(CentroArray)
	sys.stdout.write('Centroid matrix recorded\n\n')
	
	# Calculating mean euclidean distance between centroids
	sys.stdout.write('Calculating mean distance between two centroids\n')
	Values = calc_euclid_dist(CentroArray, CentroArray)
	ValueShape = Values.shape
	if ValueShape[0] == ValueShape[0]:
		pass
	else:
		sys.exit("There is a bug the centroid distance matrix. The obtained matrix is not square\n")
	Total = []
	for i in range(ValueShape[0]):
		for k in range(ValueShape[0]):
			if i < k:
				Total.append(Values[i,k])
	sys.stdout.write('Mean distance between centroids '+str(sum(Total)/len(Total))+'\n\n')
	
	if options.eval == 'y':
		sys.exit('Program finished without error\n')
	
	MAXDist = float(options.MaxDist)
	
	# Loading allele coordinates
	sys.stdout.write('Recording Matrix\n')
	Matrix = CreateNumpyArray(options.VarCoord, options.Axes, 2, "n")
	sys.stdout.write('Matrix recorded\n\n')
	
	# Calculating euclidean distances to centroids
	sys.stdout.write('Calculating allele distance to centroids\n')
	DistFromCentro = calc_euclid_dist(Matrix[0], CentroArray)
	sys.stdout.write('Allele distance to centroids calculated\n\n')
	
	# Attributing allele to centroid
	sys.stdout.write('Attributing allele to centroid. Alleles that can be attributed to several centroid were attributed to nothing\n')
	MatrixShape = DistFromCentro.shape
	DicoAttribution = {}
	if MatrixShape[1] == len(col_header):
		pass
	else:
		sys.exit("There is a bug. The number of centroid does not match...\n")
	for i in range(MatrixShape[0]):
		alleleGroup = set()
		for k in range(MatrixShape[1]):
			if DistFromCentro[i,k] <= MAXDist:
				alleleGroup.add(col_header[k])
		if len(alleleGroup) == 1:
			Marker = Matrix[1][i]
			if Marker in DicoAttribution:
				sys.exit("There is a bug. A marker is duplicated in the file. This is not allowed ...\n")
			DicoAttribution[Marker] = list(alleleGroup)[0]
		elif len(alleleGroup) > 1:
			print('Allele', Matrix[1][i], 'as more than one possible origin')
	sys.stdout.write('End of allele attribution\n\n')
	
	# Printing files
	sys.stdout.write('Printing file\n')
	outfile = open(options.prefix+'_Clust.tab', 'w')
	file = open(options.mat, 'r')
	line = file.readline()
	header = line.split()
	GroupPos = header.index('K-mean_GROUP') + 1
	outfile.write(line)
	for line in file:
		data = line.split()
		if data:
			Marker = data[0]
			if Marker in DicoAttribution:
				data[GroupPos] = DicoAttribution[Marker]
			else:
				data[GroupPos] = 'g-1'
			outfile.write('\t'.join(data))
			outfile.write('\n')
	outfile.close()
	file.close()
	sys.stdout.write('Program finished without error\n')
	
if __name__ == "__main__": __main__()