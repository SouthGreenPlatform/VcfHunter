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
sys.stdout.write('loading modules\n')
import os
import csv
import math
import time
import shutil
import random
import datetime
import optparse
import tempfile
import traceback
import fileinput
import subprocess
import multiprocessing as mp
from operator import itemgetter
from inspect import currentframe, getframeinfo

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

import numpy
from sklearn.cluster import KMeans
from sklearn.cluster import MeanShift, estimate_bandwidth
from sklearn.datasets.samples_generator import make_blobs


import matplotlib

if 'VISUALIZE_VAR_2D' in sys.argv:
	matplotlib.use('Agg')
	import matplotlib

from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.backends.backend_pdf import PdfPages

sys.stdout.write('modules loaded\n')

# def draw_color_map():
	
	# import matplotlib.patches as mpatches
	
	# fig = plt.figure(figsize=(11,11))
	# ax = fig.add_subplot(111)
	# list = numpy.arange(1000)
	# ax.set_xlim(0, 1000)
	# ax.set_ylim(0, 1000)
	# for n in list:
		# ax.add_patch(mpatches.Rectangle((n,0), 1, 1000, color=couleur(n/1000)))
	# plt.show()

def couleur(VALUE):
	
	from scipy.stats import norm
	
	v0 = norm(0, 0.25)
	v03 = norm(0.5, 0.25)
	v06 = norm(1, 0.25)
	v1 = norm(1.5, 0.25)
	
	rouge = min(1, v1.pdf(VALUE)/v1.pdf(1.5)+v03.pdf(VALUE)/v03.pdf(0.5))
	verte = min(1, v06.pdf(VALUE)/v06.pdf(1))
	bleu = min(1, v0.pdf(VALUE)/v0.pdf(0))
	
	return (rouge, verte, bleu)

def OuputMultiTread(TABLE, LABELS, DICOALNAME, SetAlname, OUT):
	outfile = open(OUT,'w', 1)
	for data in TABLE:
		if data[0] in SetAlname:
			outfile.write('\t'.join([data[0]] + ['g'+str(LABELS[DICOALNAME[data[0]]])] + data[1:]))
			outfile.write('\n')
	outfile.close()
	return 0

def CreateNumpyArray(FILE, AXES, K, AP):
	
	"""
		Fill a set with chromosomes to exclude
		
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

def ClusteringOutput(CENTROIDSGROUPS, CENTROIDSITERPOS, CENTROID, CORRESPONDANCE, LABELS, MAT, ALNAME, DICOALNAME, PROBA, THREAD, OUT):
	
	"""
		Output results
		
		:param CENTROIDSGROUPS: Centroids groups
		:type CENTROIDSGROUPS: numpy array
		:param CENTROIDSITERPOS: Centroids of centroids obtained with iteration
		:type CENTROIDSITERPOS: numpy array
		:param CENTROID: Coordinates of final centroids
		:type CENTROID: numpy array
		:param CORRESPONDANCE: correspondance between centroids of centroids and final centroids 
		:type CORRESPONDANCE: numpy array
		:param LABELS: Alleles groups
		:type LABELS: numpy array
		:param MAT: Matrix of allele coordinates
		:type MAT: numpy array
		:param ALNAME: A list of allele names
		:type ALNAME: list
		:param DICOALNAME: A dictionnary with allele names and their corresponding position in the matrix
		:type DICOALNAME: dict
		:param PROBA: Probability to be in each groups
		:type PROBA: numpy array
		:param THREAD: Number of processors availables
		:type THREAD: int
		:param OUT: A string corresponding to prefix for output
		:type OUT: str
		:return: perform the k-mean clustering
		:rtype: void
	"""
	
	SetAlname = set(ALNAME)
	
	outfile = open(OUT+'_centroid_coordinates.tab','w')
	outfile.write('\taxis'+'\taxis'.join(list(map(str, range(CENTROID.shape[1]))))+'\n')
	for i in range(CENTROIDSITERPOS.shape[0]):
		outfile.write('centroid'+str(i)+'\t'+'\t'.join(list(map(str, list(CENTROIDSITERPOS[i,:]))))+'\n')
	outfile.close()
	
	outfile = open(OUT+'_centroid_iteration_grouping.tab','w')
	outfile.write('\tK-mean_GROUP\n')
	for i in range(CENTROIDSGROUPS.shape[0]):
		outfile.write('centroid'+str(i)+'\tg'+str(CORRESPONDANCE[CENTROIDSGROUPS[i]])+'\n')
	outfile.close()
	
	outfile = open(OUT+'_group_color.tab','w')
	outfile.write('[color]\n')
	nb_color = CENTROID.shape[0]
	Interval = 1/float(nb_color-1)
	Value = 0
	for n in range(nb_color):
		coul = couleur(Value)
		outfile.write('g'+str(n)+'\t'+':'.join(['red='+str(coul[0]), 'green='+str(coul[1]), 'blue='+str(coul[2]), 'alpha=0.5\n']))
		Value += Interval
	outfile.close()
	
	# tps1 = time.time()
	
	with open(MAT, 'r') as csvfile:
		file = csv.reader(csvfile, delimiter='\t')
		table = [r for r in file]
	header = table[0]
	TotalLines = len(table[1:])
	NbLines = int(TotalLines / THREAD)
	
	i = 0
	list_job = []
	ListToCat = []
	while i+1 < THREAD:
		list_job.append(['OuputMultiTread', table[(i*NbLines)+1:(i+1)*NbLines+1], LABELS, DICOALNAME, SetAlname, OUT+'_'+str(i)+'_kMean_allele.tab'])
		ListToCat.append(OUT+'_'+str(i)+'_kMean_allele.tab')
		i+=1
	list_job.append(['OuputMultiTread', table[(i*NbLines)+1:TotalLines], LABELS, DICOALNAME, SetAlname, OUT+'_'+str(i)+'_kMean_allele.tab'])
	ListToCat.append(OUT+'_'+str(i)+'_kMean_allele.tab')
	
	pool = mp.Pool(processes=THREAD)
	result = pool.map(main, list_job)
	for n in result:
		if n:
			sys.exit('There was a bug during output...')
	pool.close()
	del result
	
	outfile = open(OUT+'_kMean_allele.tab','w', 1)
	outfile.write('\tK-mean_GROUP\t'+'\t'.join(header[1:])+'\n')
	for i in ListToCat:
		shutil.copyfileobj(open(i, 'r'), outfile)
		os.remove(i)
	outfile.close()
	
	# tps2 = time.time()
	# print('temp bloc4:', tps2-tps1)
	
	outfile = open(OUT+'_kMean_gp_prop.tab','w')
	outfile.write('\tg'+'\tg'.join(list(map(str, range(CENTROID.shape[0]))))+'\n')
	for n in range(len(ALNAME)):
		outfile.write('\t'.join(([ALNAME[n]] + list(map(str, list(PROBA[n])))))+'\n')
	outfile.close()
	
	intermediate = list(numpy.unique(LABELS))
	if -1 in intermediate:
		pass
	else:
		group_count = numpy.bincount(LABELS)
		for i in range(group_count.shape[0]):
			sys.stdout.write('Group g'+str(i)+' contained '+str(group_count[i])+' dots\n')

def CalcProbaToBeInGroup(DISTFROMCENTRO):
	
	"""
		Calculate probabilities to be in a group based on distance from centroids.
		
		:param DISTFROMCENTRO: File name of coordinates
		:type DISTFROMCENTRO: str
		:return: Array with probabilities (1: invert euclidean distance, 2: norme so that sum = 1)
		:rtype: array
	"""
	
	## This is not a clean probability calculation. This is just a calculation performed to put statistics and in case you want to filter the clusters!
	## As K-Means assumes spherical cluster, the probability of a datapoint A belonging to cluster B is inversely proportional to the euclidean distance of A to the cluster center of B. 
	## For each point x: 1. compute distance to all cluster centers. 2. invert the distances 3. norming so that the sum equals 1
	
	ListeToReturn = []
	
	for n in DISTFROMCENTRO:
		total = 0
		for k in n:
			if k == 0:
				pass
			else:
				total += 1/float(k)
		for k in n:
			if k == 0:
				ListeToReturn.append(1)
			else:
				ListeToReturn.append((1/float(k))/total)
	
	Dim = DISTFROMCENTRO.shape
	return numpy.array(ListeToReturn).reshape(Dim[0],Dim[1])

def NewMeanShift(FILE, MAT, AXES, K, AP, ITER, THREAD, NewMeanShift, BANDWIDTH, OUT):
	
	"""
		Compute MeanShift clusterization
		
		:param FILE: File name of coordinates
		:type FILE: str
		:param MAT: File name of allele coding
		:type MAT: str
		:param AXES: string of axis to work with (each axis should be separated by ":")
		:type AXES: str
		:param K: Quantile value used for bandwidth estimation
		:type K: float
		:param ITER: Number of iteration of random seed search
		:type ITER: int
		:param THREAD: Number of processors availables
		:type THREAD: int
		:param NewMeanShift: Argument to cluster all points or not
		:type NewMeanShift: str
		:param BANDWIDTH: Bandwidth value for mean shift
		:type BANDWIDTH: Boolean or string
		:param OUT: A string corresponding to prefix for output
		:type OUT: str
		:return: perform the k-mean clustering
		:rtype: void
	"""
	
	sys.stdout.write('Recording Matrix\n')
	Matrix = CreateNumpyArray(FILE, AXES, K, AP)
	
	sys.stdout.write('Performing MeanShift\n')
	if BANDWIDTH:
		bandwidth = float(BANDWIDTH)
		sys.stdout.write('Bandwidth fixed to: '+str(bandwidth)+'\n')
	else:
		bandwidth = estimate_bandwidth(Matrix[0], quantile=K, n_samples=10000, n_jobs=THREAD)
		sys.stdout.write('Bandwidth estimation: '+str(bandwidth)+'\n')
	if NewMeanShift == 'y':
		ms = MeanShift(bandwidth=bandwidth, bin_seeding=True, n_jobs=THREAD).fit(Matrix[0])
	elif NewMeanShift == 'n':
		ms = MeanShift(bandwidth=bandwidth, bin_seeding=True, n_jobs=THREAD, cluster_all=False).fit(Matrix[0])
	else:
		sys.exit('The program exited without finishing because argument "'+NewMeanShift+'" is not recognized in --MeanShiftAll option\n')
	labels = ms.labels_
	cluster_centers = ms.cluster_centers_
	
	
	intermediate = list(numpy.unique(labels))
	if -1 in intermediate:
		intermediate.remove(-1)
	labels_unique = numpy.array(intermediate)
	n_clusters_ = len(labels_unique)
	
	print("number of estimated clusters : %d" % n_clusters_)
	
	# Calculating pseudoprobability to be in a groups
	DistFromCentro = calc_euclid_dist(Matrix[0], cluster_centers)
	Proba = CalcProbaToBeInGroup(DistFromCentro)
	
	Correspondance = numpy.arange(n_clusters_)
	
	sys.stdout.write('Printing files\n')
	ClusteringOutput(Correspondance, cluster_centers, cluster_centers, Correspondance, labels, MAT, Matrix[1], Matrix[2], Proba, THREAD, OUT)

def NewKmean(FILE, MAT, AXES, K, AP, ITER, THREAD, OUT):
	
	"""
		Compute K-mean clusterization
		
		:param FILE: File name of coordinates
		:type FILE: str
		:param MAT: File name of allele coding
		:type MAT: str
		:param AXES: string of axis to work with (each axis should be separated by ":")
		:type AXES: str
		:param K: Group number
		:type K: int
		:param ITER: Number of iteration of random seed search
		:type ITER: int
		:param THREAD: Number of processors availables
		:type THREAD: int
		:param OUT: A string corresponding to prefix for output
		:type OUT: str
		:return: perform the k-mean clustering
		:rtype: void
	"""
	
	Matrix = CreateNumpyArray(FILE, AXES, K, AP)
	
	sys.stdout.write('Running K-mean\n')
	
	# Preparing n independent k-means
	list_kmean = []
	i = 0
	list_job = []
	while i < ITER:
		list_job.append(['newKMean', K, Matrix[0], i])
		i+=1
	
	# Running the N k-mean
	pool = mp.Pool(processes=THREAD)
	result = pool.map(main, list_job)
	CentroidsIterPos = []
	if 'error' in result:
		raise ErrorValue ('Bug in kMean calculation.')
	else:
		for n in result:
			for k in n:
				CentroidsIterPos.append(k)
	pool.close()
	del result
	
	CentroidsIterPos = numpy.array(CentroidsIterPos)
	
	# print(CentroidsIterPos)
	
	# K-mean clustering of centroids and getting centroid positions
	kmeans = KMeans(n_clusters=K, n_init=ITER, tol=1e-10, max_iter=10000, random_state=0, n_jobs=THREAD).fit(CentroidsIterPos)
	Centroids = kmeans.cluster_centers_
	CentroidsGroups = kmeans.labels_
	
	sys.stdout.write('K-mean on centroids\n')
	# It is time to perform the last K-mean with Centroids of Centroids as start points
	kmeans = KMeans(n_clusters=K, init=Centroids, n_init=1, tol=1e-10, max_iter=10000, random_state=0, n_jobs=THREAD).fit(Matrix[0])
	labels = kmeans.labels_
	FinalCentroids = kmeans.cluster_centers_
	
	sys.stdout.write('Calculating distance and pseudoprobabilities\n')
	# Calculating pseudoprobability to be in a groups
	DistFromCentro = calc_euclid_dist(Matrix[0], FinalCentroids)
	Proba = CalcProbaToBeInGroup(DistFromCentro)
	
	# Indentification of which Centroids of centroids correspond to which FinalCentroid
	Correspondance = kmeans.predict(Centroids)
	if len(Correspondance) != K:
		sys.stdout.write('Warning some centroids of centroids are clustered together in the final centroids data! The k-mean may not be appropriate for your data...')
		
	sys.stdout.write('Printing files\n')
	ClusteringOutput(CentroidsGroups, CentroidsIterPos, FinalCentroids, Correspondance, labels, MAT, Matrix[1], Matrix[2], Proba, THREAD, OUT)

def RecordChromToExclude(EXCLCHR):
	
	"""
		Fill a set with chromosomes to exclude
		
		:param EXCLCHR: A string listing chromosomes to exclude
		:type EXCLCHR: str
		:return: A set containing chromosomes to exclude
		:rtype: set
	"""
	
	CHR_TO_EXCLUDE = set()
	
	if not(EXCLCHR == None):
		CHR_TO_EXCLUDE = set(EXCLCHR.split('='))
	return CHR_TO_EXCLUDE

def RecordAccession(VCF, NAMES, MAT):
	
	"""
		Fill a list with accession to work with
		
		:param NAMES: Path to file containing accession names
		:type NAMES: str
		:param VCF: Path to vcf file
		:type VCF: str
		:param MAT: Path to matrix file containing grouping informations
		:type MAT: str
		:return: A list containing accession to work with
		:rtype: list
	"""
	
	ACC_TO_WORK = []
	if NAMES == None:
		sys.stdout.write('No file name was provided in --names argument. All accessions in vcf file (if provided) will be drawn. If no vcf too, all accessions in matrix will be drawn.\n')
		if VCF == None:
			file = open(MAT)
			ACC_TO_WORK = file.readline().replace('K-mean_GROUP','').replace('GROUP','').split()
			file.close()
		else:
			file = open(VCF)
			for line in file:
				data = line.split()
				if data:
					if data[0] == '#CHROM':
						ACC_TO_WORK = data[data.index('FORMAT')+1:]
						break
			file.close()
	else:
		file = open(NAMES)
		for line in file:
			data = line.split()
			if data:
				ACC_TO_WORK.append(data[0])
		file.close()
	
	return ACC_TO_WORK

def GroupToWork(DGROUP, MAT):
	
	"""
		Fill a set with group to work with
		
		:param DGROUP: A string listing group to work with
		:type DGROUP: str
		:param MAT: Path to matrix file containing grouping informations
		:type MAT: str
		:return: A set containing group to work with
		:rtype: set
	"""
	
	GROUP_TO_WORK = set()
	
	if DGROUP == None:
		sys.stdout.write('No file groups were provided in --dGroup argument. All groups defined in matrix will be drawn.\n')
		file = open(MAT)
		header = file.readline().split()
		if 'K-mean_GROUP' in header:
			for line in file:
				data = line.split()
				if data:
					GROUP_TO_WORK.add(data[1:][header.index('K-mean_GROUP')])
		elif 'GROUP' in header:
			for line in file:
				data = line.split()
				if data:
					GROUP_TO_WORK.add(data[1:][header.index('GROUP')])
		else:
			sys.stdout.write('This is embarrassing. The program exited without finishing because no groups were found in the matrix file.\n')
			
	else:
		GROUP_TO_WORK = set(DGROUP.split('='))
	
	return GROUP_TO_WORK
	
def RecordAlleleGrouping(MAT, GROUP_TO_WORK, CHR_TO_EXCLUDE, DICO_ALLELE, DICO_SITES):
	
	"""
		Fill a set with chromosomes to exclude
		
		:param MAT: Path to matrix file containing grouping informations
		:type MAT: str
		:param GROUP_TO_WORK: A set containing groups to work with
		:type GROUP_TO_WORK: set
		:param DICO_ALLELE: An empty dictionnary
		:type DICO_ALLELE: dictionnary
		:param DICO_SITES: An empty set
		:type DICO_SITES: set
		:return: Fill the set and dictionnary passed in DICO_SITES and DICO_ALLELE arguments
		:rtype: void
	"""
	
	sys.stdout.write("Recording alleles in grouping informations...\n")
	file = open(MAT)
	header = file.readline().split()
	if 'K-mean_GROUP' in header:
		for line in file:
			data = line.split()
			if data:
				variant = data[0].split(':')
				if not(variant[0] in CHR_TO_EXCLUDE):
					if not(variant[3] == 'A') and data[1:][header.index('K-mean_GROUP')] in GROUP_TO_WORK:
						DICO_ALLELE[data[0]] = data[1:][header.index('K-mean_GROUP')]
						DICO_SITES.add(':'.join(variant[0:2]))
	elif 'GROUP' in header:
		for line in file:
			data = line.split()
			if data:
				variant = data[0].split(':')
				if not(variant[0] in CHR_TO_EXCLUDE):
					if not(variant[3] == 'A') and data[1:][header.index('GROUP')] in GROUP_TO_WORK:
						DICO_ALLELE[data[0]] = data[1:][header.index('GROUP')]
						DICO_SITES.add(':'.join(variant[0:2]))
	else:
		sys.stdout.write('This is embarrassing. The program exited without finishing because no groups were found in the matrix file.\n')
	file.close()
	
	sys.stdout.write("Done\n")

def extract_function_name():
	"""Extracts failing function name from Traceback

	by Alex Martelli
	http://stackoverflow.com/questions/2380073/\
	how-to-identify-what-function-call-raise-an-exception-in-python
	"""
	tb = sys.exc_info()[-1]
	stk = traceback.extract_tb(tb, 1)
	fname = stk[0][3]
	return fname

def main(JOB):
	
	try:
		if JOB[0] == 'kMean':
			result_main = kMean(JOB[1], JOB[2], JOB[3])
			sys.stdout.write('Iteration : '+str(JOB[4])+' done\n')
			sys.stdout.flush()
		elif JOB[0] == 'newKMean':
			RESULT_MAIN = KMeans(n_clusters=JOB[1], init='random', n_init=1, tol=1e-10, max_iter=10000, random_state=0, n_jobs=1).fit(JOB[2])
			result_main = RESULT_MAIN.cluster_centers_
			sys.stdout.flush()
		elif JOB[0] == 'OuputMultiTread':
			result_main = OuputMultiTread(JOB[1], JOB[2], JOB[3], JOB[4], JOB[5])
			sys.stdout.flush()
		else:
			raise ErrorValue ('Unknown function name: '+str(JOB[0]))
	except Exception as e:
		print(e)
		result_main = 'error'
	finally:
		return result_main

def get_genotype(VCF, NAME, PREFIX):
	
	"""
		Return a file containing the genotype per accession (column) and per site (raw)
		
		:param VCF: A vcf file.
		:type VCF: vcf
		:param NAME: genotype names.
		:type NAME: string
		:param PREFIX: Prefix of the output file.
		:type PREFIX: str
		:return: A tab file containing genotype for each position.
		:rtype: void
	"""
	
	# creating output file
	outfile = open(PREFIX+'_genotype.gen','w')
	
	# creating list of accession to output
	names = []
	file = open(NAME)
	mot = ['#CHROM', 'POS']
	for line in file:
		data = line.split()
		if len(data) > 1:
			sys.exit('There is a problem in '+NAME+' file format. This file must contain per line only one accession name')
		if data:
			names.append(data[0])
			mot.append(data[0])
	outfile.write('\t'.join(mot)+'\n')
	
	# Reading and recording informations
	file = open(VCF)
	for line in file:
		data = line.split()
		if data:
			# recording header
			if data[0] == '#CHROM':
				header = list(data)
			# not working on headers of the file
			elif data[0][0] != '#':
				mot = []
				for acc in names:
					genotype = recup_geno(data, header, acc)
					if mot:
						mot.append('/'.join(genotype[2:]))
					else:
						mot.append(genotype[0])
						mot.append(genotype[1])
						mot.append('/'.join(genotype[2:]))
				outfile.write('\t'.join(mot)+'\n')
	outfile.close()

def get_genotype_and_group(VCF, MAT, DGROUP, NAME, PREFIX):
	
	"""
		Return a file containing the genotype per accession (column) and per site (raw) and corresponding group
		
		:param VCF: A vcf file.
		:type VCF: vcf
		:param MAT: Path to matrix file containing grouping informations
		:type MAT: str
		:param DGROUP: A string listing group to work with
		:type DGROUP: str
		:param NAME: genotype names.
		:type NAME: string
		:param PREFIX: Prefix of the output file.
		:type PREFIX: str
		:return: A tab file containing genotype for each position.
		:rtype: void
	"""
	
	# Recording group to work with
	GROUP_TO_WORK = set(DGROUP.split(':'))
	
	# Recording for each position and each allele, the grouping
	file = open(MAT)
	header = file.readline().split()
	dico_pos = {}
	if 'K-mean_GROUP' in header:
		column = header.index('K-mean_GROUP')+1
		for line in file:
			data = line.split()
			group = data[column]
			if group in GROUP_TO_WORK:
				SNP_ID = data[0].split(':')
				if SNP_ID[3] == 'P':
					chr = SNP_ID[0]
					pos = SNP_ID[1]
					allele = SNP_ID[2]
					if not (chr in dico_pos):
						dico_pos[chr] = {}
					if not (pos in dico_pos[chr]):
						dico_pos[chr][pos] = {}
					if allele in dico_pos[chr][pos]:
						sys.exit('There is a bug in the matrix file: The same allele is found more than once???')
					dico_pos[chr][pos][allele] = group
	
	elif 'GROUP' in header:
		column = header.index('GROUP')+1
		for line in file:
			data = line.split()
			group = data[column]
			if group in GROUP_TO_WORK:
				SNP_ID = data[0].split(':')
				if SNP_ID[3] == 'P':
					chr = SNP_ID[0]
					pos = SNP_ID[1]
					allele = SNP_ID[2]
					if not (chr in dico_pos):
						dico_pos[chr] = {}
					if not (pos in dico_pos[chr]):
						dico_pos[chr][pos] = {}
					if allele in dico_pos[chr][pos]:
						sys.exit('There is a bug in the matrix file: The same allele is found more than once???')
					dico_pos[chr][pos][allele] = group
	
	else:
		sys.exit('No groups found in the matrix file')
	file.close()
	
	# creating output file
	outfile = open(PREFIX+'_genotype.gen','w')
	
	# creating list of accession to output
	names = []
	file = open(NAME)
	mot = ['#CHROM', 'POS']
	for line in file:
		data = line.split()
		if len(data) > 1:
			sys.exit('There is a problem in '+NAME+' file format. This file must contain per line only one accession name')
		names.append(data[0])
		mot.append('Gen:'+data[0])
		mot.append('Group:'+data[0])
	outfile.write('\t'.join(mot)+'\n')
	
	# Reading and recording informations
	file = open(VCF)
	for line in file:
		data = line.split()
		if data:
			# recording header
			if data[0] == '#CHROM':
				header = list(data)
			# not working on headers of the file
			elif data[0][0] != '#':
				mot = []
				for acc in names:
					genotype = recup_geno(data, header, acc)
					grouping = []
					if genotype[0] in dico_pos:
						if genotype[1] in dico_pos[genotype[0]]:
							for n in genotype[2:]:
								if n in dico_pos[genotype[0]][genotype[1]]:
									grouping.append(dico_pos[genotype[0]][genotype[1]][n])
								else:
									grouping.append('.')
						else:
							for n in genotype[2:]:
								grouping.append('.')
					else:
						for n in genotype[2:]:
							grouping.append('.')	
					if mot:
						mot.append('/'.join(genotype[2:]))
						mot.append('/'.join(grouping))
					else:
						mot.append(genotype[0])
						mot.append(genotype[1])
						mot.append('/'.join(genotype[2:]))
						mot.append('/'.join(grouping))
				if genotype[0] in dico_pos:
					if genotype[1] in dico_pos[genotype[0]]:
						# print (mot)
						outfile.write('\t'.join(mot)+'\n')
	outfile.close()

def random_sub_set(VCF, NRAND, PREFIX):
	
	"""
		Return a matrix file containing accession identity
		
		:param VCF: A vcf file.
		:type VCF: vcf
		:param NRAND: Number of lines to get.
		:type NRAND: int
		:param PREFIX: Prefix of the output file.
		:type PREFIX: str
		:return: A vcf file subset of the original vcf.
		:rtype: void
	"""
	
	# Recording variant line in a file
	file = open(VCF)
	i = 0
	list_line = []
	for line in file:
		data = line.split()
		if data:
			i += 1
			if data[0][0] != "#":
				list_line.append(i)
	file.close()
	
	# Getting lines to keep
	dico_keep = set()
	i = 0
	while i < NRAND:
		i += 1
		valeur = random.choice(list_line)
		dico_keep.add(valeur)
		list_line.remove(valeur)
	
	# Printing vcf subset
	outfile = open(PREFIX+'_subset.vcf', 'w')
	file = open(VCF)
	i = 0
	for line in file:
		data = line.split()
		if data:
			i += 1
			if data[0][0] == "#":
				outfile.write(line)
			elif i in dico_keep:
				outfile.write(line)
	outfile.close()

def recup_geno(LINE, HEADER, ACCESSION):
	
	#doc string
	"""
		Return the genotype of an accession from a vcf line
		
		:param LINE: A list corresponding to a vcf line
		:type LINE: list
		:param HEADER: A list containing the header of the vcf
		:type HEADER: list
		:param ACCESSION: ACCESSION name to return the genotype
		:type ACCESSION: str
		:return: A list with [0] --> chromosome, [1] --> position (integer), [2] --> allele 1, [3] --> allele 2, ... [n] --> allele [n-1]
		:rtype: list
	"""
	
	# Getting possible alleles
	liste_allele = []
	# Getting reference allele
	liste_allele.append(LINE[HEADER.index('REF')])
	# Getting alternative alleles
	alt_all = LINE[HEADER.index('ALT')].split(',')
	for n in alt_all:
		liste_allele.append(n)
	
	# Getting accession format
	format_col = LINE[HEADER.index('FORMAT')].split(':')
	
	# Getting genotypes
	geno = LINE[HEADER.index(ACCESSION)].split(':')[format_col.index('GT')].replace('/', '|').split('|')
	
	# Preparing the list to return
	list2return = [LINE[HEADER.index('#CHROM')], LINE[HEADER.index('POS')]]
	for n in geno:
		if n == '.':
			list2return.append('NA')
		else:
			list2return.append(liste_allele[int(n)])
	
	# returning list
	return list2return

def get_chr_size(LINE):
	
	"""
		Get reference chromosome size from vcf file
		
		:param LINE: The ##contig line of the vcf file splitted on spaces
		:type LINE: list
		:return: A list with [0] --> chromosome name, [1] --> chromosome size
		:rtype: list
	"""
	
	split_on_chev = LINE[0].split('>')[0]
	
	split_on_eq = split_on_chev.split('=')
	
	return [split_on_eq[2].split(',')[0], int(split_on_eq[3])]

def get_chr_size_from_fasta(FASTA):
	
	"""
		Get reference chromosome size from fasta file
		
		:param FASTA: A fasta file containing reference sequence
		:type FASTA: fasta
		:return: Dictionary with key = chromosome name, value = chromosome size
		:rtype: dictionary
	"""
	
	# Loading fasta
	sequence_dict = SeqIO.index(FASTA, "fasta")
	# Recording chromosome size
	dic_size = {}
	for n in sequence_dict:
		dic_size[n] = len(str(sequence_dict[n].seq))
	
	return dic_size

def ident_autapo(DICO):
	
	"""
		Identification accession bearing autapomorphic variants (if at least 2 individuals have been genotyped)
		
		:param DICO: A dictionary containing genotype for each accession at a given site
		:type DICO: dictionary
		:return: A dictionary (set()) containing accessions with autapomorphic variant 
		:rtype: dictionary
	"""
	
	dico_to_return = set()
	for n in DICO:
		for k in DICO[n]:
			not_found = 1
			not_all_empty = 0
			for z in DICO:
				if z != n:
					if k in DICO[z]:
						not_found = 0
					if len(DICO[z]) > 0:
						not_all_empty = 1
			if not_found and not_all_empty:
				dico_to_return.add(n)
	return dico_to_return

def stat_on_vcf(VCF, PREFIX, NAMES, GFF3):
	
	"""
		Get statistics various statistics from vcf file
		
		:param VCF: A vcf file containing variant calling
		:type VCF: vcf
		:param PREFIX: Prefix of the output file
		:type PREFIX: str
		:param NAMES: Path to file containing accession to treat
		:type NAMES: str
		:param GFF3: A gff3 file.
		:type GFF3: gff3
		:return: 2 text files containing statistics, one file with general statistics and one file for per accession statistics
		:rtype: void
	"""
	
	dico_gene_gff3 = {}
	dico_gff3_position = {}
	if not(GFF3 == None):
		Record_CDS_and_UTR(GFF3, dico_gene_gff3)
		get_CDS_and_UTR_position(GFF3, dico_gff3_position)
	
	# Recording accession to treate
	file = open(NAMES)
	DICO_NAME = set()
	for line in file:
		data = line.split()
		if data:
			DICO_NAME.add(data[0])
	file.close()
	
	dico_general_stat = {}
	dico_general_stat['status'] = {}
	dico_general_stat['nb_var'] = {}
	dico_general_stat['type'] = {}
	dico_general_stat['transition_transversion'] = {}
	
	dico_general_stat['type']['SNP'] = 0
	dico_general_stat['type']['indel'] = 0
	
	dico_general_stat['transition_transversion']['A vs T'] = 0
	dico_general_stat['transition_transversion']['A vs G'] = 0
	dico_general_stat['transition_transversion']['A vs C'] = 0
	dico_general_stat['transition_transversion']['T vs C'] = 0
	dico_general_stat['transition_transversion']['T vs G'] = 0
	dico_general_stat['transition_transversion']['C vs G'] = 0
	
	
	dico_general_stat['UTR'] = 0
	dico_general_stat['CDS'] = 0
	dico_general_stat['intron'] = 0
	dico_general_stat['unannotated'] = 0
	dico_general_stat['GeneWithVar'] = set()
	
	dico_acc = {}
	for n in DICO_NAME:
		dico_acc[n] = {}
		dico_acc[n]['missing'] = 0
		dico_acc[n]['nb_homo_ref'] = 0
		dico_acc[n]['nb_homo_alt'] = 0
		dico_acc[n]['hetero'] = {}
		dico_acc[n]['autapo'] = 0
	
	# Reading and recording informations
	file = open(VCF)
	for line in file:
		data = line.split()
		if data:
			# recording header
			if data[0] == '#CHROM':
				header = list(data)
			# not working on headers of the file
			elif data[0][0] != '#':
				# calculating variant position (CDS, UTR, intron, unannotated)
				if data[0] in dico_gene_gff3:
					if int(data[1]) in dico_gene_gff3[data[0]]['CDS']:
						dico_general_stat['CDS'] += 1
					elif int(data[1]) in dico_gene_gff3[data[0]]['UTR']:
						dico_general_stat['UTR'] += 1
					elif int(data[1]) in dico_gene_gff3[data[0]]['intron']:
						dico_general_stat['intron'] += 1
					else:
						dico_general_stat['unannotated'] += 1
				else:
					dico_general_stat['unannotated'] += 1
				# calculating gene with variant number
				if data[0] in dico_gff3_position:
					if int(data[1]) in dico_gff3_position[data[0]]:
						dico_general_stat['GeneWithVar'].update(dico_gff3_position[data[0]][int(data[1])])
				# initialisation of a dictionnary to identify autapomorphy
				dico_autapo = {}
				# getting reference allele
				ref_allele = data[header.index('REF')]
				
				# getting alternative allele
				alt_allele = data[header.index('ALT')].split(',')
				
				# getting variant status
				status = data[header.index('FILTER')]
				
				# recording variant status number
				if not (status in dico_general_stat['status']):
					dico_general_stat['status'][status] = 1
				else:
					dico_general_stat['status'][status] += 1
				
				# identifying SNP and indel
				dic_snp = set(alt_allele)
				dic_snp.add(ref_allele)
				if IsIndel(dic_snp):
					dico_general_stat['type']['indel'] += 1
				else:
					dico_general_stat['type']['SNP'] += 1
				
				# calculating transition - transversion (for SNP sites only)
				if IsSNP(dic_snp):
					for alt in alt_allele:
						if alt == 'A' and ref_allele == 'T':
							dico_general_stat['transition_transversion']['A vs T'] += 1
						elif ref_allele == 'A' and alt == 'T':
							dico_general_stat['transition_transversion']['A vs T'] += 1
						elif alt == 'C' and ref_allele == 'G':
							dico_general_stat['transition_transversion']['C vs G'] += 1
						elif ref_allele == 'C' and alt == 'G':
							dico_general_stat['transition_transversion']['C vs G'] += 1
						elif alt == 'A' and ref_allele == 'C':
							dico_general_stat['transition_transversion']['A vs C'] += 1
						elif ref_allele == 'A' and alt == 'C':
							dico_general_stat['transition_transversion']['A vs C'] += 1
						elif alt == 'A' and ref_allele == 'G':
							dico_general_stat['transition_transversion']['A vs G'] += 1
						elif ref_allele == 'A' and alt == 'G':
							dico_general_stat['transition_transversion']['A vs G'] += 1
						elif alt == 'C' and ref_allele == 'T':
							dico_general_stat['transition_transversion']['T vs C'] += 1
						elif ref_allele == 'C' and alt == 'T':
							dico_general_stat['transition_transversion']['T vs C'] += 1
						elif alt == 'T' and ref_allele == 'G':
							dico_general_stat['transition_transversion']['T vs G'] += 1
						elif ref_allele == 'T' and alt == 'G':
							dico_general_stat['transition_transversion']['T vs G'] += 1
						else:
							sys.stdout.write('bug in transition transversion identification '+alt+' '+ref_allele+'\n'+line)
				
				###### working by accessions ######
				total_alleles = set()
				for n in DICO_NAME:
					genotype = recup_geno(data, header, n)
					genotype_set = set(genotype[2:])
					total_alleles.update(genotype_set)
					
					# recording genotype for autapomorphy detection
					list4autapo = set(genotype[2:])
					list4autapo.discard('NA')
					dico_autapo[n] = list4autapo
					
					# calculating missing number
					if 'NA' in genotype_set:
						dico_acc[n]['missing'] += 1
					# There is no missing data
					else:
						# calculating heterozygous site number
						if len(genotype_set) > 1:
							if len(genotype_set) in dico_acc[n]['hetero']:
								dico_acc[n]['hetero'][len(genotype_set)] += 1
							else:
								dico_acc[n]['hetero'][len(genotype_set)] = 1
						# calculating homozygous sites identical to reference
						elif ref_allele in genotype_set:
							dico_acc[n]['nb_homo_ref'] += 1
						# calculating homozygous sites different to reference
						else:
							dico_acc[n]['nb_homo_alt'] += 1
				
				# calculating variant number
				total_alleles.discard("NA")
				if not(len(total_alleles) in dico_general_stat['nb_var']):
					dico_general_stat['nb_var'][len(total_alleles)] = 1
				else:
					dico_general_stat['nb_var'][len(total_alleles)] += 1
				# if len(total_alleles) == 3:
					# print (line)
				# print (total_alleles)
				# calculating autapomorphy
				autapo = ident_autapo(dico_autapo)
				for k in autapo:
					dico_acc[k]['autapo'] += 1
					# print(k, data[0], data[1], 1, 'autapo')
	
	outfile = open(PREFIX+'_general.stat','w')
	
	outfile.write('*****General statistics*****\n')
	outfile.write('  **Variant status**\n')
	for n in dico_general_stat['status']:
		outfile.write('\t'+n+' number: '+str(dico_general_stat['status'][n])+'\n')
	outfile.write('  **Variant type**\n')
	for n in dico_general_stat['type']:
		outfile.write('\t'+n+' number: '+str(dico_general_stat['type'][n])+'\n')
	outfile.write('  **Alternative number**\n')
	for n in dico_general_stat['nb_var']:
		outfile.write('\tNumber of sites with '+str(n)+' allele(s) (in treated accessions): '+str(dico_general_stat['nb_var'][n])+'\n')
	outfile.write('  **transition-transversion**\n')
	for n in dico_general_stat['transition_transversion']:
		outfile.write('\tNumber of sites with '+n+' mutation: '+str(dico_general_stat['transition_transversion'][n])+'\n')
	outfile.write('  **Variant in annotation**\n')
	outfile.write('\tnumber of variant sites in CDSs: '+str(dico_general_stat['CDS'])+'\n')
	outfile.write('\tnumber of variant sites in UTRs: '+str(dico_general_stat['UTR'])+'\n')
	outfile.write('\tnumber of variant sites in introns: '+str(dico_general_stat['intron'])+'\n')
	outfile.write('\tnumber of variant sites in unannotated regions: '+str(dico_general_stat['unannotated'])+'\n')
	outfile.write('\tnumber of gene with variant sites: '+str(len(dico_general_stat['GeneWithVar']))+'\n')
	outfile.close()
	
	outfile = open(PREFIX+'_accession.stat','w')
	outfile.write('Statistics\Accession')
	
	accession_order = sorted(list(dico_acc.keys()))
	
	for n in accession_order:
		outfile.write('\t'+n)
	outfile.write('\n')
	
	outfile.write('Autapomorphic alleles number:')
	for n in accession_order:
		outfile.write('\t'+str(dico_acc[n]['autapo']))
	outfile.write('\n')
	
	outfile.write('Missing_call_number:')
	for n in accession_order:
		outfile.write('\t'+str(dico_acc[n]['missing']))
	outfile.write('\n')
	
	outfile.write('Number_of_homozygous_sites_identical_to_reference:')
	for n in accession_order:
		outfile.write('\t'+str(dico_acc[n]['nb_homo_ref']))
	outfile.write('\n')
	
	outfile.write('Number_of_homozygous_sites_different_to_reference:')
	for n in accession_order:
		outfile.write('\t'+str(dico_acc[n]['nb_homo_alt']))
	outfile.write('\n')
	
	dico_hetero_name = set()
	for n in dico_acc:
		for j in dico_acc[n]['hetero']:
			dico_hetero_name.add(j)
	
	for n in dico_hetero_name:
		outfile.write('Number_of_heterozygous_sites_with_'+str(n)+'_alleles:')
		for k in accession_order:
			if n in dico_acc[k]['hetero']:
				outfile.write('\t'+str(dico_acc[k]['hetero'][n]))
			else:
				outfile.write('\t0')
		outfile.write('\n')	
	outfile.close()
	
def filter_vcf(VCF, NAMES, OUTGROUP, PREFIX, RMTYPE, MINCOV, MINAL, NMISS, RMALALT):
	
	"""
		Filter vcf file and output a filtered vcf
		
		:param VCF: A vcf file containing variant calling
		:type VCF: vcf
		:param NAMES: Path to file containing accession to treat
		:type NAMES: str
		:param OUTGROUP: Path to file containing accession not to consider for filtering but to print in the output
		:type OUTGROUP: str
		:param PREFIX: Prefix of the output file
		:type PREFIX: str
		:param RMTYPE: Type of variant to remove
		:type RMTYPE: str
		:param MINCOV: Minimal coverage to keep a genotype
		:type MINCOV: int
		:param MINAL: Minimal allele coverage to keep a genotype
		:type MINAL: int
		:param NMISS: Maximal missing genotype to keep a variant site
		:type NMISS: int
		:param RMALALT: Number of alternative variant to remove the variant site.
		:type RMALALT: str
		:return: A filtered vcf file
		:rtype: void
	"""
	nb_remove = 0
	nb_kept = 0
	nb_autapo_remove = 0
	nb_allele_remove = 0
	nb_missing_remove = 0
	nb_tag_remove = 0
	nb_INDEL_remove = 0
	nb_SNP_remove = 0
	
	# Recording accession to treate
	file = open(NAMES)
	DICO_NAME = set()
	for line in file:
		data = line.split()
		if data:
			DICO_NAME.add(data[0])
	file.close()
	
	
	# Recording outgroup
	if OUTGROUP == None:
		DICO_OUTGROUP = set()
	else:
		file = open(OUTGROUP)
		DICO_OUTGROUP = set()
		for line in file:
			data = line.split()
			if data:
				DICO_OUTGROUP.add(data[0])
		file.close()
	
	# recording variant status to exclude
	if RMTYPE == None:
		exclude = []
	else:
		exclude = RMTYPE.split(':')
	
	# recording variant alternatives to exclude
	if  RMALALT == None:
		exclude_al = []
	else:
		exclude_al = [int(i) for i in RMALALT.split(':')]
	
	# Creating output file
	outfile = open(PREFIX+'_filt.vcf','w')
	
	# Reading vcf file
	PrintFilter = 1
	file = open(VCF)
	for line in file:
		data = line.split()
		if data:
			if data[0][0:8] == "##contig" and PrintFilter:
				outfile.write('##FILTER=<ID=TAGs, Description="'+' '.join(exclude)+' variant are removed">\n')
				outfile.write('##FILTER=<ID=COVERAGE, Description="Genotype having less than '+str(MINCOV)+' x coverage and less than '+str(MINAL)+' x coverage for each allele are converted to missing">\n')
				outfile.write('##FILTER=<ID=MISSING, Description="SNP with more than '+str(NMISS)+' genotype missing are removed">\n')
				PrintFilter = 0
			# Recording header
			if data[0] == '#CHROM':
				header = list(data)
				list_to_print_in_output = []
				for n in header:
					list_to_print_in_output.append(n)
					if n == 'FORMAT':
						break
				for n in header:
					if (n in DICO_OUTGROUP) or (n in DICO_NAME):
						list_to_print_in_output.append(n)
				outfile.write('\t'.join(list_to_print_in_output)+'\n')
			# Printing headers
			elif data[0][0] == '#':
				outfile.write(line)
			# Working on variant
			else:
				nb_alt = set()
				remove = 0
				# Identification we should remove the variant based on the FILTER tag
				filter_tag = data[header.index('FILTER')].split(';')
				for n in filter_tag:
					if n in exclude:
						remove = 1
				if remove == 1:
					nb_tag_remove += 1
				# Initialisation of the corrected line
				list_to_print_in_output = []
				for n in header:
					list_to_print_in_output.append(data[header.index(n)])
					if n == 'FORMAT':
						break
				# Identification genotype to convert to missing du to MINCOV and MINAL parameters
				nb_missing = 0
				ref_allele = data[header.index('REF')]
				alt_allele = data[header.index('ALT')].split(',')
				flag_format = data[header.index('FORMAT')].split(':')
				dico_autapo = {}
				for accession in header:
					# working on accession to treat
					if accession in DICO_NAME:
						filtered_accession = filter_accession(accession, data, header, flag_format, MINCOV, MINAL)
						nb_alt.update(filtered_accession[2])
						list_to_print_in_output.append(filtered_accession[0])
						if filtered_accession[1]:
							nb_missing += 1
						for z in filtered_accession[2]:
							if z in dico_autapo:
								dico_autapo[z] += 1
							else:
								dico_autapo[z] = 1
					# working on outgroup (not taken in acount for variant filtering)
					elif accession in DICO_OUTGROUP:
						list_to_print_in_output.append(filter_accession(accession, data, header, flag_format, MINCOV, MINAL)[0])
				# Validating there is not to much missing data
				if nb_missing > NMISS:
					remove = 1
					nb_missing_remove += 1
				if 'AUTAPO' in exclude:
					# Validating if it is an autapomorphy
					if cherche_autapo(dico_autapo):
						remove = 1
						nb_autapo_remove += 1
				if 'SNP' in exclude:
					# Validating if it is a SNP alone
					dic_snp = set(alt_allele)
					dic_snp.add(ref_allele)
					if IsSNP(dic_snp):
						remove = 1
						nb_SNP_remove += 1
				if 'INDELS' in exclude:
					# Validating if it is a SNP alone
					dic_snp = set(alt_allele)
					dic_snp.add(ref_allele)
					if IsIndel(dic_snp):
						remove = 1
						nb_INDEL_remove += 1
				nb_alt.discard('.')
				if len(nb_alt) in exclude_al:
					remove = 1
					nb_allele_remove += 1
				if remove:
					nb_remove += 1
				else:
					nb_kept += 1
					outfile.write('\t'.join(list_to_print_in_output)+'\n')
	sys.stdout.write('Removed variant: '+str(nb_remove)+'\n')
	sys.stdout.write('\tRemoved variant (missing): '+str(nb_missing_remove)+'\n')
	sys.stdout.write('\tRemoved variant (tag): '+str(nb_tag_remove)+'\n')
	sys.stdout.write('\tRemoved variant (autapomorphy): '+str(nb_autapo_remove)+'\n')
	sys.stdout.write('\tRemoved variant (SNP): '+str(nb_SNP_remove)+'\n')
	sys.stdout.write('\tRemoved variant (INDEL): '+str(nb_INDEL_remove)+'\n')
	sys.stdout.write('\tRemoved variant (bad allele number): '+str(nb_allele_remove)+'\n')
	sys.stdout.write('Kept variant: '+str(nb_kept)+'\n')

def filter_accession(ACCESSION, DATA, HEADER, FLAG_FORMAT, MINCOV, MINAL):
	
	"""
		Filter accession based on MINCOV and MINAL parameters
		
		:param ACCESSION: Accession to treat.
		:type ACCESSION: str
		:param DATA: A list corresponding to a vcf line.
		:type DATA: list
		:param HEADER: A list corresponding to the vcf header line.
		:type HEADER: list
		:param MINCOV: Minimal coverage to keep a genotype.
		:type MINCOV: int
		:param MINAL: Minimal allele coverage to keep a genotype.
		:type MINAL: int
		:return: A list with [0] --> gentype calling recalculated, [1] --> missing genotype (bolean), [2] --> set listing alleles found in the genotype
		:rtype: list
	"""
	
	missing = False
	to_print = ''
	accession_var = DATA[HEADER.index(ACCESSION)].split(':')
	# Looking if a genotype has been called
	geno = accession_var[FLAG_FORMAT.index('GT')].replace('/','|').split('|')
	if '.' in geno:
		missing = True
		to_print = DATA[HEADER.index(ACCESSION)]
		forautapo = set()
	else:
		# Looking if site coverage is sufficient
		convert_to_missing = False
		if not('DP' in FLAG_FORMAT):
			dp_cov = -1
		elif FLAG_FORMAT.index('DP') >= len(accession_var):
			dp_cov = -1
		elif accession_var[FLAG_FORMAT.index('DP')] == '.':
			dp_cov = -1
		else:
			dp_cov = int(accession_var[FLAG_FORMAT.index('DP')])
		# Not enough coverage
		if dp_cov < MINCOV:
			convert_to_missing = True
		# Enough coverage
		else:
			ad_cov = accession_var[FLAG_FORMAT.index('AD')].split(',')
			for n in geno:
				if int(ad_cov[int(n)]) < MINAL:
					convert_to_missing = True
		# The genotype should be converted to missing
		if convert_to_missing:
			missing = True
			missing_geno = []
			for n in geno:
				missing_geno.append('.')
			accession_var[FLAG_FORMAT.index('GT')] = '/'.join(missing_geno)
			to_print = ':'.join(accession_var)
			forautapo = set()
		else:
			to_print = DATA[HEADER.index(ACCESSION)]
			# Recording alleles found for autapomorphy search
			forautapo = set(geno)
	
	return [to_print, missing, forautapo]

def IsSNP(DICO):
	
	"""
		Identify if it is a SNP alone
		
		:param DICO: A dictionary of alleles.
		:type DICO: dictionary
		:return: A bolean with true if it is a SNP variant site only
		:rtype: bolean
	"""
	
	max_size = 0
	deleted_allele = False
	for n in DICO:
		max_size = max(len(n),max_size)
		if n == "*":
			deleted_allele = True
	
	if max_size == 1 and deleted_allele == False:
		return True
	else:
		return False

def IsIndel(DICO):
	
	"""
		Identify if it is an indel
		
		:param DICO: A dictionary of alleles.
		:type DICO: dictionary
		:return: A bolean with true if it is a SNP variant site only
		:rtype: bolean
	"""
	
	max_size = 0
	for n in DICO:
		max_size = max(len(n),max_size)
	
	if max_size == 1:
		return False
	else:
		return True

def cherche_autapo(DICO):
	
	"""
		Identify autapomorphies in the set
		
		:param DICO: A dictionary counting the number of accession bearing a given allele.
		:type DICO: dictionary
		:return: A bolean with true if the variant is an autapomorphy and false if not
		:rtype: bolean
	"""
	
	trouve = 0
	for n in DICO:
		if DICO[n] > 1:
			trouve += 1
	if trouve > 1:
		return 0
	else:
		return 1

def compare2calling(VCF, VCF2, COMP1, COMP2):
	
	"""
		Compare two variant calling
		
		:param VCF: A vcf file.
		:type VCF: vcf
		:param VCF2: A vcf file.
		:type VCF2: vcf
		:param COMP1: First accession name to compare.
		:type COMP1: str
		:param COMP2: Second accession name to compare.
		:type COMP2: str
		:return: Print to stdout number of call sites specific to comp1 and specific to comp2, number of shared sites, number of identical calls.
		:rtype: void
	"""
	
	if VCF2 == None and COMP2 == None:
		sys.exit('Please provide either a second vcf at --vcf2 argument or a second accession ID at --comp2 argument or both')
	
	# Loading genotype information of accession 1
	sites_number = 0
	dico_comp1 = {}
	file = open(VCF)
	for line in file:
		data = line.split()
		if data:
			# Recording header
			if data[0] == '#CHROM':
				header = list(data)
				if not(COMP1 in header):
					sys.exit('This is embarrassing ... the program exited without finishing because accession name passed in --comp1 was not found in vcf')
			# Printing headers
			elif data[0][0] != '#':
				sites_number += 1
				# getting genotype
				list_geno = recup_geno(data, header, COMP1)
				if list_geno[0] in dico_comp1:
					if list_geno[1] in dico_comp1[list_geno[0]]:
						sys.exit('This is embarrassing ... the program exited without finishing because two variants have the same position')
				else:
					dico_comp1[list_geno[0]] = {}
				dico_comp1[list_geno[0]][list_geno[1]] = sorted(list_geno[2:])
	file.close()
	
	# Comparing information of accession 2 with information of accession 1
	if VCF2 == None:
		VCF2 = VCF
	if COMP2 == None:
		COMP2 = COMP1
	# initializing values
	specific_comp1 = 0
	specific_comp2 = 0
	shared = 0
	identical_shared = 0
	file = open(VCF2)
	for line in file:
		data = line.split()
		if data:
			# Recording header
			if data[0] == '#CHROM':
				header = list(data)
				if not(COMP2 in header):
					if VCF2 == None:
						sys.exit('This is embarrassing ... the program exited without finishing because accession name passed in --comp2 was not found in vcf')
					else:
						sys.exit('This is embarrassing ... the program exited without finishing because accession name passed in --comp2 was not found in vcf2')
			# Printing headers
			elif data[0][0] != '#':
				# getting genotype
				list_geno = recup_geno(data, header, COMP2)
				if list_geno[0] in dico_comp1:
					if list_geno[1] in dico_comp1[list_geno[0]]:
						shared += 1
						if dico_comp1[list_geno[0]][list_geno[1]] == sorted(list_geno[2:]):
							identical_shared += 1
					else:
						specific_comp2 += 1
				else:
					specific_comp2 += 1
	file.close()
	specific_comp1 = sites_number - shared
	sys.stdout.write(COMP1+' specific variant site number in '+VCF+': '+str(specific_comp1)+'\n')
	sys.stdout.write(COMP2+' specific variant site number in '+VCF2+': '+str(specific_comp2)+'\n')
	sys.stdout.write('Shared variant site number: '+str(shared)+'\n')
	sys.stdout.write('Identical calling shared variant site number: '+str(identical_shared)+'\n')
	sys.stdout.write('Proportion of identical calling shared variant site number: '+str((float(identical_shared)/shared)*100)+'%\n')
	sys.stdout.write('Non identical calling shared variant site number: '+str(shared - identical_shared)+'\n')

def AddRefToCall(VCF, PREFIX, COV):
	
	"""
		Add a last genotype call column to the vcf.
		
		Added line have ref_silico_call name (considered as haploid)
		
		:param VCF: A vcf file.
		:type VCF: vcf
		:param PREFIX: Prefix of the output file.
		:type PREFIX: str
		:return: A vcf file with the added column at the end of the file.
		:rtype: void
	"""
	
	outfile = open(PREFIX+'_add_ref.vcf','w')
	
	# Opening vcf
	file = open(VCF)
	# Reading file line by line
	for line in file:
		data = line.split()
		if data:
			# recording header
			if data[0] == '#CHROM':
				header = list(data)
				header.append('ref_silico_call')
				outfile.write('\t'.join(header)+'\n')
			# not working on headers of the file
			elif data[0][0] == '#':
				outfile.write(line)
			else:
				# getting reference allele
				flag_format = data[header.index('FORMAT')].split(':')
				alt_allele = data[header.index('ALT')].split(',')
				new_col = []
				for n in flag_format:
					if n == 'GT':
						new_col.append('0')
					elif n == 'DP':
						new_col.append(str(COV))
					elif n == 'AD':
						mot_AD = [str(COV)]
						for k in alt_allele:
							mot_AD.append('0')
						new_col.append(','.join(mot_AD))
					else:
						new_col.append('.')
				data.append(':'.join(new_col))
				outfile.write('\t'.join(data)+'\n')
				
def Record_CDS_and_UTR(GFF3, dico_gene_gff3):
	
	"""
		Record in a dictionnary CDS and UTR positions.
		
		:param GFF3: A gff3 file.
		:type GFF3: gff3
		:param dico_gene_gff3: A dictionary.
		:type dico_gene_gff3: dictionary
		:return: Fill dictionary containing for each gene CDS, UTR and INTRON position.
		:rtype: void
	"""
	
	# initializing dictionaries
	dico_gff3 = {}
	dico_parent = {}
	# recording for each gene CDS and exon position
	file = open(GFF3)
	for line in file:
		data = line.replace('\n','').split('\t')
		if data:
			if data[0][0] != "#":
				# recording gene in dictionary
				if data[2] == 'gene':
					dico = gffline2dic(data[8])
					dico_gff3[dico['ID']] = {}
					dico_gff3[dico['ID']]['exon'] = []
					dico_gff3[dico['ID']]['CDS'] = []
					dico_gff3[dico['ID']]['Start'] = int(data[3])
					dico_gff3[dico['ID']]['End'] = int(data[4])
					dico_gff3[dico['ID']]['chr'] = data[0]
				elif data[2] == 'mRNA':
					dico = gffline2dic(data[8])
					dico_parent[dico['ID']] = dico['Parent']
				elif data[2] == 'CDS':
					dico = gffline2dic(data[8])
					nom_gene = dico_parent[dico['Parent']]
					dico_gff3[nom_gene]['CDS'].append([int(data[3]), int(data[4])])
				elif data[2] == 'exon':
					dico = gffline2dic(data[8])
					nom_gene = dico_parent[dico['Parent']]
					dico_gff3[nom_gene]['exon'].append([int(data[3]), int(data[4])])
	
	# recording each type (CDS, UTR and intron) their positions
	for n in dico_gff3:
		if not(dico_gff3[n]['chr']) in dico_gene_gff3:
			dico_gene_gff3[dico_gff3[n]['chr']] = {}
			dico_gene_gff3[dico_gff3[n]['chr']]['CDS'] = set()
			dico_gene_gff3[dico_gff3[n]['chr']]['UTR'] = set()
			dico_gene_gff3[dico_gff3[n]['chr']]['intron'] = set()
		# recording CDS positions
		for k in dico_gff3[n]['CDS']:
			debut = k[0]
			while debut <= k[1]:
				dico_gene_gff3[dico_gff3[n]['chr']]['CDS'].add(debut)
				debut += 1
		# recording UTR positions
		for k in dico_gff3[n]['exon']:
			debut = k[0]
			while debut <= k[1]:
				if not(debut in dico_gene_gff3[dico_gff3[n]['chr']]['CDS']):
					dico_gene_gff3[dico_gff3[n]['chr']]['UTR'].add(debut)
				debut += 1
		# recording intron positions
		debut = dico_gff3[n]['Start']
		while debut <=  dico_gff3[n]['End']:
			if not(debut in dico_gene_gff3[dico_gff3[n]['chr']]['CDS']) and not(debut in dico_gene_gff3[dico_gff3[n]['chr']]['UTR']):
				dico_gene_gff3[dico_gff3[n]['chr']]['intron'].add(debut)
			debut += 1
	file.close()
	
def gffline2dic(LINE):
	
	"""
		Return a dictionnary of the information contained in 9th column of a gff3 file
		
		:param LINE: The 9th column of a gff3 file.
		:type LINE: str
		:return: Dictionary containing each element contained in the line.
		:rtype: dictionary
	"""
	# initializing dictionary
	dico_line = {}
	# recording information in dictionary
	to_split = LINE.split(';')
	for n in to_split:
		data = n.split('=')
		dico_line[data[0]] = data[1]
	return dico_line

def get_CDS_and_UTR_position(GFF3, dico_gff3_position):
	
	"""
		Return a dictionnary containing for each position in UTR or CDS, the gene it correspond to
		
		:param LINE: The 9th column of a gff3 file.
		:type LINE: str
		:param dico_gene_gff3: A dictionary.
		:type dico_gene_gff3: dictionary
		:return: Dictionary containing each element contained in the line.
		:rtype: dictionary
	"""
	
	dico_parent = {}
	file = open(GFF3)
	for line in file:
		data = line.replace('\n','').split('\t')
		if data:
			if data[0][0] != "#":
				# recording gene in dictionary
				if data[2] == 'mRNA':
					dico = gffline2dic(data[8])
					dico_parent[dico['ID']] = dico['Parent']
				elif data[2] == 'CDS':
					dico = gffline2dic(data[8])
					nom_gene = dico_parent[dico['Parent']]
					debut = int(data[3])
					if not(data[0] in dico_gff3_position):
						dico_gff3_position[data[0]] = {}
					while debut <= int(data[4]):
						if not(debut in dico_gff3_position[data[0]]):
							dico_gff3_position[data[0]][debut] = set()
						dico_gff3_position[data[0]][debut].add(nom_gene)
						debut += 1
				elif data[2] == 'exon':
					dico = gffline2dic(data[8])
					nom_gene = dico_parent[dico['Parent']]
					debut = int(data[3])
					if not(data[0] in dico_gff3_position):
						dico_gff3_position[data[0]] = {}
					while debut <= int(data[4]):
						if not(debut in dico_gff3_position[data[0]]):
							dico_gff3_position[data[0]][debut] = set()
						dico_gff3_position[data[0]][debut].add(nom_gene)
						debut += 1
	file.close()

def CalcAllelicIdent(VCF, PREFIX):
	
	"""
		Return a matrix file containing accession identity
		
		:param VCF: A vcf file.
		:type VCF: vcf
		:param PREFIX: Prefix of the output file.
		:type PREFIX: str
		:return: A matrix file containing accession identity.
		:rtype: void
	"""
	
	i = 0
	# reading file
	file = open(VCF)
	for line in file:
		data = line.split()
		if data:
			# Recording header
			if data[0] == '#CHROM':
				header = list(data)
				sur_accession = 0
				list_accession = []
				for n in header:
					if n == 'FORMAT':
						sur_accession = 1
					elif sur_accession:
						list_accession.append(n)
				# creating dictionaries
				dico_nb_all = {}
				dico_nb_all_com = {}
				list_vide = []
				for n in list_accession:
					list_vide.append(0)
				for n in list_accession:
					dico_nb_all[n] = list(list_vide)
					dico_nb_all_com[n] = list(list_vide)
			elif data[0][0] != "#":
				i += 1
				print(i)
				acc_treated = set()
				for acc1 in list_accession:
					genotype1 = recup_geno(data, header, acc1)
					genotype1_set = set(genotype1[2:])
					genotype1_set.discard('NA')
					nb_all = 0
					nb_all_com = 0
					if not(len(genotype1_set) == 0):
						for acc2 in list_accession:
							if not(acc2 in acc_treated):
								genotype2 = recup_geno(data, header, acc2)
								genotype2_set = set(genotype2[2:])
								genotype2_set.discard('NA')
								if not(len(genotype2_set) == 0):
									genotype_tot = genotype1_set.union(genotype2_set)
									genotype_com = genotype1_set.intersection(genotype2_set)
									nb_all = len(genotype_tot)
									nb_all_com = len(genotype_com)
								dico_nb_all[acc1][list_accession.index(acc2)] += nb_all
								dico_nb_all[acc2][list_accession.index(acc1)] += nb_all
								dico_nb_all_com[acc1][list_accession.index(acc2)] += nb_all_com
								dico_nb_all_com[acc2][list_accession.index(acc1)] += nb_all_com
					acc_treated.add(acc1)
	file.close()
	
	# Printing results to files
	outfile = open(PREFIX+'_ident.mat','w')
	outfile.write('ID\t'+'\t'.join(list_accession)+'\n')
	for acc1 in list_accession:
		outfile.write(acc1)
		for acc2 in list_accession:
			if dico_nb_all[acc1][list_accession.index(acc2)] == 0:
				outfile.write('\tNA')
				print(acc1, acc2)
			else:
				outfile.write('\t'+str(float(dico_nb_all_com[acc1][list_accession.index(acc2)])/float(dico_nb_all[acc1][list_accession.index(acc2)])))
		outfile.write('\n')
	outfile.close()

def FormatForPCA(VCF, NAMES, PREFIX, GROUP, AXIS, MULTYPE, DGROUP, MAT):
	
	"""
		Return a matrix file containing accession identity
		
		:param VCF: A vcf file.
		:type VCF: vcf
		:param NAMES: Names of accession tob treate in the PCA.
		:type NAMES: str
		:param PREFIX: Prefix of the output file.
		:type PREFIX: str
		:param GROUP: A file containing two sections: A section[group] with in col 1 accession name ; col 2 group (UN for unknown group, Accessions with unknown group could be ommited). A section [color], that define for each group a color for pca drawing.
		:type GROUP: str
		:param AXIS: Axis number to keep at the end of pca.
		:type AXIS: int
		:param MULTYPE: Type of multivariate analysis to perform.
		:type MULTYPE: str
		:param DGROUP: A list of group to remove from the analysis.
		:type DGROUP: str
		:param MAT: File name of allele coding and grouping
		:type MAT: str
		:return: (1) A matrix file containing accession identity, (2) A script file that will be used by R to run multivariate analysis.
		:rtype: void
	"""
	
	# Recording alleles to remove from the analysis
	ListeAlleleToRemove = set()
	if DGROUP:
		if not(MAT):
			sys.exit('Please, provide a matrix of clustered alleles to --mat option as grouped alleles should be removed according to --dGroup option')
		ListGroupToRemove = DGROUP.split(':')
		file = open(MAT)
		header = file.readline().split('\t')
		TailleHeader = len(header)
		for line in file:
			data = line.split('\t')
			if data:
				if 'K-mean_GROUP' in header:
					if data[header.index('K-mean_GROUP')] in ListGroupToRemove:
						ListeAlleleToRemove.add(':'.join(data[0].split(':')[0:3]))
				elif 'GROUP' in header:
					if data[header.index('GROUP')] in ListGroupToRemove:
						ListeAlleleToRemove.add(':'.join(data[0].split(':')[0:3]))
				else:
					sys.exit('There must be a probleme in formating of the matrix of clustered alleles passed to --mat option. Please see the tutorial to have an idea of the format.')
	
	# Recording group to plot
	dico_group = {}
	dico_color = {}
	max_len_group_col = 0 # For color plotting at the end of matrix creation
	
	#111111111111111111111111111111111111111111111111122222222222222222222222222222222222222
	# dico_group_number = {}
	#111111111111111111111111111111111111111111111111122222222222222222222222222222222222222
	
	if not(GROUP == None):
		sur_group = False
		sur_color = False
		file = open(GROUP)
		for line in file:
			data = line.split()
			if data:
				if data[0] == '[group]':
					sur_group = True
					sur_color = False
				elif data[0] == '[color]':
					sur_group = False
					sur_color = True
				elif data[1] != 'UN' and sur_group:
					dico_group[data[0]] = data[1]
					
					#111111111111111111111111111111111111111111111111122222222222222222222222222222222222222
					# if not(data[1] in dico_group_number):
						# dico_group_number[data[1]] = 0
					# dico_group_number[data[1]] += 1
					#111111111111111111111111111111111111111111111111122222222222222222222222222222222222222
					
				elif sur_color:
					dico_color[data[0]] = {}
					max_len_group_col = max(len(data[0].split('_')), max_len_group_col)
					for color in data[1].split(':'):
						coding = color.split('=')
						dico_color[data[0]][coding[0]] = float(coding[1])
		file.close()
	# sys.exit(dico_group_number)
	
	# Recording accession to treate
	file = open(NAMES)
	LIST_NAME = []
	for line in file:
		data = line.split()
		if data:
			if data[0] in LIST_NAME:
				sys.exit('This is embarassing... The program exited because there is redundancy in accession name in file passed in --names argument')
			LIST_NAME.append(data[0])
	file.close()
	
	
	# reading file
	color_2_create = set()
	file = open(VCF)
	outfile = open(PREFIX+'_matrix_4_PCA.tab','w')
	outfile.write('\tGROUP')
	nb_variant = 0
	for line in file:
		data = line.split()
		if data:
			# Recording header
			if data[0] == '#CHROM':
				header = list(data)
				for acc in LIST_NAME:
					outfile.write('\t'+acc)
				outfile.write('\n')
			elif data[0][0] != "#":
				# recording variant position
				var_pos = data[header.index('POS')]
				var_chr = data[header.index('#CHROM')]
				# recording allele information
				ref_allele = [data[header.index('REF')]]
				alt_allele = data[header.index('ALT')].split(',')
				set_ref_allele = set(ref_allele)
				set_alt_allele = set(alt_allele)
				set_all_allele = set_alt_allele.union(set_ref_allele)
				AlleleToRemove = set()
				for allele in set_all_allele:
					if var_chr+':'+var_pos+':'+allele in ListeAlleleToRemove:
						AlleleToRemove.add(allele)
				set_all_allele.difference_update(AlleleToRemove)
				# recoding genotype
				dico_coded = {}
				have_missing = False
				for acc in LIST_NAME:
					geno = set(recup_geno(data, header, acc)[2:])
					if 'NA' in geno:
						have_missing = True
					for allele in set_all_allele:
						if not(allele in dico_coded):
							dico_coded[allele] = {}
							dico_coded[allele]['present'] = []
							dico_coded[allele]['absent'] = []
						if allele in geno:
							dico_coded[allele]['present'].append(1)
							dico_coded[allele]['absent'].append(0)
						else:
							dico_coded[allele]['present'].append(0)
							dico_coded[allele]['absent'].append(1)
				if have_missing:
					sys.stdout.write('Site '+var_chr+' '+var_pos+' was removed because missing data have been found\n')
				else:
					for allele in set_all_allele:
						remove = False
						if sum(dico_coded[allele]['present']) == 0:
							remove = True
							sys.stdout.write('Allele '+allele+' at position '+var_chr+' '+var_pos+' was removed because it is absent in all selected accessions\n')
						if sum(dico_coded[allele]['absent']) == 0:
							remove = True
							sys.stdout.write('Allele '+allele+' at position '+var_chr+' '+var_pos+' was removed because it is present in all selected accessions\n')
						if remove:
							del dico_coded[allele]
						else:
							# attributing group to marker
							dico_grouped = set()
							
							#1111111111111111111111111111111111111122222222222222222222222222222222222222
							# dico_count = {}
							# for n in dico_group_number:
								# dico_count[n] = 0
							#1111111111111111111111111111111111111122222222222222222222222222222222222222
							
							for n in LIST_NAME:
								if dico_coded[allele]['present'][LIST_NAME.index(n)] == 1:
									if n in dico_group:
										dico_grouped.add(dico_group[n])
										
										#1111111111111111111111111111111111111122222222222222222222222222222222222222
										# dico_count[dico_group[n]] +=1
										#1111111111111111111111111111111111111122222222222222222222222222222222222222
							
							#11111111111111111111111111111111111111
							# dico_grouped = set()
							# for n in dico_count:
								# if dico_count[n] == dico_group_number[n]:
									# dico_grouped.add(n)
							#11111111111111111111111111111111111111
							
							
							#22222222222222222222222222222222222222
							# dico_grouped2 = set()
							# for n in dico_count:
								# if dico_count[n] == dico_group_number[n]:
									# dico_grouped2.add(n)
							#22222222222222222222222222222222222222
							
							# printing line
							if len(dico_grouped) > max_len_group_col:
								outfile.write(var_chr+':'+var_pos+':'+allele+':P\t'+'UN\t'+'\t'.join(list(map(str, dico_coded[allele]['present'])))+'\n')
							elif len(dico_grouped) == 0:
								outfile.write(var_chr+':'+var_pos+':'+allele+':P\t'+'UN\t'+'\t'.join(list(map(str, dico_coded[allele]['present'])))+'\n')
								
							#22222222222222222222222222222222222222
							# elif len(dico_grouped2) > max_len_group_col:
								# outfile.write(var_chr+':'+var_pos+':'+allele+':P\t'+'UN\t'+'\t'.join(list(map(str, dico_coded[allele]['present'])))+'\n')
							# elif len(dico_grouped2) == 0:
								# outfile.write(var_chr+':'+var_pos+':'+allele+':P\t'+'UN\t'+'\t'.join(list(map(str, dico_coded[allele]['present'])))+'\n')
							#22222222222222222222222222222222222222
							
							else:
								color_2_create.add('_'.join(sorted(list(dico_grouped))))
								outfile.write(var_chr+':'+var_pos+':'+allele+':P\t'+'_'.join(sorted(list(dico_grouped)))+'\t'+'\t'.join(list(map(str, dico_coded[allele]['present'])))+'\n')
							
							# attributing group to marker
							dico_grouped = set()
							
							#1111111111111111111111111111111111111122222222222222222222222222222222222222
							# dico_count = {}
							# for n in dico_group_number:
								# dico_count[n] = 0
							#1111111111111111111111111111111111111122222222222222222222222222222222222222
							
							for n in LIST_NAME:
								if dico_coded[allele]['absent'][LIST_NAME.index(n)] == 1:
									if n in dico_group:
										dico_grouped.add(dico_group[n])
										
										#1111111111111111111111111111111111111122222222222222222222222222222222222222
										# dico_count[dico_group[n]] +=1
										#1111111111111111111111111111111111111122222222222222222222222222222222222222
							
							#11111111111111111111111111111111111111
							# dico_grouped = set()
							# for n in dico_count:
								# if dico_count[n] == dico_group_number[n]:
									# dico_grouped.add(n)
							#11111111111111111111111111111111111111
							
							
							#22222222222222222222222222222222222222
							# dico_grouped2 = set()
							# for n in dico_count:
								# if dico_count[n] == dico_group_number[n]:
									# dico_grouped2.add(n)
							#22222222222222222222222222222222222222
							
							# printing line
							if len(dico_grouped) > max_len_group_col:
								outfile.write(var_chr+':'+var_pos+':'+allele+':A\t'+'UN\t'+'\t'.join(list(map(str, dico_coded[allele]['absent'])))+'\n')
							elif len(dico_grouped) == 0:
								outfile.write(var_chr+':'+var_pos+':'+allele+':A\t'+'UN\t'+'\t'.join(list(map(str, dico_coded[allele]['absent'])))+'\n')
							
							#22222222222222222222222222222222222222
							# elif len(dico_grouped2) > max_len_group_col:
								# outfile.write(var_chr+':'+var_pos+':'+allele+':A\t'+'UN\t'+'\t'.join(list(map(str, dico_coded[allele]['absent'])))+'\n')
							# elif len(dico_grouped2) == 0:
								# outfile.write(var_chr+':'+var_pos+':'+allele+':A\t'+'UN\t'+'\t'.join(list(map(str, dico_coded[allele]['absent'])))+'\n')
							#22222222222222222222222222222222222222
							
							else:
								color_2_create.add('_'.join(sorted(list(dico_grouped))))
								outfile.write(var_chr+':'+var_pos+':'+allele+':A\t'+'_'.join(sorted(list(dico_grouped)))+'\t'+'\t'.join(list(map(str, dico_coded[allele]['absent'])))+'\n')
							nb_variant += 1
							del dico_coded[allele]
	sys.stdout.write('*****\nA total of '+str(nb_variant)+' alleles were kept for analysis\n*****\n')
	
	# Creating color to plot
	list_color_2_plot = ['UN']
	dico_color['UN'] = {}
	dico_color['UN']['red'] = 0.6
	dico_color['UN']['green'] = 0.6
	dico_color['UN']['blue'] = 0.6
	dico_color['UN']['alpha'] = 0.5
	for n in color_2_create:
		list_color_2_plot.append(n)
		if not(n in dico_color):
			dico_color[n] = {}
			dico_color[n]['red'] = 0.6
			dico_color[n]['green'] = 0.6
			dico_color[n]['blue'] = 0.6
			dico_color[n]['alpha'] = 0.5
	list_color_2_plot.sort()
	
	# Creating the R script for PCA analysis
	outfile = open(PREFIX+'_multivariate.R','w')
	
	# Creating the vector for color printing
	for n in list_color_2_plot:
		outfile.write('color_'+n+' = rgb(red='+str(dico_color[n]['red'])+', green='+str(dico_color[n]['green'])+', blue='+str(dico_color[n]['blue'])+', alpha='+str(dico_color[n]['alpha'])+')\n')
	outfile.write('gcol = c(color_'+', color_'.join(list_color_2_plot)+')\n')
	# Loading ade4 library
	outfile.write('library(ade4)\n')
	# Loading table
	outfile.write('data = read.table(file="'+PREFIX+'_matrix_4_PCA.tab")\n')
	# Record group type in the group_type variable
	outfile.write('group_type = factor(as.matrix(data[,1]))\n')
	# Removing the group line from the matrix
	outfile.write('tableau = data[,-1]\n')
	# Transpose matrix and recoding in dataframe
	outfile.write('tableau = as.data.frame(t(tableau))\n')
	# Doing PCA analysis
	if MULTYPE == 'coa':
		outfile.write('multivariate_analysis=dudi.coa(tableau, scannf = F, nf = '+str(AXIS)+')\n')
	elif MULTYPE == 'pca':
		outfile.write('multivariate_analysis=dudi.pca(tableau, scannf = F, center = F, scale = F, nf = '+str(AXIS)+')\n')
	elif MULTYPE == 'pca_normed':
		outfile.write('multivariate_analysis=dudi.pca(tableau, scannf = F, center = T, scale = T, nf = '+str(AXIS)+')\n')
	else:
		sys.exit('The program exited without finishing because multivariate analysis provided in --mulType argument is not recognised: '+MULTYPE)
	# Print barplot of synthetic variables inertia and print inertia value
	outfile.write('inertia.dudi(multivariate_analysis)$TOT$ratio*100\n')
	outfile.write('pdf("'+PREFIX+'_inertia.pdf")\n')
	outfile.write('barplot(multivariate_analysis$eig/sum(multivariate_analysis$eig)*100, xlab="Synthetic variables", ylab="% inertia")\n')
	outfile.write('dev.off()\n')
	# Print individuals and variables on new axis
	valeur = range(AXIS)
	combinaisons = []
	n_valeur = set()
	for n in valeur:
		n_valeur.add(n)
		j_valeur = set()
		for j in valeur:
			j_valeur.add(j)
			if not(j in n_valeur):
				combinaisons.append(sorted([n,j]))
	for comb in combinaisons:
		outfile.write('pdf("'+PREFIX+'_axis_'+str(comb[0]+1)+'_vs_'+str(comb[1]+1)+'.pdf")\n')
		outfile.write('par(mfrow=c(2,2))\n')
		outfile.write('s.label(multivariate_analysis$li, xax='+str(comb[0]+1)+', yax='+str(comb[1]+1)+')\n')
		outfile.write('s.class(dfxy = multivariate_analysis$co, fac = group_type, col = gcol, xax='+str(comb[0]+1)+', yax='+str(comb[1]+1)+')\n')
		outfile.write('dev.off()\n')
	for comb in combinaisons:
		outfile.write('pdf("'+PREFIX+'_axis_'+str(comb[0]+1)+'_vs_'+str(comb[1]+1)+'_accessions.pdf")\n')
		outfile.write('s.label(multivariate_analysis$li, xax='+str(comb[0]+1)+', yax='+str(comb[1]+1)+')\n')
		outfile.write('dev.off()\n')
	# Printing variables coordinates
	outfile.write('write.table(multivariate_analysis$co, file = "'+PREFIX+'_variables_coordinates.tab")\n')
	# Printing normalized variable coordinates
	outfile.write('scaled.dat <- scale(multivariate_analysis$co)\n')
	outfile.write('write.table(scaled.dat, file = "'+PREFIX+'_variables_coordinates_scaled.tab")\n')
	# Printing individuals coordinates
	outfile.write('write.table(multivariate_analysis$li, file = "'+PREFIX+'_individuals_coordinates.tab")\n')
	outfile.close()

def Draw3dPlot(VARCOORD, DAXES, MAT, GROUP, DGROUP):
	
	"""
		Perform a 3d plot of the variables in the requested axis
		
		:param VARCOORD: The tabulated file of variables coordinates in new axis.
		:type VARCOORD: str
		:param DAXES: Axes to draw in 3d plot. Axis are separated by =.
		:type DAXES: str
		:param MAT: The tabulated file containing variant allele encoded with: in column 1, variant name and column 2, variant group.
		:type MAT: str
		:param GROUP: A file containing two sections: A section[group] with in col 1 accession name ; col 2 group (UN for unknown group, Accessions with unknown group could be ommited). A section [color], that define for each group a color for pca drawing.
		:type GROUP: str
		:param DGROUP: A string containing groups to draw.
		:type DGROUP: str
		:return: Perform 3d plod.
		:rtype: void
	"""
	
	dico_group = {}
	dico_color = {}
	
	# Recording group to plot
	dico_color['UN'] = {}
	dico_color['UN']['color'] = (0.5, 0.5, 0.5, 0.5)
	dico_color['UN']['group'] = set()
	dico_color['UN']['to_plot'] = []
	max_len_group_col = 0 # For color plotting at the end of matrix creation
	if not(GROUP == None):
		sur_group = False
		sur_color = False
		file = open(GROUP)
		for line in file:
			data = line.split()
			if data:
				if data[0] == '[group]':
					sur_group = True
					sur_color = False
				elif data[0] == '[color]':
					sur_group = False
					sur_color = True
				elif sur_color:
					dico_color[data[0]] = {}
					max_len_group_col = max(len(data[0].split('_')), max_len_group_col)
					for color in data[1].split(':'):
						coding = color.split('=')
						dico_color[data[0]][coding[0]] = float(coding[1])
					# creating color tuple
					dico_color[data[0]]['color'] = (dico_color[data[0]]['red'], dico_color[data[0]]['green'], dico_color[data[0]]['blue'], dico_color[data[0]]['alpha'])
					dico_color[data[0]]['to_plot'] = []
		file.close()
	
	if not(MAT == None):
		# Opening matrix file
		file = open(MAT)
		# recording header
		header = file.readline().split()
		# Recording markers to draw and their group to a dictionnary
		for line in file:
			data = line.split()
			if data:
				dico_group[data[0]] = data[1]
		file.close()
	
	# recording axis to draw
	axis = list(map(int,DAXES.split(':')))
	if len(axis) != 3:
		sys.exit('This is embarrassing... The program exited without finishing because exacly 3 axis should be provided in --dAxes argument and it is not the case')
	
	# Opening coordinate file
	file = open(VARCOORD)
	# recording header
	header = file.readline().replace('"','').split()
	# recording coordinates
	for line in file:
		data = line.split()
		if data:
			allele = data[0].replace('"','').replace('.',':')
			if allele in dico_group:
				if dico_group[allele] in dico_color:
					dico_color[dico_group[allele]]['to_plot'].append(list(map(float, [data[axis[0]],data[axis[1]], data[axis[2]]])))
				else:
					dico_color['UN']['to_plot'].append(list(map(float, [data[axis[0]],data[axis[1]], data[axis[2]]])))
			elif len(dico_group) == 0:
				dico_color['UN']['to_plot'].append(list(map(float, [data[axis[0]],data[axis[1]], data[axis[2]]])))
	file.close()
	
	# Removing dots that should not be drawn if --dGroup is filled
	if not(DGROUP == None):
		group2draw = DGROUP.split(':')
		group_2_remove = set()
		for n in dico_color:
			if not(n in group2draw):
				group_2_remove.add(n)
		for n in group_2_remove:
			del dico_color[n]
	else:
		if len(dico_color['UN']['to_plot']) == 0:
			del dico_color['UN']
	
	# Just some statistics
	for n in dico_color:
		sys.stdout.write('Group '+n+' contains '+str(len(dico_color[n]['to_plot']))+' dots\n')
		dico_color[n]['to_plot'] = numpy.array(dico_color[n]['to_plot'])
	
	# Drawing plot
	Draw_Plot(dico_color, [header[axis[0]-1], header[axis[1]-1], header[axis[2]-1]])

def Draw2dPlot(VARCOORD, DAXES, MAT, GROUP, DGROUP, PREFIX):
	
	"""
		Perform a 3d plot of the variables in the requested axis
		
		:param VARCOORD: The tabulated file of variables coordinates in new axis.
		:type VARCOORD: str
		:param DAXES: Axes to draw. Axis are separated by =.
		:type DAXES: str
		:param MAT: The tabulated file containing variant allele encoded with: in column 1, variant name and column 2, variant group.
		:type MAT: str
		:param GROUP: A file containing at least a section [color], that define for each group a color for pca drawing (RGB+alpha).
		:type GROUP: str
		:param DGROUP: A string containing groups to draw.
		:type DGROUP: str
		:return: Perform 3d plod.
		:rtype: void
	"""
	
	# recording axis to draw
	axis = list(map(int,DAXES.split(':')))
	
	dico_group = {}
	dico_color = {}
	
	# Recording group to plot
	dico_color['UN'] = {}
	dico_color['UN']['color'] = (0.5, 0.5, 0.5, 0.5)
	dico_color['UN']['group'] = set()
	dico_color['UN']['to_plot'] = []
	max_len_group_col = 0 # For color plotting at the end of matrix creation
	if not(GROUP == None):
		sur_group = False
		sur_color = False
		file = open(GROUP)
		for line in file:
			data = line.split()
			if data:
				if data[0] == '[group]':
					sur_group = True
					sur_color = False
				elif data[0] == '[color]':
					sur_group = False
					sur_color = True
				elif sur_color:
					dico_color[data[0]] = {}
					max_len_group_col = max(len(data[0].split('_')), max_len_group_col)
					for color in data[1].split(':'):
						coding = color.split('=')
						dico_color[data[0]][coding[0]] = float(coding[1])
					# creating color tuple
					dico_color[data[0]]['color'] = (dico_color[data[0]]['red'], dico_color[data[0]]['green'], dico_color[data[0]]['blue'], dico_color[data[0]]['alpha'])
					dico_color[data[0]]['to_plot'] = []
		file.close()
	
	if not(MAT == None):
		# Opening matrix file
		file = open(MAT)
		# recording header
		header = file.readline().split()
		# Recording markers to draw and their group to a dictionnary
		for line in file:
			data = line.split()
			if data:
				dico_group[data[0]] = data[1]
		file.close()
	
	# Opening coordinate file
	file = open(VARCOORD)
	# recording header
	header = file.readline().replace('"','').split()
	# recording coordinates
	for line in file:
		data = line.split()
		if data:
			allele = data[0].replace('"','').replace('.',':')
			if allele in dico_group:
				if dico_group[allele] in dico_color:
					dico_color[dico_group[allele]]['to_plot'].append(list(map(float, data[1:])))
				else:
					dico_color['UN']['to_plot'].append(list(map(float, data[1:])))
			elif len(dico_group) == 0:
				dico_color['UN']['to_plot'].append(list(map(float, data[1:])))
	file.close()
	
	# Removing dots that should not be drawn if --dGroup is filled
	if not(DGROUP == None):
		group2draw = DGROUP.split(':')
		group_2_remove = set()
		for n in dico_color:
			if not(n in group2draw):
				group_2_remove.add(n)
		for n in group_2_remove:
			del dico_color[n]
	else:
		if len(dico_color['UN']['to_plot']) == 0:
			del dico_color['UN']
	
	for n in dico_color:
		sys.stdout.write('Group '+n+' contains '+str(len(dico_color[n]['to_plot']))+' dots\n')
		dico_color[n]['to_plot'] = numpy.array(dico_color[n]['to_plot'])
	
	# Drawing 2d plot
	fait = set()
	for axis1 in axis:
		if not(axis1 in fait):
			fait.add(axis1)
			for axis2 in axis:
				if not(axis2 in fait):
					Draw_2d_Plot(dico_color, [axis1-1, axis2-1], PREFIX+'_axis'+'_vs_axis'.join(list(map(str,[axis1, axis2])))+'.png', 'Grouping axis'+' vs axis'.join(list(map(str,[axis1, axis2]))))

def Draw_Plot(DICO, HEAD):
	
	"""
		Draw the plot based on a dictionnary
		
		:param DICO: A dictionary containing for each group (key) a sub dictionary with an array corresponding to coordinates in 'to_plot' key and group color in 'color' key.
		:type DICO: dictionary
		:param HEAD: A list with three elements corresponding to axes names.
		:type HEAD: list
		:return: return a 3d plot (interactive).
		:rtype: void
	"""
	
	# Drawing plot
	fig = plt.figure(figsize=(11,11))
	ax = fig.add_subplot(111, projection='3d')
	legend = []
	legend_name = []
	for n in DICO:
		dot_legend, = ax.plot(DICO[n]['to_plot'][:,0], DICO[n]['to_plot'][:,1], DICO[n]['to_plot'][:,2], 'o', markersize=5, color=DICO[n]['color'])
		legend.append(dot_legend)
		legend_name.append('group '+str(n))
	ax.set_xlabel(HEAD[0])
	ax.set_ylabel(HEAD[1])
	ax.set_zlabel(HEAD[2])
	plt.title('Alleles on synthetics axes')
	ax.legend(legend,legend_name)
	plt.show()

def Draw_2d_Plot(DICO, HEAD, NAME, TITLE):
	
	"""
		Draw the plot based on a dictionnary
		
		:param DICO: A dictionary containing for each group (key) a sub dictionary with an array corresponding to coordinates in 'to_plot' key and group color in 'color' key.
		:type DICO: dictionary
		:param HEAD: A list with 2 elements corresponding to axes number to draw.
		:type HEAD: list
		:param NAME: Output file name.
		:type NAME: str
		:param TITLE: Graphe title.
		:type TITLE: str
		:return: draw a 2 plot png file.
		:rtype: void
	"""
	
	# Drawing plot
	fig = plt.figure(figsize=(13,10))
	fig.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05)
	ax = plt.subplot2grid((1,13),(0,0), colspan=10, rowspan=1)
	legend = []
	legend_name = []
	ListeOrdered = sorted(list(DICO.keys()))
	for n in ListeOrdered:
		dot_legend, = ax.plot(DICO[n]['to_plot'][:,HEAD[0]], DICO[n]['to_plot'][:,HEAD[1]], 'o', markersize=5, color=DICO[n]['color'])
		legend.append(dot_legend)
		legend_name.append('group '+str(n))
	ax.set_xlabel('Axis'+str(HEAD[0]+1))
	ax.set_ylabel('Axis'+str(HEAD[1]+1))
	if TITLE != None:
		plt.title(TITLE)
	ax.legend(legend,legend_name, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
	fig.savefig(NAME)
	plt.close(fig)

def get_max_centro_move(DIST):
	
	"""
		return the maximal value found in the diagonal of the matrix
		
		:param DIST: An array containing by line dot coordinate.
		:type DIST: array
		:return: return a float corresponding to the maximal movment in the matrix.
		:rtype: float
	"""
	
	if DIST.shape[0] != DIST.shape[1]:
		sys.exit('A problem has been encountered when clustering, a cluster has been lost during computation')
	maxi_move = 0
	for i in range(DIST.shape[0]):
		maxi_move = max(maxi_move, DIST[i,i])
	return maxi_move

def Calc_centroid(MY_ARRAY, CLUSTER):
	
	"""
		Calculate new centroid position
		
		:param MY_ARRAY: An array containing by line dot coordinate.
		:type MY_ARRAY: array
		:param CLUSTER: An array with presence absence of dots in centroids.
		:type CLUSTER: array
		:return: return an array with new centroid coordinates.
		:rtype: array
	"""
	
	# recalculating the cluster matrix (each row is divided by its sum)
	normalized_matrix = []
	
	for i in range(CLUSTER.shape[0]):
		normalized_matrix.append(CLUSTER[i,:]/float(numpy.sum(CLUSTER[i,:])))
	normalized_matrix = numpy.matrix(numpy.array(normalized_matrix))
		
	# Calculating new centroid coordinates
	return numpy.array(normalized_matrix * numpy.matrix(MY_ARRAY))

def cluster_around_centroid(DIS_MAT):
	"""
		Create a (nb(dot),nb(centroid)) array with presence absence of dots in centroids.
		
		:param DIS_MAT: The distance matrix between dot and centroids.
		:type DIS_MAT: array
		:return: return an array with presence absence of dots in centroids.
		:rtype: array
	"""
	
	# Creating empty list
	grouped_array = []
		
	# dots clustering
	for i in range(DIS_MAT.shape[0]):
		list = [0]*DIS_MAT.shape[1]
		list[numpy.random.choice(numpy.where(DIS_MAT[i,:] == numpy.amin(DIS_MAT[i,:]))[0])] = 1
		grouped_array.append(list)
	
	# converting list to array
	grouped_array = numpy.transpose(numpy.array(grouped_array))
	return grouped_array

def calc_euclid_dist(INDI, CENTROID):
	
	"""
		Calculate euclidiane distance between individuals and centroids
		
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

def kMean(ARRAY, NGROUP, RANDOM):
	
	"""
		k-mean clustering
		
		:param ARRAY: An array containing by line dot coordinate.
		:type ARRAY: array
		:param NGROUP: Group number to perform.
		:type NGROUP: str
		:param RANDOM: An array of random centroids. If value == 'empty', the random step is performed in this fonction
		:type RANDOM: array
		:return: return a list containing as first element an array with final centroid position (row = centroid, col = centroid coordinates) and second element is an array with presence absence of dots in centroids (row = centroid, lines = dot).
		:rtype: void
	"""
	
	if RANDOM == 'empty':
		#---------------
		#Initialisation
		#---------------
		# Picking k value at random
		centroid = []
		randomed_set = set()
		while len(randomed_set) < NGROUP:
			random_dot = random.randint(0,len(ARRAY)-1)
			if not(':'.join(list(map(str, list(ARRAY[random_dot,:])))) in randomed_set):
				randomed_set.add(':'.join(list(map(str, list(ARRAY[random_dot,:])))))
				centroid.append(ARRAY[random_dot,:])
		centroid = numpy.array(centroid)
	else:
		centroid = RANDOM
	
	# Calculating euclidian distance
	dist_matrix = calc_euclid_dist(ARRAY, centroid)
	
	# Clustering
	cluster = cluster_around_centroid(dist_matrix)
	# print cluster
	
	# Calculating centroid position
	centroid_coordinates = Calc_centroid(ARRAY, cluster)
	
	# Calculating distance between centroids
	centroid_distance_matrix = calc_euclid_dist(centroid, centroid_coordinates)
	
	# Getting maximal movment between centroids
	max_move = get_max_centro_move(centroid_distance_matrix)
	
	#----------------------------
	#Iterrating until convergence
	#----------------------------
	while max_move > 0:
		centroid = centroid_coordinates
		# Calculating euclidian distance
		dist_matrix = calc_euclid_dist(ARRAY, centroid)
		
		# Clustering
		cluster = cluster_around_centroid(dist_matrix)
		# print cluster
		
		# Calculating centroid position
		centroid_coordinates = Calc_centroid(ARRAY, cluster)
		
		# Calculating distance between centroids
		centroid_distance_matrix = calc_euclid_dist(centroid, centroid_coordinates)
		
		# Getting maximal movment between centroids
		max_move = get_max_centro_move(centroid_distance_matrix)
	
	return [centroid_coordinates, cluster]

def look4corres(ARRAY):
	
	"""
		Look for line and column correspondance in an distance matrix (return the minimal value position for each row in the array)
		
		:param ARRAY: An array containing by line dot coordinate.
		:type ARRAY: array
		:return: return list with N elements. The ieme element correspond to position of the minimal value in the ieme raw in the array.
		:rtype: list
	"""
	
	list_cor = []
	for i in range(ARRAY.shape[0]):
		list_cor.append(numpy.random.choice(numpy.where(ARRAY[i,:] == numpy.amin(ARRAY[i,:]))[0]))
	
	return list_cor

def group_centroids_pseudo_kMean(DICO_CENTROIDS, PREFIX):
	
	"""
		Group centroids of each iteration based on kMean approach
		
		:param DICO_CENTROIDS: A dictionary containing each final centroids of each iterations coordinates.
		:type DICO_CENTROIDS: dictionary
		:param PREFIX: Prefix of the output file.
		:type PREFIX: str
		:return: return a list with [0] = list of final centroid position and [2] a dictionary with key = iteration number, list of group correspondance.
		:rtype: dictionary
	"""
	
	# Selecting first centroids before iterations
	# 1 - First iterartion
	mean_centroid = False
	for n in DICO_CENTROIDS:
		if not(mean_centroid):
			cluster_used = 1
			mean_centroid = {}
			for i in range(DICO_CENTROIDS[n].shape[0]):
				mean_centroid[i] = numpy.array(DICO_CENTROIDS[n][i,:])
		else:
			new_centroid = []
			for i in range(len(mean_centroid)):
				new_centroid.append(mean_centroid[i]/float(cluster_used))
			new_centroid = numpy.array(new_centroid)
			# Calculating distance between centroids
			centroid_distance_matrix = calc_euclid_dist(new_centroid, DICO_CENTROIDS[n])
			# Getting which centroid correspond to which centroids (raw -central array, col n[0])
			centro_corres = look4corres(centroid_distance_matrix)
			if len(set(centro_corres)) != len(centro_corres):
				if cluster_used == 1:
					cluster_used = 1
					mean_centroid = {}
					for i in range(DICO_CENTROIDS[n].shape[0]):
						mean_centroid[i] = numpy.array(DICO_CENTROIDS[n][i,:])
			else:
				cluster_used += 1
				# Recording additional centroid positions
				for i in range(len(mean_centroid)):
					mean_centroid[i] = mean_centroid[i] + DICO_CENTROIDS[n][centro_corres[i]]
	new_centroid = []
	for i in range(len(mean_centroid)):
		new_centroid.append(mean_centroid[i]/float(cluster_used))
	new_centroid = numpy.array(new_centroid)
	
	if float(cluster_used)/len(DICO_CENTROIDS) < 1:
		sys.stdout.write('The kMean approache may not be appropiate with the group number provided as only '+str(float(cluster_used)/len(DICO_CENTROIDS)*100)+' percent of centroids where used for random final centroids identification\n')
	
	# 2 - Creation an array containing all centroids of all iterations
	centroid_array = []
	for iter in DICO_CENTROIDS:
		for i in range(DICO_CENTROIDS[iter].shape[0]):
			centroid_array.append(DICO_CENTROIDS[iter][i,:])
	centroid_array = numpy.array(centroid_array)
	
	# 3 - K-mean step on centroids
	centro_cluster = kMean(centroid_array, DICO_CENTROIDS[0].shape[0], new_centroid)
	cluster = centro_cluster[1]
	
	# 4 - Recording in dictionary group correspondance with final clustering
	dico_2_return = {}
	for iter in range(cluster.shape[1]):
		if numpy.sum(cluster[:,iter]) != 1:
			sys.exit('Bug in final estimation of centroids.')
		if not(int(iter/cluster.shape[0])) in dico_2_return:
			dico_2_return[int(iter/cluster.shape[0])] = []
		dico_2_return[int(iter/cluster.shape[0])].append(int(numpy.where(cluster[:,iter] == 1)[0]))
	
	# Creating dico for 3d plot
	dico_plot = {}
	for i in range(cluster.shape[0]):
		dico_plot[i] = {}
		dico_plot[i]['to_plot'] = []
		for allele in range(cluster.shape[1]):
			if cluster[i,allele] == 1:
				dico_plot[i]['to_plot'].append(centroid_array[allele,:])
	
	# Printing centroids grouping in a file
	outfile = open(PREFIX+'_centroid_iteration_grouping.tab','w')
	outfile_coor = open(PREFIX+'_centroid_coordinates.tab','w')
	outfile_coor.write('\taxis'+'\taxis'.join(list(map(str, range(centroid_array.shape[1]))))+'\n')
	outfile.write('\tK-mean_GROUP\n')
	centro_number = 0
	for i in dico_plot:
		dico_plot[i]['to_plot'] = numpy.array(dico_plot[i]['to_plot'])
		sys.stdout.write('Group g'+str(i)+' contained '+str(dico_plot[i]['to_plot'].shape[0])+' centroids\n')
		for j in range(dico_plot[i]['to_plot'].shape[0]):
			centro_number += 1
			outfile_coor.write('centroid'+str(centro_number)+'\t'+'\t'.join(list(map(str, list(dico_plot[i]['to_plot'][j,:]))))+'\n')
			outfile.write('centroid'+str(centro_number)+'\tg'+str(i)+'\n')
	outfile.close()
	outfile_coor.close()
	
	# Creating color
	axes = []
	for n in range(centro_cluster[0].shape[1]):
		axes.append(n)
	nb_color = len(dico_plot)
	color = [(0,0,0,0.5),(1,0,0,0.5),(0,1,0,0.5),(0,0,1,0.5),(1,1,0,0.5),(0,1,1,0.5),(1,0,1,0.5)]
	if nb_color > len(color):
		sys.stdout.write('In the graphe groups will have the same color because there is not enough color defined in the script.\nPlease contact Guillaume MARTIN (guillaume.martin@cirad.fr) for adding colors.\n')
	for n in dico_plot:
		dico_plot[n]['color'] = color[n%len(color)]
		if int(n/len(color)) > 0:
			sys.stdout.write('In the graphe groups '+str(n)+' will have the same color as the '+str(n%len(color))+' group.\n')
	
	
	return [centro_cluster[0], dico_2_return]

def kMean_clust(VARCOORD, DAXES, NGROUP, MAT, PREFIX, THREAD, ITER):
	
	"""
		Perform a k-mean clustering and outpout results
		
		:param VARCOORD: The tabulated file of variables coordinates in new axis.
		:type VARCOORD: str
		:param DAXES: Axes to use for clustering.
		:type DAXES: str
		:param NGROUP: Group number to perform.
		:type NGROUP: str
		:param MAT: The tabulated file containing variant allele encoded.
		:type MAT: str
		:param PREFIX: Prefix of the output file.
		:type PREFIX: str
		:param THREAD: Thread available for calculation.
		:type THREAD: int
		:param ITER: Iteration number for kMean calculation.
		:type ITER: str
		:return: return a file containing clustered variables (alleles).
		:rtype: void
	"""
	
	# recording axis that will be used
	axes = list(map(int, DAXES.split(':')))
	if len(axes) < 2:
		sys.exit('This is embarrassing... The program exited without finishing because at least 2 axes should be filled in --dAxes argument and it is not the case')
	
	#-----------------------------------------------------------
	# Loading coordinates in an array (Absent lines are removed)
	#-----------------------------------------------------------
	my_array = []
	# Opening coordinate file
	file = open(VARCOORD)
	# recording header
	header = file.readline().replace('"','').split()
	# recording coordinates
	col_header = []
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
	my_array = numpy.array(my_array)
	
	#----------------------------------
	#Runing kMean clustering iterations
	#----------------------------------
	centrale_array = False
	dico_centroids = {}
	dico_group = {}
	i = 0
	repet = 0
	thread = 4
	list_job = []
	while i < ITER:
		list_job.append(['kMean', my_array, NGROUP, 'empty', i])
		i += 1
	if list_job:
		pool = mp.Pool(processes=THREAD)
		result = pool.map(main, list_job)
		if 'error' in result:
			raise ErrorValue ('Bug in kMean calculation.')
		else:
			for n in result:
				dico_centroids[repet] = n[0]
				dico_group[repet] = n[1]
				repet += 1
		pool.close()
	del result
	
	#------------------------------------
	#Calculating mean centroids positions
	#------------------------------------
	final_clustering = group_centroids_pseudo_kMean(dico_centroids, PREFIX)
	
	#--------------------------------
	#Grouping dots on final centroids
	#--------------------------------
	# Calculating euclidian distance
	dist_matrix = calc_euclid_dist(my_array, final_clustering[0])
	
	# Clustering
	cluster = cluster_around_centroid(dist_matrix)
	
	#------------------------------------------
	#Calculating dots grouping during iteration
	#------------------------------------------
	# Dictionary creation
	dico_dot_grouping = {}
	list_empty = list(0 for i in range(cluster.shape[0]))
	for i in range(cluster.shape[1]):
		dico_dot_grouping[col_header[i]] = list(list_empty)
	# Recording in dictionary
	for iter in dico_group:
		for group in range(dico_group[iter].shape[0]):
			coresponding_group = final_clustering[1][iter][group]
			in_the_group = list(numpy.where(dico_group[iter][group,:] == 1)[0])
			for dot in in_the_group:
				dico_dot_grouping[col_header[dot]][coresponding_group] += 1
	
	# Recording max group probability in an array
	dico_ambiguous = []
	for dot in dico_dot_grouping:
		dico_ambiguous.append(max(dico_dot_grouping[dot])/float(sum(dico_dot_grouping[dot])))
	min_prob = min(dico_ambiguous)
	max_prob = max(dico_ambiguous)
	dico_ambiguous = numpy.array(dico_ambiguous)
	
	# Printing dot grouping proportion in a file
	# Opening coordinate file
	file = open(MAT)
	# Recording header (not used here)
	header = file.readline().replace('"','').split()
	# Writting to file
	outfile = open(PREFIX+'_kMean_gp_prop.tab','w')
	mot = []
	for i in range(cluster.shape[0]):
		mot.append('g'+str(i))
	outfile.write('\t'+'\t'.join(mot)+'\n')
	for line in file:
		data = line.split()
		if data:
			if data[0] in dico_dot_grouping:
				outfile.write(data[0]+'\t'+'\t'.join(list(map(str, [i/float(sum(dico_dot_grouping[data[0]])) for i in dico_dot_grouping[data[0]]])))+'\n')
	outfile.close()
	
	#------------------------------------------
	#Outputing groups and preparing for 2d plot
	#------------------------------------------
	grouped = {}
	# Creating dico for 2d plot
	dico_plot = {}
	for i in range(cluster.shape[0]):
		dico_plot[i] = {}
		dico_plot[i]['to_plot'] = []
		for allele in range(cluster.shape[1]):
			if cluster[i,allele] == 1:
				nom_allele = col_header[allele]
				if col_header[allele] in grouped:
					sys.exit('bug')
				dico_plot[i]['to_plot'].append(my_array[allele,:])
				grouped[nom_allele] = i
	
	sys.stdout.write('\n**********\n')
	for i in dico_plot:
		dico_plot[i]['to_plot'] = numpy.array(dico_plot[i]['to_plot'])
		sys.stdout.write('Group g'+str(i)+' contained '+str(dico_plot[i]['to_plot'].shape[0])+' dots\n')
	
	# Creating color and creating a color group file
	outfile = open(PREFIX+'_group_color.tab','w')
	outfile.write('[color]\n')
	nb_color = len(dico_plot)
	color = [(0,0,0,0.5),(1,0,0,0.5),(0,1,0,0.5),(0,0,1,0.5),(1,1,0,0.5),(0,1,1,0.5),(1,0,1,0.5)]
	if nb_color > len(color):
		sys.stdout.write('In the graphe groups will have the same color because there is not enough color defined in the script.\nPlease contact Guillaume MARTIN (guillaume.martin@cirad.fr) for adding colors.\n')
	for n in dico_plot:
		dico_plot[n]['color'] = color[n%len(color)]
		outfile.write('g'+str(n)+'\t'+':'.join(['red='+str(color[n%len(color)][0]), 'green='+str(color[n%len(color)][1]), 'blue='+str(color[n%len(color)][2]), 'alpha='+str(color[n%len(color)][3])])+'\n')
		if int(n/len(color)) > 0:
			sys.stdout.write('In the graphe groups '+str(n)+' will have the same color as the '+str(n%len(color))+' group.\n')
	outfile.close()
	
	# Opening coordinate file
	file = open(MAT)
	# Recording header
	header = file.readline().replace('"','').split()
	# Writting to file
	outfile = open(PREFIX+'_kMean_allele.tab','w')
	outfile.write('\tK-mean_GROUP\t'+'\t'.join(header)+'\n')
	for line in file:
		data = line.split()
		if data:
			if data[0] in grouped:
				outfile.write(data[0]+'\tg'+str(grouped[data[0]])+'\t'+'\t'.join(data[1:])+'\n')
	file.close()
	outfile.close()

def FltrOnMaxGpProp(GPPROPFILE, MAT, GPPROPVALUE, PREFIX):
	
	"""
		Filter grouping matrix based on the maximal grouping proportion calculated in K-mean step
		
		:param GPPROPFILE: The table file with allele grouping information.
		:type GPPROPFILE: str
		:param MAT: The table file containing variant allele encoded plus grouping information.
		:type MAT: str
		:param GPPROPVALUE: The minimal value to keep an allele in the matrix file.
		:type GPPROPVALUE: float
		:param PREFIX: Prefix of the output file.
		:type PREFIX: str
		:return: A filtered table file containing variant allele encoded plus grouping information.
		:rtype: void
	"""
	
	# 1 - Recording variant to keep
	dico_var_to_keep = {}
	file = open(GPPROPFILE)
	header = file.readline().split()
	for line in file:
		data = line.split()
		if data:
			values = list(map(float, data[1:]))
			max_val = max(values)
			if max_val >= GPPROPVALUE:
				dico_var_to_keep[data[0]] = header[values.index(max_val)]
	file.close()
	
	# 2 - Printing variant kept
	outfile = open(PREFIX+'_kMean_allele_filtered_with_'+str(GPPROPVALUE)+'_value.tab', 'w')
	file = open(MAT)
	header = file.readline()
	outfile.write(header)
	for line in file:
		data = line.split()
		if data:
			if data[0] in dico_var_to_keep:
				if data[1] == dico_var_to_keep[data[0]]:
					outfile.write(line)
				else:
					sys.stdout.write('The group '+data[1]+' of variant '+data[0]+' does not correspond to the maximal group '+dico_var_to_keep[data[0]]+'\n')
	outfile.close()
	file.close()

def DrawAllele(VCF, MAT, NAMES, DGROUP, PREFIX, GCOL, FASTA, EXCLCHR):
	
	"""
		Draw circos representation of clustered allele
		
		:param VCF: A vcf file.
		:type VCF: vcf
		:param MAT: The tabulated file containing variant allele encoded plus grouping information.
		:type MAT: str
		:param NAMES: Path to file containing accession to treat.
		:type NAMES: str
		:param DGROUP: A string containing groups to draw.
		:type DGROUP: str
		:param PREFIX: Prefix of the output file.
		:type PREFIX: str
		:param GCOL: A file containing at least a section [color], that define for each group a color for pca drawing (RGB+alpha).
		:type GCOL: str
		:param EXCLCHR: A string containing chromosomes to exclude.
		:type EXCLCHR: str
		:return: Draw circos picture.
		:rtype: void
	"""
	
	
	# Recording chromosomes to exclude
	chr_to_exclude = set()
	if not(EXCLCHR == None):
		chr_to_exclude = set(EXCLCHR.split('='))
		
	# Creation of the Karyotype file
	cumul_size = create_karyotype_file(PREFIX, VCF, FASTA, chr_to_exclude)
	
	#Creating the housekeeping.conf file
	create_houskeeping(PREFIX)
	
	# Creation of the dictionary that will contain information for circos configuration file
	dico_conf = {}
	dico_conf['misc'] = []
	dico_conf['ideogram'] = []
	dico_conf['image'] = []
	dico_conf['colors'] = []
	dico_conf['plots'] = []
	
	dico_conf['misc'].append('karyotype = '+PREFIX+'_Karyotype.tab')
	dico_conf['misc'].append('<<include etc/colors_fonts_patterns.conf>>')
	dico_conf['misc'].append('<<include '+PREFIX+'_housekeeping.conf>>')
	
	dico_conf['ideogram'].append('<spacing>')
	dico_conf['ideogram'].append('default = 0.005r')
	dico_conf['ideogram'].append('</spacing>')
	dico_conf['ideogram'].append('show_label       = yes')
	dico_conf['ideogram'].append('label_font       = bold')
	dico_conf['ideogram'].append('label_with_tag   = yes')
	dico_conf['ideogram'].append('label_size       = 40')
	dico_conf['ideogram'].append('label_parallel   = yes')
	dico_conf['ideogram'].append('label_case       = upper')
	dico_conf['ideogram'].append('label_radius     = dims(ideogram,radius_inner)+ 10p')
	dico_conf['ideogram'].append('radius           = 0.80r')
	dico_conf['ideogram'].append('thickness        = 60p')
	dico_conf['ideogram'].append('fill             = yes')
	dico_conf['ideogram'].append('stroke_color     = dgrey')
	dico_conf['ideogram'].append('stroke_thickness = 2p')
	
	dico_conf['image'].append('<<include etc/image.conf>>')
	dico_conf['image'].append('file*         = '+PREFIX+'_circos.png')
	dico_conf['image'].append('radius*       = 3000p')
	dico_conf['image'].append('angle_offset* = -98.1')
	dico_conf['image'].append('svg* = no')
	
	# For variant sites pritting
	dico_conf['plots'].append('<plot>')
	dico_conf['plots'].append('type             = tile')
	dico_conf['plots'].append('file             = '+PREFIX+'_var_pos.tab')
	dico_conf['plots'].append('r1               = 1.025r')
	dico_conf['plots'].append('r0               = 1.005r')
	dico_conf['plots'].append('orientation      = out')
	dico_conf['plots'].append('margin           = 0u')
	dico_conf['plots'].append('thickness        = 50r')
	dico_conf['plots'].append('padding          = 0')
	dico_conf['plots'].append('stroke_thickness = 1')
	dico_conf['plots'].append('stroke_color     = black')
	dico_conf['plots'].append('color            = black')
	dico_conf['plots'].append('</plot>')
	
	# recording accession names to draw
	acc_to_draw = []
	if NAMES == None:
		sys.stdout.write('No file name was provided in --names argument. All accessions in vcf file (if provided) will be drawn. If no vcf too, all accessions in matrix will be drawn.\n')
		if VCF == None:
			file = open(MAT)
			acc_to_draw = file.readline().replace('K-mean_GROUP','').replace('GROUP','').split()
			file.close()
		else:
			file = open(VCF)
			for line in file:
				data = line.split()
				if data:
					if data[0] == '#CHROM':
						acc_to_draw = data[data.index('FORMAT')+1:]
						break
			file.close()
	else:
		file = open(NAMES)
		for line in file:
			data = line.split()
			if data:
				acc_to_draw.append(data[0])
		file.close()
	
	# recording groups to draw
	group_to_draw = set()
	if DGROUP == None:
		sys.stdout.write('No file groups were provided in --dGroup argument. All groups defined in matrix will be drawn.\n')
	else:
		group_to_draw = set(DGROUP.split('='))
	
	# recording allele grouping
	dico_allele = {}
	file = open(MAT)
	header = file.readline().split()
	if DGROUP == None:
		if 'K-mean_GROUP' in header:
			for line in file:
				data = line.split()
				if data:
					if not(data[0].split(':')[3] == 'A'):
						dico_allele[data[0]] = data[1:][header.index('K-mean_GROUP')]
						group_to_draw.add(data[1:][header.index('K-mean_GROUP')])
		elif 'GROUP' in header:
			for line in file:
				data = line.split()
				if data:
					if not(data[0].split(':')[3] == 'A'):
						dico_allele[data[0]] = data[1:][header.index('GROUP')]
						group_to_draw.add(data[1:][header.index('K-mean_GROUP')])
		else:
			sys.stdout.write('This is embarrassing. The program exited without finishing because no groups were found in the matrix file.\n')
	else:
		if 'K-mean_GROUP' in header:
			for line in file:
				data = line.split()
				if data:
					if not(data[0].split(':')[3] == 'A') and data[1:][header.index('K-mean_GROUP')] in group_to_draw:
						dico_allele[data[0]] = data[1:][header.index('K-mean_GROUP')]
		elif 'GROUP' in header:
			for line in file:
				data = line.split()
				if data:
					if not(data[0].split(':')[3] == 'A') and data[1:][header.index('K-mean_GROUP')] in group_to_draw:
						dico_allele[data[0]] = data[1:][header.index('GROUP')]
		else:
			sys.stdout.write('This is embarrassing. The program exited without finishing because no groups were found in the matrix file.\n')
	file.close()
	group_to_draw = list(group_to_draw)
	
	# Generating tile file for variant sites position
	outfile = open(PREFIX+'_var_pos.tab','w')
	for var in dico_allele:
		pos = var.split(':')
		if not(pos[0] in chr_to_exclude):
			outfile.write('\t'.join([pos[0], pos[1], pos[1]])+'\n')
	outfile.close()
	
	# Generating tiles files of variant allele sites groups and recording group for circos plotting
	i = 1.015
	sys.stdout.write('**********\nColor order in the circos (inner to outer):\n')
	for gp in group_to_draw:
		i += 0.02
		outfile = open(PREFIX+'_var_gp_'+gp+'.tab','w')
		for var in dico_allele:
			if dico_allele[var] == gp:
				pos = var.split(':')
				outfile.write('\t'.join([pos[0], pos[1], pos[1]])+'\n')
		outfile.close()
		# for the circos.conf file
		dico_conf['plots'].append('<plot>')
		dico_conf['plots'].append('type             = tile')
		dico_conf['plots'].append('file             = '+PREFIX+'_var_gp_'+gp+'.tab')
		dico_conf['plots'].append('r1               = '+str(i+0.02)+'r')
		dico_conf['plots'].append('r0               = '+str(i)+'r')
		dico_conf['plots'].append('orientation      = out')
		dico_conf['plots'].append('margin           = 0u')
		dico_conf['plots'].append('thickness        = 48r')
		dico_conf['plots'].append('padding          = 0')
		dico_conf['plots'].append('stroke_thickness = 1')
		dico_conf['plots'].append('stroke_color     = '+gp)
		dico_conf['plots'].append('color            = '+gp)
		dico_conf['plots'].append('</plot>')
		sys.stdout.write(gp+'\n')
	sys.stdout.write('**********\n')
		
	
	# recording color to draw
	dico_color = {}
	if GCOL == None:
		sys.stdout.write('No group color file was provided in --gCol argument. Color will be chosen randomly. If there is more than 7 groups, severals groups will have the same color.\n')
		color = [(0,0,0,0.7),(1,0,0,0.7),(0,1,0,0.7),(0,0,1,0.7),(1,1,0,0.7),(0,1,1,0.7),(1,0,1,0.7)]
		if len(group_to_draw) > len(color):
			sys.stdout.write('In the graphe some groups will have the same color because there is not enough color defined in the script.\nPlease contact Guillaume MARTIN (guillaume.martin@cirad.fr) for adding colors.\n')
		for i in range(len(group_to_draw)):
			dico_color[group_to_draw[i]] = color[i%len(color)]
			if int(i/len(color)) > 0:
				sys.stdout.write('In the graphe groups '+group_to_draw[i]+' will have the same color as the '+str(group_to_draw[i%len(color)])+' group.\n')
	else:
		file = open(GCOL)
		sur_group = False
		sur_color = False
		for line in file:
			data = line.split()
			if data:
				if data[0] == '[group]':
					sur_group = True
					sur_color = False
				elif data[0] == '[color]':
					sur_group = False
					sur_color = True
				elif sur_color:
					if data[0] in group_to_draw:
						color_list = data[1].split(':')
						dic = {}
						for k in color_list:
							dic[k.split('=')[0]] = float(k.split('=')[1])
						if not('red' in dic):
							sys.exit('There is a problem in the file passed to --gCol argument: red color proportion is missing')
						if not('green' in dic):
							sys.exit('There is a problem in the file passed to --gCol argument: green color proportion is missing')
						if not('blue' in dic):
							sys.exit('There is a problem in the file passed to --gCol argument: blue color proportion is missing')
						if not('alpha' in dic):
							sys.exit('There is a problem in the file passed to --gCol argument: alpha color proportion is missing')
						dico_color[data[0]] = (dic['red'], dic['green'], dic['blue'], dic['alpha'])
		# checking all groups have a color
		for n in group_to_draw:
			if not(n in dico_color):
				sys.exit('There is a problem. The group '+n+' have no color. the program exited without finishing')
	
	# recording colors
	for n in group_to_draw:
		dico_conf['colors'].append(n+' = '+str(int(255*dico_color[n][0]))+','+str(int(255*dico_color[n][1]))+','+str(int(255*dico_color[n][2]))+','+str(dico_color[n][3]))
	
	# Creating the tiles files for each accessions
	i = 1
	sys.stdout.write('**********\nAccession order in the circos (outer to inner):\n')
	for acc in acc_to_draw:
		dico_toto = {}
		i -= 0.005
		if i <= 0.01:
			sys.stdout.write('All accessions alleles will not be drawn in the circos because there is to much accession to draw')
			break
		if VCF == None:
			print('toto')
		else:
			# getting ploidy informations
			file = open(VCF)
			for line in file:
				data = line.split()
				if data:
					if data[0] == '#CHROM':
						header = list(data)
						ploidy = False
					elif not(data[0][0] == '#'):
						if ploidy == False:
							genotype = recup_geno(data, header, acc)
							acc_allele = genotype[2:]
							ploidy = len(acc_allele)
							precedent_group = []
							haplo_representative_number = []
							for n in range(ploidy):
								# initialization
								precedent_group.append('')
								haplo_representative_number.append(0)
								outfile = open(PREFIX+'_var_'+acc+'_'+str(n)+'.tab','w')
								outfile.close()
								# for the circos.conf file
								dico_conf['plots'].append('<plot>')
								dico_conf['plots'].append('type             = tile')
								dico_conf['plots'].append('file             = '+PREFIX+'_var_'+acc+'_'+str(n)+'.tab')
								dico_conf['plots'].append('r1               = '+str(i)+'r')
								dico_conf['plots'].append('r0               = '+str(i-0.01)+'r')
								dico_conf['plots'].append('orientation      = out')
								dico_conf['plots'].append('margin           = 0u')
								dico_conf['plots'].append('thickness        = 24r')
								dico_conf['plots'].append('padding          = 0')
								dico_conf['plots'].append('stroke_thickness = 0')
								dico_conf['plots'].append('stroke_color     = white')
								dico_conf['plots'].append('</plot>')
								i -= 0.01
								sys.stdout.write(acc+', allele '+str(n)+'\n')
								if i <= 0.01:
									sys.stdout.write('All accessions alleles will not be drawn in the circos because there is to much accession to draw')
									break
							if i <= 0.01:
								break
							# record for each group number of representatives by haplotypes
							dico_allele_by_gp = {}
							for group in group_to_draw:
								dico_allele_by_gp[group] = []
								for n in range(ploidy):
									dico_allele_by_gp[group].append(0)
						
						# Now working on genotypes
						genotype = recup_geno(data, header, acc)
						acc_allele = genotype[2:]
						
						avail_haplo = []
						for n in range(ploidy):
							avail_haplo.append(n)
						
						# For alleles corresponding to groups just before
						group_found = True
						while group_found:
							group_found = False
							grouped_allele = False
							for allele in acc_allele:
								variant = ':'.join([genotype[0], genotype[1], allele, 'P'])
								if variant in dico_allele:
									group = dico_allele[variant]
									recup_possible_pos = []
									for n in avail_haplo:
										if group == precedent_group[n]:
											recup_possible_pos.append(n)
									# The group is found just before
									if not(len(recup_possible_pos) == 0):
										position = random.choice(recup_possible_pos)
										group_found = True
										grouped_allele = allele
										avail_haplo.remove(position)
										precedent_group[position] = group
										dico_allele_by_gp[group][position] += 1
										haplo_representative_number[position] += 1
										outfile = open(PREFIX+'_var_'+acc+'_'+str(position)+'.tab','a')
										outfile.write('\t'.join([genotype[0], genotype[1], genotype[1], 'color='+group+'\n']))
										outfile.close()
										break
							if grouped_allele:
								acc_allele.remove(grouped_allele)
						
						# For alleles corresponding to groups already found but not immediately before
						list_grouped_allele = []
						for allele in acc_allele:
							variant = ':'.join([genotype[0], genotype[1], allele, 'P'])
							if variant in dico_allele:
								group = dico_allele[variant]
								# recording max group haplotype value
								max_pos_value = 0
								for n in avail_haplo:
									if max_pos_value < dico_allele_by_gp[group][n]:
										max_pos_value = dico_allele_by_gp[group][n]
								# recording haplotypes having the max value if max_pos_value > 0
								if max_pos_value > 0:
									recup_possible_pos = []
									for n in avail_haplo:
										if max_pos_value == dico_allele_by_gp[group][n]:
											recup_possible_pos.append(n)
									# The group has already been found before
									if not(len(recup_possible_pos) == 0):
										position = random.choice(recup_possible_pos)
										list_grouped_allele.append(allele)
										avail_haplo.remove(position)
										precedent_group[position] = group
										dico_allele_by_gp[group][position] += 1
										haplo_representative_number[position] += 1
										outfile = open(PREFIX+'_var_'+acc+'_'+str(position)+'.tab','a')
										outfile.write('\t'.join([genotype[0], genotype[1], genotype[1], 'color='+group+'\n']))
										outfile.close()
						for n in list_grouped_allele:
							acc_allele.remove(n)
						
						# For alleles corresponding which groups were not already found
						list_grouped_allele = []
						for allele in acc_allele:
							variant = ':'.join([genotype[0], genotype[1], allele, 'P'])
							if variant in dico_allele:
								group = dico_allele[variant]
								# recording max haplotype value
								min_haplo_value = False
								for n in avail_haplo:
									if not(min_haplo_value):
										min_haplo_value = haplo_representative_number[n]
									elif min_haplo_value > haplo_representative_number[n]:
										min_haplo_value = haplo_representative_number[n]
								# recording haplotypes having the minimal value
								recup_possible_pos = []
								for n in avail_haplo:
									if min_haplo_value == haplo_representative_number[n]:
										recup_possible_pos.append(n)
								if len(recup_possible_pos) == 0:
									sys.exit('There is a bug when trying to group alleles on haplotypes')
								else:
									position = random.choice(recup_possible_pos)
									list_grouped_allele.append(allele)
									avail_haplo.remove(position)
									precedent_group[position] = group
									dico_allele_by_gp[group][position] += 1
									haplo_representative_number[position] += 1
									outfile = open(PREFIX+'_var_'+acc+'_'+str(position)+'.tab','a')
									outfile.write('\t'.join([genotype[0], genotype[1], genotype[1], 'color='+group+'\n']))
									outfile.close()
						for n in list_grouped_allele:
							acc_allele.remove(n)
			file.close()
	sys.stdout.write('**********\n')
	
	# Generating the configuration file
	outfile = open(PREFIX+'_circos.conf','w')
	for n in dico_conf:
		if n == 'misc':
			outfile.write('\n'.join(dico_conf[n]))
		else:
			outfile.write('\n\n<'+n+'>\n\n')
			outfile.write('\n'.join(dico_conf[n]))
			outfile.write('\n\n</'+n+'>\n\n')
	outfile.close()
	
	# Drawing circos
	os.system('circos --conf '+PREFIX+'_circos.conf')

def	create_karyotype_file(PREFIX, VCF, FASTA, EXCLCHR):
	
	"""
		Create a karyotype file that will be used for circos
		
		:param PREFIX: Prefix of the output file.
		:type PREFIX: str
		:param VCF: A vcf file.
		:type VCF: vcf
		:param FASTA: Fasta file containing reference sequences.
		:type FASTA: fasta
		:param EXCLCHR: A dictionary containing chromosomes to exclude.
		:type EXCLCHR: dictionary
		:return: Cumuled chromosome size
		:rtype: int
	"""
	
	# Recording information
	if FASTA == None:
		dico_chr = {}
		file = open(VCF)
		for line in file:
			data = line.split()
			if data:
				if data[0][0:9] == '##contig=':
					info = get_chr_size(data)
					dico_chr[info[0]] = info[1]
				elif data[0] == '#CHROM':
					break
	else:
		dico_chr = get_chr_size_from_fasta(FASTA)
	
	# calculating chromosomes cumulated size
	cumul_size = 0
	for n in sorted(dico_chr.keys()):
		if not(n in EXCLCHR):
			cumul_size += dico_chr[n]
	
	# Creating Karyotype file
	outfile = open(PREFIX+'_Karyotype.tab','w')
	outfile.write('chr - legend legend 0 '+str(int(cumul_size*5/100))+' white\n')
	for n in sorted(dico_chr.keys()):
		if not(n in EXCLCHR):
			outfile.write('chr - '+n+' '+n+' 0 '+str(dico_chr[n])+' white\n')
	outfile.close()
	
	return cumul_size

def create_houskeeping(PREFIX):
	
	"""
		Create the housekepping.conf file (needed if many point are drawn)
		
		:param PREFIX: Prefix of the output file.
		:type PREFIX: str
		:return: the housekeeping.conf file
		:rtype: void
	"""
	
	outfile = open(PREFIX+'_housekeeping.conf','w')
	outfile.write('anglestep       = 0.5\n')
	outfile.write('minslicestep    = 10\n')
	outfile.write('beziersamples   = 40\n')
	outfile.write('debug           = no\n')
	outfile.write('warnings        = no\n')
	outfile.write('imagemap        = no\n')
	outfile.write('paranoid        = yes\n')
	outfile.write('units_ok        = bupr\n')
	outfile.write('units_nounit    = n\n')
	outfile.write('file_delim = \s\n')
	outfile.write('file_delim_collapse = yes\n')
	outfile.write('list_record_delim = \s*[;,]\s*\n')
	outfile.write('list_field_delim  = \s*[:=]\s*\n')
	outfile.write('options_record_delim = [,;]\n')
	outfile.write('options_field_delim  = =\n')
	outfile.write('skip_missing_expression_vars = no\n')
	outfile.write('legacy_underline_expression_syntax = no\n')
	outfile.write('svg_font_scale = 1.3\n')
	outfile.write('sup_baseline_shift = 40\n')
	outfile.write('sub_baseline_shift = -40\n')
	outfile.write('sup_fontsize = 90\n')
	outfile.write('sub_fontsize = 90\n')
	outfile.write('default_font   = default\n')
	outfile.write('default_font_name  = Arial\n')
	outfile.write('default_font_color = black\n')
	outfile.write('default_color  = black\n')
	outfile.write('<guides>\n')
	outfile.write('thickness      = 1\n')
	outfile.write('size           = 5\n')
	outfile.write('type           = outline\n')
	outfile.write('<object>\n')
	outfile.write('all            = no\n')
	outfile.write('ideogram       = no\n')
	outfile.write('ideogram_label = no\n')
	outfile.write('</object>\n')
	outfile.write('<color>\n')
	outfile.write('default = lblue\n')
	outfile.write('text    = red\n')
	outfile.write('</color>\n')
	outfile.write('</guides>\n')
	outfile.write('debug_group = summary,output\n')
	outfile.write('debug_auto_timer_report = 30\n')
	outfile.write('debug_word_separator = " "\n')
	outfile.write('debug_undef_text     = _undef_\n')
	outfile.write('debug_empty_text     = _emptylist_\n')
	outfile.write('debug_validate       = yes\n')
	outfile.write('debug_output_tidy    = no\n')
	outfile.write('text_pixel_subsampling = 1\n')
	outfile.write('text_snuggle_method    = array\n')
	outfile.write('restrict_parameter_names = no\n')
	outfile.write('case_sensitive_parameter_names = no\n')
	outfile.write('calculate_track_statistics = yes\n')
	outfile.write('color_cache_static = yes\n')
	outfile.write('color_cache_file   = circos.colorlist\n')
	outfile.write('color_lists_use    = yes\n')
	outfile.write('memoize = yes\n')
	outfile.write('quit_on_dump = yes\n')
	outfile.write('offsets = 0,0\n')
	outfile.write('max_ticks            = 5000\n')
	outfile.write('max_ideograms        = 200\n')
	outfile.write('max_links            = 25000\n')
	outfile.write('max_points_per_track = 1000000\n')
	outfile.write('undefined_ideogram = skip\n')
	outfile.write('relative_scale_iterations = 10\n')
	outfile.write('relative_scale_spacing    = mode\n')
	outfile.write('data_out_of_range = trim,warn #\n')
	outfile.write('track_defaults = etc/tracks\n')
	outfile.write('round_brush_use           = yes\n')
	outfile.write('round_brush_min_thickness = 5\n')
	outfile.write('anti_aliasing = yes\n')
	outfile.write('housekeeping = yes\n')
	outfile.write('auto_eval = no\n')
	outfile.close()

def comb(liste, rest, nb, liste_total):
	"""Create all possible combinations of n elements contained in a list"""
	if nb == 0:
		liste_total.append(liste)
	else:
		for i in range(len(rest)):
			sub_list = list(liste)
			sub_list.append(rest[i])
			comb(sub_list, rest[i:], nb-1, liste_total)

def comb_ss_repet(liste, rest, liste_total, nb):
	"""Create all possible combinations of n elements contained in a list without repetition"""
	if nb == 0:
		liste_total.append(liste)
	else:
		for i in range(len(rest)):
			sub_list = list(liste)
			sub_list.append(rest[i])
			comb_ss_repet(sub_list, rest, liste_total, nb-1)

def MergeVcf(VCF, VCF2, COMP1, PREFIX):
	
	"""
		Add variant calling of an individual contained in VCF2 to VCF1
		
		:param VCF: A VCF file.
		:type VCF: str
		:param VCF2: A second VCF file.
		:type VCF2: str
		:param COMP1: Accession name to add to VCF1.
		:type COMP1: str
		:param PREFIX: Prefix of the output file.
		:type PREFIX: str
		:return: A vcf file
		:rtype: void
	"""
	
	# loading informations to add to vcf1
	dico_acc = {}
	file = open(VCF2)
	for line in file:
		data = line.split()
		if data:
			# recording header
			if data[0] == '#CHROM':
				header = list(data)
			# not working on headers of the file
			elif data[0][0] != '#':
				format = data[header.index('FORMAT')].split(':')
				allele_cov = data[header.index(COMP1)].split(':')[format.index('AD')].split(',')
				variants_all = [data[header.index('REF')]]
				variants_all = variants_all + data[header.index('ALT')].split(',')
				genotype = recup_geno(data, header, COMP1)
				chrom = genotype[0]
				if not(chrom in dico_acc):
					dico_acc[chrom] = {}
				pos = genotype[1]
				dico_acc[chrom][pos] = {}
				for i in range(len(variants_all)):
					dico_acc[chrom][pos][variants_all[i]] = allele_cov[i]
				dico_acc[chrom][pos]['geno'] = genotype[2:]
				dico_acc[chrom][pos]['value'] = data[header.index(COMP1)].split(':')[2:]
	ploidy = len(genotype[2:])
	file.close()
	
	# printing the new vcf file
	found = 0
	not_found = 0
	outfile = open(PREFIX+'_merged.vcf','w')
	file = open(VCF)
	for line in file:
		data = line.split()
		if data:
			# recording header
			if data[0] == '#CHROM':
				header = list(data)
				outfile.write('\t'.join(header)+'\t'+COMP1+'\n')
			# not working on headers of the file
			elif data[0][0] == '#':
				outfile.write(line)
			else:
				variants_all = [data[header.index('REF')]]
				variants_all = variants_all + data[header.index('ALT')].split(',')
				alternate_allele = data[header.index('ALT')].split(',')
				chrom = data[header.index('#CHROM')]
				pos = data[header.index('POS')]
				if not(chrom in dico_acc):
					variant_acc = []
					for n in data[header.index('FORMAT')].split(':'):
						if n == 'GT':
							mot = []
							for i in range(ploidy):
								mot.append('.')
							variant_acc.append('/'.join(mot))
						elif n == 'AD':
							mot = []
							for i in range(variants_all):
								mot.append('.')
							variant_acc.append(','.join(mot))
						elif n == 'PL':
							mot = []
							for i in range(int(math.factorial(len(variants_all)+ploidy-1) / (math.factorial(ploidy)*(len(variants_all)-1)))):
								mot.append('.')
							variant_acc.append(','.join(mot))
						else:
							variant_acc.append('.')
					data.append(':'.join(variant_acc))
					outfile.write('\t'.join(data)+'\n')
					not_found += 1
				elif not(pos in dico_acc[chrom]):
					variant_acc = []
					for n in data[header.index('FORMAT')].split(':'):
						if n == 'GT':
							mot = []
							for i in range(ploidy):
								mot.append('.')
							variant_acc.append('/'.join(mot))
						elif n == 'AD':
							mot = []
							for i in range(len(variants_all)):
								mot.append('.')
							variant_acc.append(','.join(mot))
						elif n == 'PL':
							mot = []
							for i in range(int(math.factorial(len(variants_all)+ploidy-1) / (math.factorial(ploidy)*(len(variants_all)-1)))):
								mot.append('.')
							variant_acc.append(','.join(mot))
						else:
							variant_acc.append('.')
					data.append(':'.join(variant_acc))
					outfile.write('\t'.join(data)+'\n')
					not_found += 1
				else:
					if 'NA' in dico_acc[chrom][pos]['geno']:
						final_variant = dico_acc[chrom][pos]['value']
						mot = []
						for allele in variants_all:
							if allele in dico_acc[chrom][pos]:
								mot.append(dico_acc[chrom][pos][allele])
							else:
								mot.append('0')
						final_variant.insert(0,','.join(mot))
						mot = []
						for i in range(ploidy):
							mot.append('.')
						final_variant.insert(0,'/'.join(mot))
						not_found += 1
					else:
						for allele in dico_acc[chrom][pos]['geno']:
							if not(allele in variants_all):
								variants_all.append(allele)
								alternate_allele.append(allele)
						final_variant = dico_acc[chrom][pos]['value']
						mot = []
						for allele in variants_all:
							if allele in dico_acc[chrom][pos]:
								mot.append(dico_acc[chrom][pos][allele])
							else:
								mot.append('0')
						final_variant.insert(0,','.join(mot))
						mot = []
						for acc_allele in dico_acc[chrom][pos]['geno']:
							mot.append(str(variants_all.index(acc_allele)))
						final_variant.insert(0,'/'.join(mot))
						found += 1
					data[header.index('ALT')] = ','.join(alternate_allele)
					data.append(':'.join(final_variant))
					outfile.write('\t'.join(data)+'\n')
	outfile.close()
	file.close()
	
	sys.stdout.write('Number of sites with genotype calling: '+str(found)+'\n')
	sys.stdout.write('Number of sites with no genotype calling: '+str(not_found)+'\n')

def get_color(GCOL, GROUP_TO_DRAW, dico_color):
	if GCOL == None:
		sys.stdout.write('No group color file was provided in --gCol argument. Color will be chosen randomly. If there is more than 7 groups, severals groups will have the same color.\n')
		color = [(0,0,0,0.7),(1,0,0,0.7),(0,1,0,0.7),(0,0,1,0.7),(1,1,0,0.7),(0,1,1,0.7),(1,0,1,0.7)]
		if len(GROUP_TO_DRAW) > len(color):
			sys.stdout.write('In the graphe some groups will have the same color because there is not enough color defined in the script.\nPlease contact Guillaume MARTIN (guillaume.martin@cirad.fr) for adding colors.\n')
		for i in range(len(GROUP_TO_DRAW)):
			dico_color[GROUP_TO_DRAW[i]] = color[i%len(color)]
			if int(i/len(color)) > 0:
				sys.stdout.write('In the graphe groups '+GROUP_TO_DRAW[i]+' will have the same color as the '+str(GROUP_TO_DRAW[i%len(color)])+' group.\n')
	else:
		file = open(GCOL)
		sur_group = False
		sur_color = False
		for line in file:
			data = line.split()
			if data:
				if data[0] == '[group]':
					sur_group = True
					sur_color = False
				elif data[0] == '[color]':
					sur_group = False
					sur_color = True
				elif sur_color:
					if data[0] in GROUP_TO_DRAW:
						color_list = data[1].split(':')
						dic = {}
						for k in color_list:
							dic[k.split('=')[0]] = float(k.split('=')[1])
						if not('red' in dic):
							sys.exit('There is a problem in the file passed to --gCol argument: red color proportion is missing')
						if not('green' in dic):
							sys.exit('There is a problem in the file passed to --gCol argument: green color proportion is missing')
						if not('blue' in dic):
							sys.exit('There is a problem in the file passed to --gCol argument: blue color proportion is missing')
						if not('alpha' in dic):
							sys.exit('There is a problem in the file passed to --gCol argument: alpha color proportion is missing')
						dico_color[data[0]] = (dic['red'], dic['green'], dic['blue'], dic['alpha'])
		# checking all groups have a color
		for n in GROUP_TO_DRAW:
			if not(n in dico_color):
				sys.exit('There is a problem. The group '+n+' have no color. the program exited without finishing')

def	get_allele_group_prop(VCF, MAT, PREFIX, NAMES, DGROUP, EXCLCHR):
	
	"""
		Identify genome structure
		
		:param VCF: A vcf file.
		:type VCF: vcf
		:param MAT: The tabulated file containing variant allele encoded plus grouping information.
		:type MAT: str
		:param NAMES: Path to file containing accession to treat.
		:type NAMES: str
		:param DGROUP: A string containing groups to draw.
		:type DGROUP: str
		:param PREFIX: Prefix of the output file.
		:type PREFIX: str
		:param EXCLCHR: A string containing chromosomes to exclude.
		:type EXCLCHR: str
		:return: A file containing the allele group number for each accessions.
		:rtype: void
	"""
	
	# Recording chromosomes to exclude
	chr_to_exclude = set()
	if not(EXCLCHR == None):
		chr_to_exclude = set(EXCLCHR.split(':'))
	
	# recording accession names to work with
	acc_to_draw = []
	if NAMES == None:
		sys.stdout.write('No file name was provided in --names argument. All accessions in vcf file (if provided) will be treated. If no vcf too, all accessions in matrix will be treated.\n')
		if VCF == None:
			file = open(MAT)
			acc_to_draw = file.readline().replace('K-mean_GROUP','').replace('GROUP','').split()
			file.close()
		else:
			file = open(VCF)
			for line in file:
				data = line.split()
				if data:
					if data[0] == '#CHROM':
						acc_to_draw = data[data.index('FORMAT')+1:]
						break
			file.close()
	else:
		file = open(NAMES)
		for line in file:
			data = line.split()
			if data:
				acc_to_draw.append(data[0])
		file.close()
	
	# recording groups to work with
	group_to_draw = set()
	if DGROUP == None:
		sys.stdout.write('No file groups were provided in --dGroup argument. All groups defined in matrix will be treated.\n')
	else:
		group_to_draw = set(DGROUP.split(':'))
	
	
	sys.stdout.write("Recording alleles in grouping informations...\n")
	# recording allele grouping
	dico_allele = {}
	dico_sites = set()
	file = open(MAT)
	header = file.readline().split()
	if DGROUP == None:
		if 'K-mean_GROUP' in header:
			for line in file:
				data = line.split()
				if data:
					variant = data[0].split(':')
					if not(variant[0] in chr_to_exclude):
						if not(variant[3] == 'A'):
							dico_allele[data[0]] = data[1:][header.index('K-mean_GROUP')]
							group_to_draw.add(data[1:][header.index('K-mean_GROUP')])
							dico_sites.add(':'.join(variant[0:2]))
		elif 'GROUP' in header:
			for line in file:
				data = line.split()
				if data:
					variant = data[0].split(':')
					if not(variant[0] in chr_to_exclude):
						if not(variant[3] == 'A'):
							dico_allele[data[0]] = data[1:][header.index('GROUP')]
							group_to_draw.add(data[1:][header.index('GROUP')])
							dico_sites.add(':'.join(variant[0:2]))
		else:
			sys.stdout.write('This is embarrassing. The program exited without finishing because no groups were found in the matrix file.\n')
	else:
		if 'K-mean_GROUP' in header:
			for line in file:
				data = line.split()
				if data:
					variant = data[0].split(':')
					if not(variant[0] in chr_to_exclude):
						if not(variant[3] == 'A') and data[1:][header.index('K-mean_GROUP')] in group_to_draw:
							dico_allele[data[0]] = data[1:][header.index('K-mean_GROUP')]
							dico_sites.add(':'.join(variant[0:2]))
		elif 'GROUP' in header:
			for line in file:
				data = line.split()
				if data:
					variant = data[0].split(':')
					if not(variant[0] in chr_to_exclude):
						if not(variant[3] == 'A') and data[1:][header.index('GROUP')] in group_to_draw:
							dico_allele[data[0]] = data[1:][header.index('GROUP')]
							dico_sites.add(':'.join(variant[0:2]))
		else:
			sys.stdout.write('This is embarrassing. The program exited without finishing because no groups were found in the matrix file.\n')
	file.close()
	group_to_draw = list(group_to_draw)
	
	sys.stdout.write("Done\n")
	################################
	# At this point we have:
	#	group_to_draw 	= a list of groups that will be used for analysis
	#	dico_allele 	= a dictionary with key --> variant allele; value --> group
	#	dico_sites		= a set with variant sites to work with
	#	acc_to_draw		= a list of accession to work with
	#
	######################################################################
	#    Now we are recording informations contained in the vcf files    #
	######################################################################
	
	
	# a dictionary with allele group number observed
	dico_acc_freq = {}
	for n in acc_to_draw:
		dico_acc_freq[n] = {}
		for j in group_to_draw:
			dico_acc_freq[n][j] = 0
	
	sys.stdout.write("Recording VCF informations...\n")
	file = open(VCF)
	for line in file:
		data = line.split()
		if data:
			if data[0] == '#CHROM':
				# recording accession name
				header = list(data)
			elif not(data[0][0] == '#'):
				# getting variant site position
				chrom = data[header.index('#CHROM')]
				pos = data[header.index('POS')]
				if ':'.join([chrom, pos]) in dico_sites:
					for acc in acc_to_draw:
						genotype = recup_geno(data, header, acc)
						for allele in genotype[2:]:
							variant = ':'.join([chrom,pos,allele,'P'])
							if variant in dico_allele:
								group = dico_allele[variant]
								if group in dico_acc_freq[acc]:
									dico_acc_freq[acc][group] += 1
	
	outfile = open(PREFIX+'_gp_prop.tab','w')
	outfile.write('\t'.join(['group-accessions']+acc_to_draw)+'\n')
	for gp in group_to_draw:
		mot = str(gp)
		for acc in acc_to_draw:
			mot = mot+'\t'+str(dico_acc_freq[acc][gp])
		outfile.write(mot+'\n')
	outfile.close()


def __main__():
	#Parse Command Line
	parser = optparse.OptionParser(usage="python %prog [options]\n\nProgram designed by Guillaume MARTIN : guillaume.martin@cirad.fr")
	# Wrapper options. 
	parser.add_option( '',	'--vcf',			dest='vcf',				default=None,			help='The vcf file. [Default: %default]')
	parser.add_option( '',	'--vcf2',			dest='vcf2',			default=None,			help='A second vcf file. If COMPARE is passed in --type argument: can be used to compare two variant calling in two vcf file. [Default: %default]')
	parser.add_option( '',	'--names',			dest='names',			default=None,			help='A one column file containing accession names to treat. [Default: %default]')
	parser.add_option( '',	'--outgroup',		dest='outgroup',		default=None,			help='A one column file containing accession names that will not be used for filtering but will remain in the output file. [Default: %default]')
	parser.add_option( '',	'--type',			dest='type',			default=None,			help='Type of treatment to perform: STAT, GENE_STAT, FILTER, COMPARE, ADD_REF, AL_IDENTITY, FACTORIAL, RANDOM_SUB_SET, VISUALIZE_VAR_3D, VISUALIZE_VAR_2D, SNP_CLUST-Kmean, SNP_CLUST-MeanShift, FILTER_ON_MAX_GP_PROP, ALLELIC_STRUCT, GENOME_BLOCS, MERGE_VCF, GET_GENOTYPE, ALL_PROP. [Default: %default]')
	parser.add_option( '',	'--fasta',			dest='fasta',			default=None,			help='A fasta file containing reference sequence. Not needed if standard vcf file (with sequence length). [Default: %default]')
	parser.add_option( '',	'--gff3',			dest='gff3',			default=None,			help='A gff3 file containing gene annotation. [Default: %default]')
	parser.add_option( '',	'--RmType',			dest='RmType',			default=None,			help='Variant status to filter out (several values can be passed in this case they should be separated by :). Values: PASS, DP_FILTER, QD_FILTER, SnpCluster, INDELS, SNP, AUTAPO [Default: %default]')
	parser.add_option( '',	'--RmAlAlt',		dest='RmAlAlt',			default=None,			help='Number of alleles at the site to remove the variant site (several values can be passed and should be sepatated by :). Values: 1,2,3,...,n [Default: %default]')
	parser.add_option( '',	'--MinCov',			dest='MinCov',			default='10',			help='Minimal coverage by accession to keep genotype calling (integer). If the value is lower, genotype will be converted to unknown for the concerned accession. [Default: %default]')
	parser.add_option( '',	'--MinAl',			dest='MinAl',			default='3',			help='Minimal allele coverage by accession to keep genotype calling (integer). If the value is lower for at least one allele, genotype will be converted to unknown for the concerned accession. [Default: %default]')
	parser.add_option( '',	'--nMiss',			dest='nMiss',			default='0',			help='Maximal number of missing genotype in a line to keep the line (integer). [Default: %default]')
	parser.add_option( '',	'--nRand',			dest='nRand',			default='1000',			help='Number of variant site to get randomly from the vcf file (integer). [Default: %default]')
	parser.add_option( '',	'--mulType',		dest='mulType',			default='coa',			help='Multivariate analysis type. Possible values: coa, pca and pca_normed [Default: %default]')
	parser.add_option( '',	'--nAxes',			dest='nAxes',			default='4',			help='Axes number to keep after PCA analysis (integer). [Default: %default]')
	parser.add_option( '',	'--dAxes',			dest='dAxes',			default=None,			help='Axes to draw/use in 3d plot/kmean clustering. Axis should be separated by ":". Three values are required for 3d plot. [Default: %default]')
	parser.add_option( '',	'--group',			dest='group',			default=None,			help='A file containing two sections: A section[group] with in col 1 accession name ; col 2 group (UN for unknown group, accessions with unknown group could be ommited). A section [color], that define for each group a color for pca drawing (in RGB+alpha percentage, ex: red=1:green=0:blue=0:alpha=0.1). [Default: %default]')
	parser.add_option( '',	'--dGroup',			dest='dGroup',			default=None,			help='Groups ids to draw. Groups names should be separated by ":". Groups names will be searched in the file provided in --mat argument in the "K-mean_GROUP" column, and if such column is not present, they will be searched in "GROUP" columns. If none of these columns are found the program will exit without finishing. [Default: %default]')
	parser.add_option( '',	'--gCol',			dest='gCol',			default=None,			help='A file containing color attributed to group. This is a table file with column 1 : group name, column 2 : rgb color (ex: group1 red=1:green=0:blue=0) [Default: %default]')
	parser.add_option( '',	'--prefix',			dest='prefix',			default='WorkOnVcf', 	help='The prefix for output files. [Default: %default]')
	parser.add_option( '',	'--comp1',			dest='comp1',			default=None, 			help='First accession name to compare. [Default: %default]')
	parser.add_option( '',	'--comp2',			dest='comp2',			default=None, 			help='Second accession name to compare. If not provided, the same accession name will be searched to compare calling vcf passed in --vcf and --vcf2 arguments. If --comp2 and --vcf2 arguments are passed --comp2 accession name will be searched in vcf2. If --vcf2 argument is not passed, --comp2 accession name will be searched in vcf. [Default: %default]')
	parser.add_option( '',	'--ref_cov',		dest='ref_cov',			default='1000', 		help='Value to put to AD and DP flags if "ADD_REF" is passed to --type. [Default: %default]')
	parser.add_option( '',	'--VarCoord',		dest='VarCoord',		default=None, 			help='The tabulated file of variables coordinates in new axis. (The --prefix + _variables_coordinates.tab file generated when runing this script with PCA type) [Default: %default]')
	parser.add_option( '',	'--mat',			dest='mat',				default=None, 			help='The tabulated file containing variant allele encoded. (The --prefix + _matrix_4_PCA.tab file generated when runing this script with PCA type) [Default: %default]')
	parser.add_option( '',	'--nGroup',			dest='nGroup',			default='2', 			help='Group number for the k-mean algorithm that will cluster variables (alleles) based on their coordinates. [Default: %default]')
	parser.add_option( '',	'--quantile',		dest='quantile',		default='0.2', 			help='The quantile value to estimate de bandwidth parameters used in the MeanShift. Value should be in [0:1]. [Default: %default]')
	parser.add_option( '',	'--ExclChr',		dest='ExclChr',			default=None, 			help='Chromosome names to exclude from analysis. Each chromosomes should be separated by ":". [Default: %default]')
	parser.add_option( '',	'--thread',			dest='thread',			default='1', 			help='Number of processor available. [Default: %default]')
	parser.add_option( '',	'--iter',			dest='iter',			default='100', 			help='Parallele k-mean clustering attempt to try. [Default: %default]')
	parser.add_option( '',	'--gpPropFile',		dest='gpPropFile',		default=None, 			help='The group file proportion (*_kMean_gp_prop.tab file generated with SNP_CLUST). [Default: %default]')
	parser.add_option( '',	'--gpPropValue',	dest='gpPropValue',		default='0.95',			help='The minimal value the max group proportion file to keep the dot (float between 0 and 1). [Default: %default]')
	parser.add_option( '',	'--win',			dest='win',				default='25', 			help='Half window size around a variant site to evaluate the structure at the site. [Default: %default]')
	parser.add_option( '',	'--MeanShiftAll',	dest='MeanShiftAll',	default='y', 			help='Cluster all point in the MeanShift. Possible values, "y" or "n" [Default: %default]')
	parser.add_option( '',	'--bandwidth',		dest='bandwidth',		default=None, 			help='Bandwidth value used for mean shift. If filled, the --quantile parameter is ignored. [Default: %default]')
	parser.add_option( '',	'--AP',				dest='AP',				default='n', 			help='Cluster absent (A) and present (P) lines. Possible values, "y" or "n" [Default: %default]')

	(options, args) = parser.parse_args()
	
	# draw_color_map()
	# sys.exit()
	
	if options.type == None:
		sys.exit('Please provide an analysis type to --type argument')
	
	# Extract a random sub-set of the VCF
	if options.type == 'RANDOM_SUB_SET':
		sys.stdout.write('Associated parameters:\n\t--vcf\n\t--nRand\n\t--prefix\n')
		if options.vcf == None:
			sys.exit('Please provide a vcf file to --vcf argument')
		random_sub_set(options.vcf, int(options.nRand), options.prefix)
	
	# Calculating statistics
	if options.type == 'STAT':
		sys.stdout.write('Associated parameters:\n\t--vcf\n\t--prefix\n\t--names\n\t--gff3\n')
		if options.vcf == None:
			sys.exit('Please provide a vcf file to --vcf argument')
		if options.names == None:
			sys.exit('Please provide a name file to --names argument')
		stat_on_vcf(options.vcf, options.prefix, options.names, options.gff3)
	
	# Filtering vcf file
	if options.type == 'FILTER':
		sys.stdout.write('Associated parameters:\n\t--vcf\n\t--names\n\t--outgroup\n\t--prefix\n\t--RmType\n\t--MinCov\n\t--MinAl\n\t--nMiss\n\t--RmAlAlt\n')
		if options.vcf == None:
			sys.exit('Please provide a vcf file to --vcf argument')
		if options.names == None:
			sys.exit('Please provide a name file to --names argument')
		filter_vcf(options.vcf, options.names, options.outgroup, options.prefix, options.RmType, int(options.MinCov), int(options.MinAl), int(options.nMiss), options.RmAlAlt)
	
	# Comparing two variant calling
	if options.type == 'COMPARE':
		sys.stdout.write('Associated parameters:\n\t--vcf\n\t--vcf2\n\t--comp1\n\t--comp2\n')
		if options.vcf == None:
			sys.exit('Please provide a vcf file to --vcf argument')
		if options.comp1 == None:
			sys.exit('Please provide an accession name to --comp1 argument')
		compare2calling(options.vcf, options.vcf2, options.comp1, options.comp2)
	
	# Add reference to calling
	if options.type == 'ADD_REF':
		sys.stdout.write('Associated parameters:\n\t--vcf\n\t--prefix\n\t--ref_cov\n')
		if options.vcf == None:
			sys.exit('Please provide a vcf file to --vcf argument')
		AddRefToCall(options.vcf, options.prefix, int(options.ref_cov))
	
	# Callculate genotype identity
	if options.type == 'AL_IDENTITY':
		sys.stdout.write('Associated parameters:\n\t--vcf\n\t--prefix\n')
		if options.vcf == None:
			sys.exit('Please provide a vcf file to --vcf argument')
		sys.stdout.write('Warning, the statistic used here has been designed to work with data with different ploidy level and does not correpondond to an identity as you may expect\n')
		CalcAllelicIdent(options.vcf, options.prefix)
	
	# Do ACP
	if options.type == 'FACTORIAL':
		sys.stdout.write('Associated parameters:\n\t--vcf\n\t--names\n\t--prefix\n\t--group\n\t--nAxes\n\t--mulType\n')
		sys.stdout.write('Optional and dependant parameters:\n\t--dGroup: If passed, all alleles belonging to groups passed to this option will be removed.\n\t--mat: Matrix of grouped alleles (with either a GROUP or a K-mean_GROUP column).'
		' If a K-mean_GROUP column is found, the filter will be performed on this column, else it will be performed on the GROUP one\n')
		if options.vcf == None:
			sys.exit('Please provide a vcf file to --vcf argument')
		FormatForPCA(options.vcf, options.names, options.prefix, options.group, int(options.nAxes), options.mulType, options.dGroup, options.mat)
		os.system('R CMD BATCH '+options.prefix+'_multivariate.R')
	
	# Draw 3d plot
	if options.type == 'VISUALIZE_VAR_3D':
		sys.stdout.write('Associated parameters:\n\t--VarCoord\n\t--dAxes\n\t--mat\n\t--group\n\t--dGroup\n')
		if options.VarCoord == None:
			sys.exit('Please provide a file name to --VarCoord argument')
		Draw3dPlot(options.VarCoord, options.dAxes, options.mat, options.group, options.dGroup)
	
	# Draw 2d plot
	if options.type == 'VISUALIZE_VAR_2D':
		sys.stdout.write('Associated parameters:\n\t--VarCoord\n\t--dAxes\n\t--mat\n\t--group\n\t--dGroup\n\t--prefix\n')
		if options.VarCoord == None:
			sys.exit('Please provide a file name to --VarCoord argument')
		Draw2dPlot(options.VarCoord, options.dAxes, options.mat, options.group, options.dGroup, options.prefix)
	
	# Cluster SNP
	if options.type == 'SNP_CLUST-Kmean':
		sys.stdout.write('Associated parameters:\n\t--VarCoord\n\t--mat\n\t--dAxes\n\t--nGroup\n\t--AP\n\t--iter\n\t--thread\n\t--prefix\n')
		if options.dAxes == None:
			sys.exit('Please provide value to --dAxes argument')
		if options.mat == None:
			sys.exit('Please provide a matrix file to --mat argument')
		# kMean_clust(options.VarCoord, options.dAxes, int(options.nGroup), options.mat, options.prefix, int(options.thread), int(options.iter))
		NewKmean(options.VarCoord, options.mat, options.dAxes, int(options.nGroup), options.AP, int(options.iter), int(options.thread), options.prefix)
	
	# Cluster SNP
	if options.type == 'SNP_CLUST-MeanShift':
		sys.stdout.write('Associated parameters:\n\t--VarCoord\n\t--mat\n\t--dAxes\n\t--quantile\n\t--AP\n\t--iter\n\t--thread\n\t--MeanShiftAll\n\t--bandwidth\n\t--prefix\n')
		if options.dAxes == None:
			sys.exit('Please provide value to --dAxes argument')
		if options.mat == None:
			sys.exit('Please provide a matrix file to --mat argument')
		NewMeanShift(options.VarCoord, options.mat, options.dAxes, float(options.quantile), options.AP, int(options.iter), int(options.thread), options.MeanShiftAll, options.bandwidth, options.prefix)
	
	# Filtering on max grouping proportion
	if options.type == 'FILTER_ON_MAX_GP_PROP':
		sys.stdout.write('Associated parameters:\n\t--gpPropFile\n\t--mat\n\t--gpPropValue\n\t--prefix\n')
		if options.gpPropFile == None:
			sys.exit('Please provide a file to --gpPropFile argument')
		if options.mat == None:
			sys.exit('Please provide a matrix file to --mat argument')
		if options.gpPropValue == None:
			sys.exit('Please provide value to --gpPropValue argument')
		FltrOnMaxGpProp(options.gpPropFile, options.mat, float(options.gpPropValue), options.prefix)
	
	# Merging vcf based on names
	if options.type == 'MERGE_VCF':
		sys.stdout.write('Associated parameters:\n\t--vcf\n\t--vcf2\n\t--comp1\n\t--prefix\n')
		if options.vcf == None:
			sys.exit('Please provide a vcf file to --vcf argument')
		if options.vcf2 == None:
			sys.exit('Please provide a vcf file to --vcf2 argument')
		if options.comp1 == None:
			sys.exit('Please provide an accession name to --comp1 argument')
		MergeVcf(options.vcf, options.vcf2, options.comp1, options.prefix)
	
	# Get genotype from VCF
	if options.type == 'GET_GENOTYPE':
		sys.stdout.write('Associated parameters:\n\t--vcf\n\t--names\n\t--prefix\n')
		if options.vcf == None:
			sys.exit('Please provide a vcf file to --vcf argument')
		get_genotype(options.vcf, options.names, options.prefix)
	
	# Calculate allele group proportions
	if options.type == 'ALL_PROP':
		sys.stdout.write('Associated parameters:\n\t--vcf\n\t--mat\n\t--prefix\n\t--names\n\t--dGroup\n\t--ExclChr\n')
		if options.vcf == None:
			sys.exit('Please provide a vcf file to --vcf argument')
		get_allele_group_prop(options.vcf, options.mat, options.prefix, options.names, options.dGroup, options.ExclChr)
	
	# Get genotype from VCF
	if options.type == 'GET_GENOTYPE_AND_GROUP':
		sys.stdout.write('Associated parameters:\n\t--vcf\n\t--mat\n\t--dGroup\n\t--names\n\t--prefix\n')
		if options.vcf == None:
			sys.exit('Please provide a vcf file to --vcf argument')
		get_genotype_and_group(options.vcf, options.mat, options.dGroup, options.names, options.prefix)
		
if __name__ == "__main__": __main__()