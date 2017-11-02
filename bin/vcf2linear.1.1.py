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
import optparse
import os
import shutil
import subprocess
import tempfile
import fileinput
import time
import random
import math
import datetime
import traceback
import operator
import multiprocessing as mp
from inspect import currentframe, getframeinfo
from operator import itemgetter

import numpy
import matplotlib
matplotlib.use('Agg')
import matplotlib
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.patches as mpatches
from scipy.interpolate import spline
from textwrap import wrap

from scipy.stats import norm

sys.stdout.write('modules loaded\n')

def remove_values_from_list(the_list, val):
   return [value for value in the_list if value != val]

def moyenne(L):
	if len(L) == 0:
		moyenne = 0
	else:
		moyenne = sum(L)/float(len(L))
	return moyenne

def mediane(L):
	L.sort()
	N = len(L)
	n = N/2.0
	p = int(n)
	if n == 0:
		return 0
	if n == 1:
		return (L[0])
	elif n == p:
		return (L[p-1]+L[p])/2.0
	else:
		return L[p]
		
def intervalle(L, P):
	L.sort()
	N = len(L)
	return [L[int((N-(N*P))/2.0)],L[int(N-((N-(N*P))/2.0))]]

def variance(L) :
	n = len(L)
	mq = moyenne(L)**2
	s = sum([x**2 for x in L])
	variance = (n/(n-1))*(s / n - mq)
	return variance
	
def ecart_type(L) :
	VAR = variance(L)
	ecart_type = math.sqrt(VAR)
	return ecart_type

def Normal_prob(MU, MOY, SD):
	
	k= 1/(SD * (math.sqrt(2*math.pi)))
	u = (MU-MOY)/SD
	p = k*math.exp(-0.5*u*u)
	return p

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

def Get_chr_size(VCF):
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
	file.close()
	return dico_chr

def RecordAlleleGrouping(DICO, GROUP, DICO_CHR, GROUPS):
	
	for n in DICO_CHR:
		DICO[n] = {}
	
	file = open(GROUP)
	header = []
	for line in file:
		data = line.split()
		if data:
			if header:
				SNP = data[0].split(':')
				CHR = SNP[0]
				POS = int(SNP[1])
				ALLELE = SNP[2]
				PRESENT = SNP[3]
				if PRESENT == 'P':
					if CHR in DICO:
						if "K-mean_GROUP" in header:
							groupe = data[header.index("K-mean_GROUP")+1]
						elif "GROUP" in header:
							groupe = data[header.index("GROUP")+1]
						else:
							sys.exit("Wrong format for file "+GROUP+". Neither GROUP and K-mean columns were found")
						if groupe in GROUPS:
							if not(POS in DICO[CHR]):
								DICO[CHR][POS] = {}
							DICO[CHR][POS][ALLELE] = groupe
			else:
				header = data
	file.close()

def RecordGrouping(VCF, ACC_TO_DRAW, DICO_CHR, DICO_GROUP, DICO):
	
	file = open(VCF)
	for line in file:
		data = line.split()
		if data:
			if data[0] == "#CHROM":
				header = data
			elif data[0][0] != "#":
				FORMAT = data[header.index("FORMAT")].split(':')
				CHR = data[header.index("#CHROM")]
				POS = int(data[header.index("POS")])
				REF = [data[header.index("REF")]]
				ALT = data[header.index("ALT")].split(',')
				VARIANT = REF + ALT
				if CHR in DICO_GROUP:
					if POS in DICO_GROUP[CHR]:
						for ACC in ACC_TO_DRAW:
							ACCESSION = data[header.index(ACC)].split(':')
							GENOTYPE = ACCESSION[FORMAT.index("GT")].replace('|','/').split('/')
							HETERO = set(GENOTYPE)
							if not ('.' in GENOTYPE):
								if not(ACC in DICO):
									DICO[ACC] = {}
								if not(CHR in DICO[ACC]):
									DICO[ACC][CHR] = {}
								if not(POS in DICO[ACC][CHR]):
									DICO[ACC][CHR][POS] = {}
									DICO[ACC][CHR][POS]['grouping'] = []
									if len(HETERO) > 1:
										DICO[ACC][CHR][POS]['hetero'] = 1
									else:
										DICO[ACC][CHR][POS]['hetero'] = 0
								for allele in GENOTYPE:
									INDEX = int(allele)
									ALLELE = VARIANT[INDEX]
									if ALLELE in DICO_GROUP[CHR][POS]:
										groupe = DICO_GROUP[CHR][POS][ALLELE]
									else:
										groupe = 'NA'
									DICO[ACC][CHR][POS]['grouping'].append(groupe)
	file.close()

def RecordGroup2Hybidise(NAMES):
	# Recording accessions to work with and groups
	DICO_ACC_GP = {}
	file = open(NAMES)
	for line in file:
		data = line.split()
		if data:
			if not (data[1] in DICO_ACC_GP):
				DICO_ACC_GP[data[1]] = []
			DICO_ACC_GP[data[1]].append(data[0])
	file.close()
	return DICO_ACC_GP

def combinliste(seq, k): #  Code from http://python.jpvweb.com/mesrecettespython/doku.php?id=combinaisons
	# Code from http://python.jpvweb.com/mesrecettespython/doku.php?id=combinaisons
	"""Renvoie la liste des combinaisons avec répétition des objets de seq pris k à k"""
	p = []
	i, imax = 0, 2**len(seq)-1
	while i<=imax:
		s = []
		j, jmax = 0, len(seq)-1
		while j<=jmax:
			if (i>>j)&1==1:
				s.append(seq[j])
			j += 1
		if len(s)==k:
			p.append(s)
		i += 1 
	return p

def combinlisterep(seq, k): #  Code from http://python.jpvweb.com/mesrecettespython/doku.php?id=combinaisons
	# Code from http://python.jpvweb.com/mesrecettespython/doku.php?id=combinaisons
	"""Renvoie la liste des combinaisons avec répétition des objets de seq pris k à k"""
	# ajoute chaque objet de seq pour quils apparaissent chacun k fois
	seq2 = []
	for elem in seq:
		if elem not in seq2:
			for i in range(0,k):
				seq2.append(elem)
	# calcule la liste "normale" des combinaisons
	p = combinliste(seq2, k)
	# élimine de cette liste les éléments identiques (comme [1,2] et [1,2])
	p2 = []
	for x in p:
		if x not in p2:
			p2.append(x)
	# et renvoie le résultat
	return p2

def do_cross(DICO_ACC_GP, PLOIDY):
	liste_group = list(DICO_ACC_GP.keys())
	combinaisons = combinlisterep(liste_group, PLOIDY)
	return combinaisons

def count_allele(DICO_CHR_POS, LIST_POSITIONS, GROUPS):
	
	dico_count = {}
	liste_value = []
	for n in LIST_POSITIONS:
		liste_value = liste_value + DICO_CHR_POS[n]['grouping']
	for n in GROUPS:
		dico_count[n] = liste_value.count(n)
	return dico_count

def DoTheHybrid(DICO_gp, CROSS, DICO_ACC_GP, DICO, NB_INDIVIDUALS):

	# recording position shared by all accessions
	dic_pos = {}
	for gp in DICO_ACC_GP:
		for acc in DICO_ACC_GP[gp]:
			for chr in DICO_gp[acc].keys():
				if not(chr in dic_pos):
					dic_pos[chr] = set(DICO_gp[acc][chr].keys())
				else:
					TempoSet = set(DICO_gp[acc][chr].keys())
					dic_pos[chr].intersection_update(TempoSet)
	
	for gpcross in CROSS:
		TheCross = ':'.join(gpcross)
		ind = 0
		while ind < NB_INDIVIDUALS:
			ind += 1
			simul_id = ":".join([TheCross,str(ind)])
			DICO[simul_id] = {}
		
		for chr in dic_pos:
			ind = 0
			while ind < NB_INDIVIDUALS:
				ind += 1
				simul_id = ":".join([TheCross,str(ind)])
				DICO[simul_id][chr] = {}
			
			for pos in dic_pos[chr]:
				ind = 0
				while ind < NB_INDIVIDUALS:
					ind += 1
					simul_id = ":".join([TheCross,str(ind)])
					DICO[simul_id][chr][pos] = {}
					DICO[simul_id][chr][pos]['grouping'] = []
					for gp in gpcross:
						parent = DICO_ACC_GP[gp][random.randint(0,len(DICO_ACC_GP[gp])-1)]
						allele = random.choice(DICO_gp[parent][chr][pos]['grouping'])
						DICO[simul_id][chr][pos]['grouping'].append(allele)
		print(TheCross)

def CalculateHybridProb(DICO_gp, DICO_ACC_GP, DICO, NB_INDIVIDUALS):

	# recording position shared by all accessions
	dic_pos = {}
	for gp in DICO_ACC_GP:
		for acc in DICO_ACC_GP[gp]:
			for chr in DICO_gp[acc].keys():
				if not(chr in dic_pos):
					dic_pos[chr] = set(DICO_gp[acc][chr].keys())
				else:
					TempoSet = set(DICO_gp[acc][chr].keys())
					dic_pos[chr].intersection_update(TempoSet)
	
	
			
	# For a simple verification (of bugs...)
	TOTAL = 0
	for gpcross in DICO_ACC_GP:
		TOTAL += len(DICO_ACC_GP[gpcross])
	
	# Recording informations
	for chr in dic_pos:
		DICO[chr] = {}
		for pos in dic_pos[chr]:
			DICO[chr][pos] = {}
			for gpcross in DICO_ACC_GP:
				DICO[chr][pos][gpcross] = [0,0,0,0,0,0,0,0]
				## [ P; Q; p; q; Pn;Qn;pn;qn] (P;Q;Pn;Qn ==> count, p;q;pn;qn ==> proportions, n ==> noise)
				## [ 0; 1; 2; 3; 4 ;5 ;6 ;7 ]
			for gpcross in DICO_ACC_GP:
				for acc in DICO_ACC_GP[gpcross]:
					group = DICO_gp[acc][chr][pos]['grouping']
					for gp in group:
						if gp == gpcross:
							DICO[chr][pos][gpcross][0] += 1
						else:
							DICO[chr][pos][gpcross][1] += 1
							if gp != 'NA':
								DICO[chr][pos][gp][4] += 1
						for noiseGp in DICO_ACC_GP:
							if noiseGp != gpcross:
								if gp != noiseGp:
									DICO[chr][pos][noiseGp][5] += 1
			for gpcross in DICO_ACC_GP:
				DICO[chr][pos][gpcross][2] = float(DICO[chr][pos][gpcross][0])/(DICO[chr][pos][gpcross][0]+DICO[chr][pos][gpcross][1])
				DICO[chr][pos][gpcross][3] = float(DICO[chr][pos][gpcross][1])/(DICO[chr][pos][gpcross][0]+DICO[chr][pos][gpcross][1])
				DICO[chr][pos][gpcross][6] = float(DICO[chr][pos][gpcross][4])/(DICO[chr][pos][gpcross][4]+DICO[chr][pos][gpcross][5])
				DICO[chr][pos][gpcross][7] = float(DICO[chr][pos][gpcross][5])/(DICO[chr][pos][gpcross][4]+DICO[chr][pos][gpcross][5])
				
				# Verification of bugs
				if len(DICO_ACC_GP[gpcross])*2 != sum(DICO[chr][pos][gpcross][0:2]):
					sys.exit('bug in '+gpcross+' count!!!')
				elif (TOTAL - len(DICO_ACC_GP[gpcross]))*2 != sum(DICO[chr][pos][gpcross][4:6]):
					sys.exit('bug in '+gpcross+' noise count!!!')

def EstimateHybridValueFromSimulation(DICO_HYBRIDS, GROUPS, PLOIDY, CURRENT_POSITION, NB_INDIVIDUALS, CROSS_2_DO, DO_DICO_HYBRID, CHR):
	for gp in GROUPS:
		DICO_HYBRIDS[gp] = {}
		for plo in range(PLOIDY):
			DICO_HYBRIDS[gp]["H"+str(plo+1)] = []
		DICO_HYBRIDS[gp]['noise'] = []
	
	gp_done = set()
	for cross in CROSS_2_DO:
		ind = 0
		while ind < NB_INDIVIDUALS:
			ind += 1
			simul_id = ":".join(cross+[str(ind)])
			SetCross = set(cross)
			for gp in SetCross:
				hybrid_count = count_allele(DO_DICO_HYBRID[simul_id][CHR], CURRENT_POSITION, GROUPS)
				DICO_HYBRIDS[gp]["H"+str(cross.count(gp))].append(hybrid_count[gp])
			for gp_noise in GROUPS:
				if not(gp_noise in SetCross):
					DICO_HYBRIDS[gp_noise]['noise'].append(hybrid_count[gp_noise])

def EstimateHybridValueFromBinomial(GROUPS, DICO_HYBRIDS, PLOIDY, CURRENT_POSITION, DICO_HYBRID_MEAN_AND_VAR, CHR):
	for gp in GROUPS:
		DICO_HYBRIDS[gp] = {}
		for plo in range(PLOIDY):
			DICO_HYBRIDS[gp]["H"+str(plo+1)] = [0,0,0]
			for position in CURRENT_POSITION:
				DICO_HYBRIDS[gp]["H"+str(plo+1)][0] += (DICO_HYBRID_MEAN_AND_VAR[CHR][position][gp][2]*(plo+1))
				DICO_HYBRIDS[gp]["H"+str(plo+1)][1] += (DICO_HYBRID_MEAN_AND_VAR[CHR][position][gp][2]*(plo+1)*DICO_HYBRID_MEAN_AND_VAR[CHR][position][gp][3])
		
		DICO_HYBRIDS[gp]['noise'] = [0,0,0]
		for position in CURRENT_POSITION:
			DICO_HYBRIDS[gp]['noise'][0] += (DICO_HYBRID_MEAN_AND_VAR[CHR][position][gp][6]*PLOIDY)
			DICO_HYBRIDS[gp]['noise'][1] += (DICO_HYBRID_MEAN_AND_VAR[CHR][position][gp][6]*(plo+1)*DICO_HYBRID_MEAN_AND_VAR[CHR][position][gp][7])

def CalcProbFromSimulated(ACC, dico_prop, PREFIX, groups, do_dico_hybrid, WINDOW, cross_2_do, PLOIDY, NB_INDIVIDUALS, DICO_HYBRID_MEAN_AND_VAR, PROPORTION):
	
	# For header display
	HeaderDiplay = []
	for n in groups:
		HeaderDiplay.append(n)
	for i in range(PLOIDY):
		for n in groups:
			HeaderDiplay.append("Loc-mu-H"+str(i+1)+"-"+n)
		for n in groups:
			HeaderDiplay.append("Loc-sd-H"+str(i+1)+"-"+n)
		HeaderDiplay.append("Loc-max-sd-H"+str(i+1))
	HeaderDiplay.append('hetero')
	for i in range(PLOIDY):
		for n in groups:
			HeaderDiplay.append("Prob-H"+str(i+1)+"-"+n)
	for n in groups:
		HeaderDiplay.append("Loc-mu-noise-"+n)
	for n in groups:
		HeaderDiplay.append("Loc-sd-noise-"+n)
	HeaderDiplay.append("Loc-max-sd-noise")
	for n in groups:
		HeaderDiplay.append("Prob-noise-"+n)
	
	# Working chromosome per chromosome
	for chr in dico_prop[ACC]:
		outfile = open(PREFIX+'/'+PREFIX+'_'+ACC+'_'+chr+'.tab', 'w')
		outfile.write('\t'.join(['#CHROM','POS']+HeaderDiplay)+'\n')
		
		current_position = []
		current_hetero = []
		list_position = sorted(list(dico_prop[ACC][chr].keys()))
		for pos in list_position:
			current_position.append(pos)
			current_hetero.append(dico_prop[ACC][chr][pos]['hetero'])
			if len(current_position) == (WINDOW*2 + 1):
				
				# We need to count for each group the number of grouped allele
				acc_gp_count = count_allele(dico_prop[ACC][chr], current_position, groups)
				
				# We need to estimate values on simulated accessions
				dico_hybrids = {}
				EstimateHybridValueFromSimulation(dico_hybrids, groups, PLOIDY, current_position, NB_INDIVIDUALS, cross_2_do, do_dico_hybrid, chr)
				
				
				# Calculating statistics
				dico_stat = {}
				mot_final = []
				for gp in groups:
					dico_stat[gp] = {}
					dico_stat[gp]['count'] = acc_gp_count[gp]
					mot_final.append(dico_stat[gp]['count'])
					
				# Calculating expected values
				for i in range(PLOIDY):
					sd_liste = []
					for gp in groups:
						dico_stat[gp]["Loc-mu-H"+str(i+1)] = moyenne(dico_hybrids[gp]["H"+str(i+1)])
						mot_final.append(dico_stat[gp]["Loc-mu-H"+str(i+1)])
					for gp in groups:
						if PROPORTION:
							dico_stat[gp]["Loc-sd-H"+str(i+1)] = dico_stat[gp]["Loc-mu-H"+str(i+1)]*PROPORTION
						else:
							dico_stat[gp]["Loc-sd-H"+str(i+1)] = ecart_type(dico_hybrids[gp]["H"+str(i+1)])
						sd_liste.append(dico_stat[gp]["Loc-sd-H"+str(i+1)])
						mot_final.append(dico_stat[gp]["Loc-sd-H"+str(i+1)])
					dico_stat["Loc-max-sd-H"+str(i+1)] = max(sd_liste)
					mot_final.append(dico_stat["Loc-max-sd-H"+str(i+1)])
				
				# Calculating heterozygosity value
				dico_stat["hetero"] = sum(current_hetero)/len(current_hetero)
				mot_final.append(dico_stat["hetero"])
				
				# Calculating probabilities
				for i in range(PLOIDY):
					for gp in groups:
						maxi_sd = dico_stat["Loc-max-sd-H"+str(i+1)]
						mu = dico_stat[gp]["Loc-mu-H"+str(i+1)]
						mini = mu-maxi_sd
						maxi = mu+maxi_sd
						value = dico_stat[gp]['count']
						if value == 0:
							dico_stat[gp]["Prob-H"+str(i+1)] = 0
						elif value > maxi:
							dico_stat[gp]["Prob-H"+str(i+1)] =1
						elif value < mini:
							dico_stat[gp]["Prob-H"+str(i+1)] = round(Normal_prob(value, mini, maxi_sd)/Normal_prob(mini, mini, maxi_sd), 7)
						else:
							dico_stat[gp]["Prob-H"+str(i+1)] = 1
						mot_final.append(dico_stat[gp]["Prob-H"+str(i+1)])
				
				# Calculating mean group noise
				for gp in groups:
					dico_stat[gp]["Loc-mu-noise"] = moyenne(dico_hybrids[gp]['noise'])
					mot_final.append(dico_stat[gp]["Loc-mu-noise"])
				
				# Calculating sd group noise
				sd_liste = []
				for gp in groups:
					if PROPORTION:
						dico_stat[gp]["Loc-sd-noise"] = dico_stat[gp]["Loc-mu-noise"]*PROPORTION
					else:
						dico_stat[gp]["Loc-sd-noise"] = ecart_type(dico_hybrids[gp]['noise'])
					mot_final.append(dico_stat[gp]["Loc-sd-noise"])
					sd_liste.append(dico_stat[gp]["Loc-sd-noise"])
				dico_stat["Loc-max-sd-noise"] = max(sd_liste)
				mot_final.append(dico_stat["Loc-max-sd-noise"])
				
				# Calculating Noise probabilities
				for gp in groups:
					maxi_sd = dico_stat["Loc-max-sd-noise"]
					mu = dico_stat[gp]["Loc-mu-noise"]
					mini = mu-maxi_sd
					maxi = mu+maxi_sd
					value = dico_stat[gp]['count']
					if value < mini:
						dico_stat[gp]["Prob-noise"] = 1
					elif value > maxi:
						dico_stat[gp]["Prob-noise"] = Normal_prob(value, maxi, maxi_sd)/Normal_prob(maxi, maxi, maxi_sd)
					else:
						dico_stat[gp]["Prob-noise"] = 1
					mot_final.append(dico_stat[gp]["Prob-noise"])
				
				outfile.write('\t'.join([chr,str(current_position[WINDOW])]+list(map(str, mot_final)))+'\n')
				outfile.flush()
				
				del current_hetero[0]
				del current_position[0]
		outfile.close()
	return 0

def CalcProbFromBinomial(ACC, dico_prop, PREFIX, groups, do_dico_hybrid, WINDOW, cross_2_do, PLOIDY, NB_INDIVIDUALS, DICO_HYBRID_MEAN_AND_VAR, PROPORTION):
	
	# For header display
	HeaderDiplay = []
	for n in groups:
		HeaderDiplay.append(n)
	for i in range(PLOIDY):
		for n in groups:
			HeaderDiplay.append("Loc-mu-H"+str(i+1)+"-"+n)
		for n in groups:
			HeaderDiplay.append("Loc-sd-H"+str(i+1)+"-"+n)
		HeaderDiplay.append("Loc-max-sd-H"+str(i+1))
	HeaderDiplay.append('hetero')
	for i in range(PLOIDY):
		for n in groups:
			HeaderDiplay.append("Prob-H"+str(i+1)+"-"+n)
	for n in groups:
		HeaderDiplay.append("Loc-mu-noise-"+n)
	for n in groups:
		HeaderDiplay.append("Loc-sd-noise-"+n)
	HeaderDiplay.append("Loc-max-sd-noise")
	for n in groups:
		HeaderDiplay.append("Prob-noise-"+n)
	
	# Working chromosome per chromosome
	for chr in dico_prop[ACC]:
		outfile = open(PREFIX+'/'+PREFIX+'_'+ACC+'_'+chr+'.tab', 'w')
		outfile.write('\t'.join(['#CHROM','POS']+HeaderDiplay)+'\n')
		
		current_position = []
		current_hetero = []
		list_position = sorted(list(dico_prop[ACC][chr].keys()))
		for pos in list_position:
			current_position.append(pos)
			current_hetero.append(dico_prop[ACC][chr][pos]['hetero'])
			if len(current_position) == (WINDOW*2 + 1):
				
				# We need to count for each group the number of grouped allele
				acc_gp_count = count_allele(dico_prop[ACC][chr], current_position, groups)
				
				# We need to estimate values on simulated accessions
				dico_hybrids_bis = {}
				EstimateHybridValueFromBinomial(groups, dico_hybrids_bis, PLOIDY, current_position, DICO_HYBRID_MEAN_AND_VAR, chr)
				
				# Calculating statistics
				dico_stat = {}
				mot_final = []
				for gp in groups:
					dico_stat[gp] = {}
					dico_stat[gp]['count'] = acc_gp_count[gp]
					mot_final.append(dico_stat[gp]['count'])
					
				# Calculating expected values
				for i in range(PLOIDY):
					sd_liste = []
					for gp in groups:
						dico_stat[gp]["Loc-mu-H"+str(i+1)] = dico_hybrids_bis[gp]["H"+str(i+1)][0]
						mot_final.append(dico_stat[gp]["Loc-mu-H"+str(i+1)])
					for gp in groups:
						if PROPORTION:
							dico_stat[gp]["Loc-sd-H"+str(i+1)] = dico_stat[gp]["Loc-mu-H"+str(i+1)]*PROPORTION
						else:
							dico_stat[gp]["Loc-sd-H"+str(i+1)] = math.sqrt(dico_hybrids_bis[gp]["H"+str(i+1)][1])
						sd_liste.append(dico_stat[gp]["Loc-sd-H"+str(i+1)])
						mot_final.append(dico_stat[gp]["Loc-sd-H"+str(i+1)])
					if max(sd_liste) == 0:
						sd_liste.append(1e-16)
					dico_stat["Loc-max-sd-H"+str(i+1)] = max(sd_liste)
					mot_final.append(dico_stat["Loc-max-sd-H"+str(i+1)])
					
				# Calculating heterozygosity value
				dico_stat["hetero"] = sum(current_hetero)/len(current_hetero)
				mot_final.append(dico_stat["hetero"])
				
				# Calculating probabilities
				for i in range(PLOIDY):
					for gp in groups:
						maxi_sd = dico_stat["Loc-max-sd-H"+str(i+1)]
						mu = dico_stat[gp]["Loc-mu-H"+str(i+1)]
						mini = mu-maxi_sd
						maxi = mu+maxi_sd
						value = dico_stat[gp]['count']
						if value == 0:
							dico_stat[gp]["Prob-H"+str(i+1)] = 0
						elif value > maxi:
							dico_stat[gp]["Prob-H"+str(i+1)] =1
						elif value < mini:
							dico_stat[gp]["Prob-H"+str(i+1)] = round(Normal_prob(value, mini, maxi_sd)/Normal_prob(mini, mini, maxi_sd), 7)
						else:
							dico_stat[gp]["Prob-H"+str(i+1)] = 1
						mot_final.append(dico_stat[gp]["Prob-H"+str(i+1)])
						
				# Calculating mean group noise
				for gp in groups:
					dico_stat[gp]["Loc-mu-noise"] = dico_hybrids_bis[gp]['noise'][0]
					mot_final.append(dico_stat[gp]["Loc-mu-noise"])
					
				# Calculating sd group noise
				sd_liste = []
				for gp in groups:
					if PROPORTION:
						dico_stat[gp]["Loc-sd-noise"] = dico_stat[gp]["Loc-mu-noise"]*PROPORTION
					else:
						dico_stat[gp]["Loc-sd-noise"] = math.sqrt(dico_hybrids_bis[gp]['noise'][1])
					mot_final.append(dico_stat[gp]["Loc-sd-noise"])
					sd_liste.append(dico_stat[gp]["Loc-sd-noise"])
				if max(sd_liste) == 0:
					sd_liste.append(1e-16)
				dico_stat["Loc-max-sd-noise"] = max(sd_liste)
				mot_final.append(dico_stat["Loc-max-sd-noise"])
				
				# Calculating Noise probabilities
				for gp in groups:
					maxi_sd = dico_stat["Loc-max-sd-noise"]
					mu = dico_stat[gp]["Loc-mu-noise"]
					mini = mu-maxi_sd
					maxi = mu+maxi_sd
					value = dico_stat[gp]['count']
					if value < mini:
						dico_stat[gp]["Prob-noise"] = 1
					elif value > maxi:
						dico_stat[gp]["Prob-noise"] = Normal_prob(value, maxi, maxi_sd)/Normal_prob(maxi, maxi, maxi_sd)
					else:
						dico_stat[gp]["Prob-noise"] = 1
					mot_final.append(dico_stat[gp]["Prob-noise"])
					
				outfile.write('\t'.join([chr,str(current_position[WINDOW])]+list(map(str, mot_final)))+'\n')
				outfile.flush()
				
				del current_hetero[0]
				del current_position[0]
				# sys.exit()
		outfile.close()
	return 0

def CalcProb(ACC, dico_prop, PREFIX, groups, do_dico_hybrid, WINDOW, cross_2_do, PLOIDY, NB_INDIVIDUALS, DICO_HYBRID_MEAN_AND_VAR, PROPORTION):
	
	# For header display
	HeaderDiplay = []
	for n in groups:
		HeaderDiplay.append(n)
	for i in range(PLOIDY):
		for n in groups:
			HeaderDiplay.append("Loc-mu-H"+str(i+1)+"-"+n)
		for n in groups:
			HeaderDiplay.append("Loc-sd-H"+str(i+1)+"-"+n)
		HeaderDiplay.append("Loc-max-sd-H"+str(i+1))
	HeaderDiplay.append('hetero')
	for i in range(PLOIDY):
		for n in groups:
			HeaderDiplay.append("Prob-H"+str(i+1)+"-"+n)
	for n in groups:
		HeaderDiplay.append("Loc-mu-noise-"+n)
	for n in groups:
		HeaderDiplay.append("Loc-sd-noise-"+n)
	HeaderDiplay.append("Loc-max-sd-noise")
	for n in groups:
		HeaderDiplay.append("Prob-noise-"+n)
	
	# Working chromosome per chromosome
	for chr in dico_prop[ACC]:
		outfile = open(PREFIX+'/'+PREFIX+'_'+ACC+'_'+chr+'.tab', 'w')
		outfile.write('\t'.join(['#CHROM','POS']+HeaderDiplay)+'\n')
		
		current_position = []
		current_hetero = []
		list_position = sorted(list(dico_prop[ACC][chr].keys()))
		for pos in list_position:
			current_position.append(pos)
			current_hetero.append(dico_prop[ACC][chr][pos]['hetero'])
			if len(current_position) == (WINDOW*2 + 1):
				
				tps1 = time.time()
				# We need to count for each group the number of grouped allele
				acc_gp_count = count_allele(dico_prop[ACC][chr], current_position, groups)
				tps2 = time.time()
				print('temp bloc1:', tps2-tps1)
				
				
				tps1 = time.time()
				# We need to estimate values on simulated accessions
				dico_hybrids = {}
				EstimateHybridValueFromSimulation(dico_hybrids, groups, PLOIDY, current_position, NB_INDIVIDUALS, cross_2_do, do_dico_hybrid, chr)
				tps2 = time.time()
				print('temp bloc2:', tps2-tps1)
				
				
				###############################################################################################################
				
				tps1 = time.time()
				# We need to estimate values on simulated accessions
				dico_hybrids_bis = {}
				EstimateHybridValueFromBinomial(groups, dico_hybrids_bis, PLOIDY, current_position, DICO_HYBRID_MEAN_AND_VAR, chr)
				tps2 = time.time()
				print('temp bloc2bis:', tps2-tps1)
				
				#--------------------------------------------------------------------------------------------------------------
				
				
				tps1 = time.time()
				# Calculating statistics
				dico_stat = {}
				mot_final = []
				for gp in groups:
					dico_stat[gp] = {}
					dico_stat[gp]['count'] = acc_gp_count[gp]
					mot_final.append(dico_stat[gp]['count'])
				tps2 = time.time()
				print('temp bloc3:', tps2-tps1)
				
				tps1 = time.time()
				# Calculating expected values
				for i in range(PLOIDY):
					sd_liste = []
					for gp in groups:
						dico_stat[gp]["Loc-mu-H"+str(i+1)] = moyenne(dico_hybrids[gp]["H"+str(i+1)])
						mot_final.append(dico_stat[gp]["Loc-mu-H"+str(i+1)])
					for gp in groups:
						if PROPORTION:
							dico_stat[gp]["Loc-sd-H"+str(i+1)] = dico_stat[gp]["Loc-mu-H"+str(i+1)]*PROPORTION
						else:
							dico_stat[gp]["Loc-sd-H"+str(i+1)] = ecart_type(dico_hybrids[gp]["H"+str(i+1)])
						sd_liste.append(dico_stat[gp]["Loc-sd-H"+str(i+1)])
						mot_final.append(dico_stat[gp]["Loc-sd-H"+str(i+1)])
					for gp in groups:
						print(gp, str(i+1), dico_stat[gp]["Loc-mu-H"+str(i+1)], dico_stat[gp]["Loc-sd-H"+str(i+1)], dico_hybrids_bis[gp]["H"+str(i+1)][0], math.sqrt(dico_hybrids_bis[gp]["H"+str(i+1)][1]), dico_hybrids_bis[gp]["H"+str(i+1), dico_hybrids_bis[gp]["H"+str(i+1)][0]])
					dico_stat["Loc-max-sd-H"+str(i+1)] = max(sd_liste)
					mot_final.append(dico_stat["Loc-max-sd-H"+str(i+1)])
				tps2 = time.time()
				print('temp bloc4:', tps2-tps1)
				
				tps1 = time.time()
				# Calculating heterozygosity value
				dico_stat["hetero"] = sum(current_hetero)/len(current_hetero)
				mot_final.append(dico_stat["hetero"])
				tps2 = time.time()
				print('temp bloc5:', tps2-tps1)
				
				tps1 = time.time()
				# Calculating probabilities
				for i in range(PLOIDY):
					for gp in groups:
						maxi_sd = dico_stat["Loc-max-sd-H"+str(i+1)]
						mu = dico_stat[gp]["Loc-mu-H"+str(i+1)]
						mini = mu-maxi_sd
						maxi = mu+maxi_sd
						value = dico_stat[gp]['count']
						if value == 0:
							dico_stat[gp]["Prob-H"+str(i+1)] = 0
						elif value > maxi:
							dico_stat[gp]["Prob-H"+str(i+1)] =1
						elif value < mini:
							dico_stat[gp]["Prob-H"+str(i+1)] = round(Normal_prob(value, mini, maxi_sd)/Normal_prob(mini, mini, maxi_sd), 7)
						else:
							dico_stat[gp]["Prob-H"+str(i+1)] = 1
						mot_final.append(dico_stat[gp]["Prob-H"+str(i+1)])
				tps2 = time.time()
				print('temp bloc6:', tps2-tps1)
				
				tps1 = time.time()
				# Calculating mean group noise
				for gp in groups:
					dico_stat[gp]["Loc-mu-noise"] = moyenne(dico_hybrids[gp]['noise'])
					mot_final.append(dico_stat[gp]["Loc-mu-noise"])
				tps2 = time.time()
				print('temp bloc7:', tps2-tps1)
				
				tps1 = time.time()
				# Calculating sd group noise
				sd_liste = []
				for gp in groups:
					if PROPORTION:
						dico_stat[gp]["Loc-sd-noise"] = dico_stat[gp]["Loc-mu-noise"]*PROPORTION
					else:
						dico_stat[gp]["Loc-sd-noise"] = ecart_type(dico_hybrids[gp]['noise'])
					mot_final.append(dico_stat[gp]["Loc-sd-noise"])
					sd_liste.append(dico_stat[gp]["Loc-sd-noise"])
				for gp in groups:
					print('noise', gp, str(i+1), dico_stat[gp]["Loc-mu-noise"], dico_stat[gp]["Loc-sd-noise"], dico_hybrids_bis[gp]['noise'][0], math.sqrt(dico_hybrids_bis[gp]['noise'][1]), dico_hybrids_bis[gp]['noise'])
				dico_stat["Loc-max-sd-noise"] = max(sd_liste)
				mot_final.append(dico_stat["Loc-max-sd-noise"])
				tps2 = time.time()
				print('temp bloc8:', tps2-tps1)
				
				tps1 = time.time()
				# Calculating Noise probabilities
				for gp in groups:
					maxi_sd = dico_stat["Loc-max-sd-noise"]
					mu = dico_stat[gp]["Loc-mu-noise"]
					mini = mu-maxi_sd
					maxi = mu+maxi_sd
					value = dico_stat[gp]['count']
					if value < mini:
						dico_stat[gp]["Prob-noise"] = 1
					elif value > maxi:
						dico_stat[gp]["Prob-noise"] = Normal_prob(value, maxi, maxi_sd)/Normal_prob(maxi, maxi, maxi_sd)
					else:
						dico_stat[gp]["Prob-noise"] = 1
					mot_final.append(dico_stat[gp]["Prob-noise"])
				tps2 = time.time()
				print('temp bloc9:', tps2-tps1)
				
				tps1 = time.time()
				outfile.write('\t'.join([chr,str(current_position[WINDOW])]+list(map(str, mot_final)))+'\n')
				outfile.flush()
				tps2 = time.time()
				print('temp bloc10:', tps2-tps1)
				
				tps1 = time.time()
				del current_hetero[0]
				del current_position[0]
				tps2 = time.time()
				print('temp bloc11:', tps2-tps1)
				# sys.exit()
		outfile.close()

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

def AttributeGpOnProb(DICO, PLOIDY, CHROM, ACC, DICO_GROUPED, GROUP):
			
		for pos in range(len(DICO[ACC][CHROM]['pos'])):
			liste_values = []
			liste_groupes = []
			ListeBestGroup = []
			for gp in GROUP:
				ProbGroupNoise = DICO[ACC][CHROM]['noise']['Prob'+gp][pos]
				for i in range(PLOIDY):
					liste_groupes.append(gp+':'+str(i+1))
					ProbValue = DICO[ACC][CHROM][gp]['Prob'+str(i+1)][pos]
					if ProbValue > ProbGroupNoise:
						liste_values.append(ProbValue)
					else:
						liste_values.append(0)
			if len(liste_values) != len(liste_groupes):
				sys.exit('Fuck1')
			while len(ListeBestGroup) < PLOIDY:
				if len(liste_values) != len(liste_groupes):
					sys.exit('Fuck2')
				NbAvailable = PLOIDY - len(ListeBestGroup)
				max_value = 0
				max_value = max(liste_values)
				if max_value > 0.1: 
					indices = [i for i, x in enumerate(liste_values) if x == max_value]
					PossibleGroups = []
					to_remove = []
					for Ind in indices:
						ProbGp = liste_groupes[Ind]
						to_remove.append(ProbGp)
						GpToAdd = ProbGp.split(':')[0]
						# To manage already grouped from previous step
						ToAddNumber = int(ProbGp.split(':')[1]) - ListeBestGroup.count(GpToAdd)
						if ToAddNumber > 0:
						# To manage already grouped in this step
							ToAddFinal = ToAddNumber - PossibleGroups.count(GpToAdd)
							if ToAddFinal > 0:
								for k in range(ToAddFinal):
									PossibleGroups.append(GpToAdd)
					if len(PossibleGroups) > NbAvailable:
						print ('Not attributes because ambiguous ', pos, NbAvailable, PossibleGroups, ListeBestGroup)
						for toto in indices:
							print(toto, max_value, liste_groupes[toto])
						while len(ListeBestGroup) < PLOIDY:
							ListeBestGroup.append('un')
					else:
						for k in PossibleGroups:
							ListeBestGroup.append(k)
					# revoving grouped informations
					for rem in to_remove:
						liste_groupes = remove_values_from_list(liste_groupes, rem)
					liste_values = remove_values_from_list(liste_values, max_value)
				else:
					ListeBestGroup.append('un')
			DICO_GROUPED.append(ListeBestGroup)

def GroupGroups(DICO_GROUPED, GROUPS, PLOIDY, DICO, ACC, CHROM, DICO_ORDER, PREFIX):
		
		#Verficication
		if len(DICO[ACC][CHROM]['pos']) != len(DICO_GROUPED):
			sys.exit('Fuck3')
		
		# Counting group values
		ListGlobCount = []
		for gp in GROUPS:
			total = 0
			if gp == 'un':
				total = 0
			else:
				for pos in DICO_GROUPED:
					total += pos.count(gp)
			ListGlobCount.append(total)
		print(GROUPS)
		print(ListGlobCount)
		
		# Initiating HaploCount
		DicoLocCount = {}
		for gp in GROUPS:
			DicoLocCount[gp] = [0] * PLOIDY
		
		# Initiating max Grouped
		ListLocGrouped = [0] * len(GROUPS)
		
		# Grouping ....
		for pos in range(len(DICO[ACC][CHROM]['pos'])):
			list_value = []
			for gp in DICO_GROUPED[pos]:
				list_value.append([gp, ListGlobCount[GROUPS.index(gp)], ListLocGrouped[GROUPS.index(gp)]])
			if len(DICO_ORDER) == 0:
				DICO_ORDER.append([])
				liste_sorted = sorted(list_value, key=operator.itemgetter(1), reverse=True)
				if len(liste_sorted) != PLOIDY:
					sys.exit('Fuck4')
				for k in range(PLOIDY):
					gp = liste_sorted[k][0]
					DICO_ORDER[-1].append(gp)
					if gp != 'un':
						DicoLocCount[gp][k] += 1
						ListLocGrouped[GROUPS.index(gp)] += 1
			else:
				DICO_ORDER.append(['']*PLOIDY)
				availablePos = set(list(range(PLOIDY)))
				liste_sorted = sorted(list_value, key=operator.itemgetter(1), reverse=True)
				if len(liste_sorted) != PLOIDY:
					sys.exit('Fuck5')
				# We need to decide who should be attributed first
				AttributionOrderFirst = []
				AttributionOrderSecond = []
				for k in range(PLOIDY):
					gp = liste_sorted[k][0]
					if gp == 'un':
						AttributionOrderSecond.append(gp)
					else:
						# getting group position in preceding position
						indicesPrec = set([i for i, x in enumerate(DICO_ORDER[-2]) if x == gp])
						# getting final intersection
						FinalInterSect = indicesPrec.intersection(availablePos)
						if len(FinalInterSect) > 0:
							AttributionOrderFirst.append(gp)
							ChosenPos = random.sample(FinalInterSect, 1)[0]
							availablePos.difference_update(set([ChosenPos]))
						else:
							AttributionOrderSecond.append(gp)
						
				AttributionOrder = AttributionOrderFirst+AttributionOrderSecond
				availablePos = set(list(range(PLOIDY)))
				for gp in AttributionOrder:
					TmpDicoLocCount = list(DicoLocCount[gp])
					ChosenPos = ''
					if gp in DICO_ORDER[-2]:
						# getting group position in preceding position
						indicesPrec = set([i for i, x in enumerate(DICO_ORDER[-2]) if x == gp])
						# getting best position acording to the local grouping
						indicesLocal = set([i for i, x in enumerate(DicoLocCount[gp]) if x == max(DicoLocCount[gp])])
						# getting intersection
						intersect = indicesPrec.intersection(indicesLocal)
						# getting final intersection
						FinalInterSect = intersect.intersection(availablePos)
						# If single or multiple possible value we sample one at random
						if len(FinalInterSect) > 0:
							ChosenPos = random.sample(FinalInterSect, 1)[0]
						# Else, we choose based on first precedent position and second based on haplotype containing max of the group and then at random
						else:
							# getting group position in preceding position
							indicesPrec = set([i for i, x in enumerate(DICO_ORDER[-2]) if x == gp])
							# getting final intersection
							FinalInterSect1 = indicesPrec.intersection(availablePos)
							# If single or multiple possible value we sample one at random
							if len(FinalInterSect1) == 1:
								ChosenPos = random.choice(list(FinalInterSect1))
							elif len(FinalInterSect1) > 1:
								max_value = max(TmpDicoLocCount)
								while max_value > 0:
									# getting best position acording to the local grouping
									indicesLocal = set([i for i, x in enumerate(DicoLocCount[gp]) if x == max_value])
									# getting intermediate final intersection
									FinalInterSect2 = indicesLocal.intersection(availablePos)
									# getting final intersection
									FinalInterSect = FinalInterSect1.intersection(FinalInterSect2)
									# If single or multiple possible value we sample one at random
									if len(FinalInterSect) > 0:
										ChosenPos = random.choice(list(FinalInterSect))
										break
									# print(TmpDicoLocCount)
									TmpDicoLocCount.remove(max_value)
									if len(TmpDicoLocCount) == 0:
										break
									max_value = max(TmpDicoLocCount)
							# Else, we choose based on haplotype containing max of the group
							else:
								max_value = max(TmpDicoLocCount)
								while max_value > 0:
									# getting best position acording to the local grouping
									indicesLocal = set([i for i, x in enumerate(DicoLocCount[gp]) if x == max_value])
									# getting final intersection
									FinalInterSect = indicesLocal.intersection(availablePos)
									# If single or multiple possible value we sample one at random
									if len(FinalInterSect) > 0:
										ChosenPos = random.choice(list(FinalInterSect))
										break
									# print(TmpDicoLocCount)
									TmpDicoLocCount.remove(max_value)
									if len(TmpDicoLocCount) == 0:
										break
									max_value = max(TmpDicoLocCount)
						if ChosenPos == '':
							# print('toto', pos, availablePos, len(availablePos))
							ChosenPos = random.choice(list(availablePos))
						availablePos.difference_update(set([ChosenPos]))
					else:
						# getting best position acording to the local grouping
						indicesLocal = set([i for i, x in enumerate(DicoLocCount[gp]) if x == max(DicoLocCount[gp])])
						# getting final intersection
						FinalInterSect = indicesLocal.intersection(availablePos)
						# If single or multiple possible value we sample one at random
						if len(FinalInterSect) > 0:
							ChosenPos = random.choice(list(FinalInterSect))
						# Else, we choose at random...
						else:
							ChosenPos = random.choice(list(availablePos))
						availablePos.difference_update(set([ChosenPos]))
					DICO_ORDER[-1][ChosenPos] = gp
				
				# for counting: not added directly because counting has an impact on rodering at a position
				for k in range(PLOIDY):
					gp = DICO_ORDER[-1][k]
					if gp != 'un':
						DicoLocCount[gp][k] += 1
						ListLocGrouped[GROUPS.index(gp)] += 1
			# print (DICO_ORDER[-1])
			if '' in DICO_ORDER[-1]:
				sys.exit()

def IdentBlocs(DICO_ORDER, DICO, PLOIDY, ACC, CHROM, PREFIX):
	
	for plo in range(PLOIDY):
		outfile = open(PREFIX+'/'+ACC+'_'+CHROM+'_haplo'+str(plo+1)+'.tab','w')
		outfile.close()
		
	ListGpPrec = [0]*PLOIDY
	ListPosStart = []
	ListPosPrec = []
	for i in range(len(DICO[ACC][CHROM]['pos'])):
		pos = DICO[ACC][CHROM]['pos'][i]
		if ListPosPrec == []:
			for plo in range(PLOIDY):
				ListGpPrec[plo] = DICO_ORDER[i][plo]
				ListPosStart.append(pos)
				ListPosPrec.append(pos)
		else:
			for plo in range(PLOIDY):
				if ListGpPrec[plo] != DICO_ORDER[i][plo]:
					outfile = open(PREFIX+'/'+ACC+'_'+CHROM+'_haplo'+str(plo+1)+'.tab','a')
					outfile.write('\t'.join([ACC, CHROM, str(ListPosStart[plo]), str(ListPosPrec[plo]), ListGpPrec[plo]]))
					outfile.write('\n')
					outfile.close()
					ListPosStart[plo] = pos
					ListGpPrec[plo] = DICO_ORDER[i][plo]
				ListPosPrec[plo] = pos
	for plo in range(PLOIDY):
		outfile = open(PREFIX+'/'+ACC+'_'+CHROM+'_haplo'+str(plo+1)+'.tab','a')
		outfile.write('\t'.join([ACC, CHROM, str(ListPosStart[plo]), str(ListPosPrec[plo]), ListGpPrec[plo]]))
		outfile.write('\n')

def draw_plot(DICO, PREFIX, GCOL, GROUPS, ACC, CHROM, CHR_LEN, PLOIDY):
	
	# color for the plot
	list_color = [(1,0,0), (0,1,0), (0,0,1), (0.87,0.43,0.08), (0.5,0,0.5)]
	
	# line type
	list_line = ['--', '-', ':']
	
	
	# recording color to draw
	dico_color = {}
	get_color(GCOL, GROUPS, dico_color)
	
	X = 6 + len(GROUPS)*10 + 10 + 3
	POSSPAN = 0
	
	# Drawing heterozygosity plot
	fig = plt.figure(figsize=(10.5, 14.85))
	fig.subplots_adjust(left=0.05, right=0.95, top=0.98, bottom=0.02)
	ax = plt.subplot2grid((X,15),(POSSPAN,0), colspan=13, rowspan=6)
	ax.set_ylim(0, 1.05)
	ax.set_xlim(0, CHR_LEN)
	ax.plot(DICO[ACC][CHROM]['pos'], DICO[ACC][CHROM]['hetero'], linewidth=1, color='black', linestyle='-')
	ax.axes.xaxis.set_ticklabels([])
	ax.tick_params(axis='y', labelsize=5)
	ax.set_title(ACC+' '+CHROM, fontweight='bold', position=(0.5, 1))
	
	ax = plt.subplot2grid((X,15),(0,13), colspan=2, rowspan=3)
	plt.axis([0, 1, 0, 1])
	ax.axis('off')
	ax.text(0.0, 0.5, "\n".join(wrap('Heterozygosity level', 20)), size=9, va='center')
	POSSPAN += 6
	
	# Drawing groups
	for gp in GROUPS:
		ax = plt.subplot2grid((X,15),(POSSPAN,0), colspan=13, rowspan=10)
		ax.set_xlim(0, CHR_LEN)
		ax.set_ylim(0, max([max(DICO[ACC][CHROM][gp]['max'+str(PLOIDY)]),max(DICO[ACC][CHROM]['count'][gp])])*1.05)
		for i in range(PLOIDY):
			ax.plot([],[], color=list_color[i], alpha=0.5, label=gp+' H'+str(i+1)+' interval', linewidth=10)
			ax.fill_between(DICO[ACC][CHROM]['pos'], DICO[ACC][CHROM][gp]['min'+str(i+1)],	DICO[ACC][CHROM][gp]['max'+str(i+1)],		color=list_color[i],	alpha=0.5, label=gp+' H'+str(i+1)+' interval')
		ax.plot([],[], color='grey',	alpha=0.5, label=gp+' Noise interval',	linewidth=10)
		
		ax.fill_between(DICO[ACC][CHROM]['pos'], DICO[ACC][CHROM]['noise']['min'+gp],	DICO[ACC][CHROM]['noise']['max'+gp],	color='grey',	alpha=0.5, label=gp+' noise interval')
		ax.plot(DICO[ACC][CHROM]['pos'], DICO[ACC][CHROM]['count'][gp], linewidth=1, linestyle='-', color='black', label=gp+' accession count')
		
		ax.axes.xaxis.set_ticklabels([])
		ax.tick_params(axis='y', labelsize=3)
		ax.legend(bbox_to_anchor=(1.005, 1), loc=2, borderaxespad=0., fontsize=9)
		POSSPAN += 10
	
	# Drawing probabilities
	ax = plt.subplot2grid((X,15),(POSSPAN,0), colspan=13, rowspan=10)
	for i in range(PLOIDY):
		for gp in GROUPS:
			ax.plot(DICO[ACC][CHROM]['pos'], DICO[ACC][CHROM][gp]['Prob'+str(i+1)], linewidth=1, linestyle=list_line[i], color=dico_color[gp][0:3], label=gp+' H'+str(i+1)+' Prob.')
	ax.set_ylim(-0.01, 1.05)
	ax.set_xlim(0, CHR_LEN)
	ax.axes.xaxis.set_ticklabels([])
	ax.tick_params(axis='y', labelsize=5)
	ax.legend(bbox_to_anchor=(1.005, 1), loc=2, borderaxespad=0., fontsize=9)
	POSSPAN += 10
	
	# Group Attribution based on maximal probabilities
	dico_grouped = []
	AttributeGpOnProb(DICO, PLOIDY, CHROM, ACC, dico_grouped, GROUPS)
	
	# It's time to order groups
	dico_order = []
	GroupGroups(dico_grouped, GROUPS+['un'], PLOIDY, DICO, ACC, CHROM, dico_order, PREFIX)
	
	# Now we identify blocs
	IdentBlocs(dico_order, DICO, PLOIDY, ACC, CHROM, PREFIX)
	
	
	# Now it is time to draw
	## Initiating graphe
	ax = plt.subplot2grid((X,15),(POSSPAN,0), colspan=13, rowspan=3)
	ax.set_ylim(0, PLOIDY)
	ax.set_xlim(0, CHR_LEN)
	ax.axes.xaxis.set_ticklabels([])
	ax.xaxis.set_ticks_position('none') 
	ax.axes.yaxis.set_ticklabels([])
	ax.yaxis.set_ticks_position('none')
	
	## getting informations and drawing
	for plo in range(PLOIDY):
		file = open(PREFIX+'/'+ACC+'_'+CHROM+'_haplo'+str(plo+1)+'.tab', 'r')
		# recording information
		for line in file:
			data = line.split()
			if data:
				x0 = int(data[2])
				y0 = plo
				gp = data[4]
				Xd = (int(data[3])-x0)+1
				Yd = 0.80
				if gp in dico_color:
					col = dico_color[gp][0:3]
				else:
					col = (0.5,0.5,0.5)
				ax.add_patch(mpatches.Rectangle((x0,y0), Xd, Yd, color=col))
		file.close()
	
	fig.savefig(PREFIX+'/'+ACC+'_'+CHROM+'_density.pdf')
	plt.close(fig)

def Record_data(GROUPS, ACC_TO_DRAW, DIC_CHR, PLOIDY, PREFIX, DICO2DRAW):
	
	for acc in ACC_TO_DRAW:
		DICO2DRAW[acc] = {}
		for chr in DIC_CHR:
			DICO2DRAW[acc][chr] = {}
			file = open(PREFIX+'/'+PREFIX+'_'+acc+'_'+chr+'.tab')
			header = file.readline().split()
			DICO2DRAW[acc][chr]['pos'] = []
			DICO2DRAW[acc][chr]['hetero'] = []
			DICO2DRAW[acc][chr]['count'] = {}
			DICO2DRAW[acc][chr]['noise'] = {}
			for n in GROUPS:
				DICO2DRAW[acc][chr]['count'][n] = []
				DICO2DRAW[acc][chr]['noise']['Prob'+n] = []
				DICO2DRAW[acc][chr]['noise']['min'+n] = []
				DICO2DRAW[acc][chr]['noise']['max'+n] = []
			for n in GROUPS:
				DICO2DRAW[acc][chr][n] = {}
				for i in range(PLOIDY):
					DICO2DRAW[acc][chr][n]['min'+str(i+1)] = []
					DICO2DRAW[acc][chr][n]['max'+str(i+1)] = []
					DICO2DRAW[acc][chr][n]['Prob'+str(i+1)] = []
			
			line = file.readline().split()
			while line:
				pos = int(line[header.index('POS')])
				DICO2DRAW[acc][chr]['pos'].append(pos)
				DICO2DRAW[acc][chr]['hetero'].append(float(line[header.index('hetero')]))
				# filling accession count
				for n in GROUPS:
					DICO2DRAW[acc][chr]['count'][n].append(float(line[header.index(n)]))
				
				# filling noise expectations
				for n in GROUPS:
					DICO2DRAW[acc][chr]['noise']['Prob'+n].append(float(line[header.index('Prob-noise-'+n)]))
					DICO2DRAW[acc][chr]['noise']['min'+n].append(float(line[header.index('Loc-mu-noise-'+n)])-float(line[header.index('Loc-max-sd-noise')]))
					DICO2DRAW[acc][chr]['noise']['max'+n].append(float(line[header.index('Loc-mu-noise-'+n)])+float(line[header.index('Loc-max-sd-noise')]))
				
				# filling grouped values expectations
				for i in range(PLOIDY):
					for n in GROUPS:
						DICO2DRAW[acc][chr][n]['min'+str(i+1)].append(float(line[header.index('Loc-mu-H'+str(i+1)+'-'+n)])-float(line[header.index('Loc-max-sd-H'+str(i+1))]))
						DICO2DRAW[acc][chr][n]['max'+str(i+1)].append(float(line[header.index('Loc-mu-H'+str(i+1)+'-'+n)])+float(line[header.index('Loc-max-sd-H'+str(i+1))]))
						DICO2DRAW[acc][chr][n]['Prob'+str(i+1)].append(float(line[header.index('Prob-H'+str(i+1)+'-'+n)]))
				line = file.readline().split()
			file.close()
			
			# Adding start to the dictionnary
			DICO2DRAW[acc][chr]['pos'].insert(0, 1)
			DICO2DRAW[acc][chr]['hetero'].insert(0, DICO2DRAW[acc][chr]['hetero'][0])
			for n in GROUPS:
				DICO2DRAW[acc][chr]['count'][n].insert(0, DICO2DRAW[acc][chr]['count'][n][0])
				DICO2DRAW[acc][chr]['noise']['Prob'+n].insert(0, DICO2DRAW[acc][chr]['noise']['Prob'+n][0])
				DICO2DRAW[acc][chr]['noise']['min'+n].insert(0, DICO2DRAW[acc][chr]['noise']['min'+n][0])
				DICO2DRAW[acc][chr]['noise']['max'+n].insert(0, DICO2DRAW[acc][chr]['noise']['max'+n][0])
			for i in range(PLOIDY):
				for n in GROUPS:
					DICO2DRAW[acc][chr][n]['min'+str(i+1)].insert(0, DICO2DRAW[acc][chr][n]['min'+str(i+1)][0])
					DICO2DRAW[acc][chr][n]['max'+str(i+1)].insert(0, DICO2DRAW[acc][chr][n]['max'+str(i+1)][0])
					DICO2DRAW[acc][chr][n]['Prob'+str(i+1)].insert(0, DICO2DRAW[acc][chr][n]['Prob'+str(i+1)][0])
			
			# adding end to the dictionnary
			DICO2DRAW[acc][chr]['pos'].append(DIC_CHR[chr])
			DICO2DRAW[acc][chr]['hetero'].append(DICO2DRAW[acc][chr]['hetero'][-1])
			for n in GROUPS:
				DICO2DRAW[acc][chr]['count'][n].append(DICO2DRAW[acc][chr]['count'][n][-1])
				DICO2DRAW[acc][chr]['noise']['Prob'+n].append(DICO2DRAW[acc][chr]['noise']['Prob'+n][-1])
				DICO2DRAW[acc][chr]['noise']['min'+n].append(DICO2DRAW[acc][chr]['noise']['min'+n][-1])
				DICO2DRAW[acc][chr]['noise']['max'+n].append(DICO2DRAW[acc][chr]['noise']['max'+n][-1])
			for i in range(PLOIDY):
				for n in GROUPS:
					DICO2DRAW[acc][chr][n]['min'+str(i+1)].append(DICO2DRAW[acc][chr][n]['min'+str(i+1)][-1])
					DICO2DRAW[acc][chr][n]['max'+str(i+1)].append(DICO2DRAW[acc][chr][n]['max'+str(i+1)][-1])
					DICO2DRAW[acc][chr][n]['Prob'+str(i+1)].append(DICO2DRAW[acc][chr][n]['Prob'+str(i+1)][-1])

def CalcGroupProp(VCF, NAMES, NAMES2, PREFIX, CHR, WINDOW, GROUP, GCOL, PLOIDY, THREAD, TYPE, PROPORTION):
	
	"""
		Identify genome structure
		
		:param VCF: A vcf file.
		:type VCF: vcf
		:param NAMES: Path to file containing accession to treat.
		:type NAMES: str
		:param NAMES2: Path to file containing accession to treat.
		:type NAMES2: str
		:param PREFIX: Prefix of the output file.
		:type PREFIX: str
		:param CHR: A string containing chromosomes to work with.
		:type CHR: str
		:param WINDOW: Half widow size to validate genome genotype.
		:type WINDOW: int
		:param GROUP: Grouping information file either using a priori or with PCA.
		:type GROUP: str
		:param GCOL: Color for groups file
		:type GCOL: str
		:param PLOIDY: Accessions ploidy level
		:type PLOIDY: int
		:param THREAD: Number of processors available
		:type THREAD: int
		:param TYPE: Type of probability calculation
		:type TYPE: str
		:param PROPORTION: Proportion of the expected value to define the probabilities
		:type PROPORTION: flaot
		:return: Draw circos picture and return a bloc file identifying genome structure.
		:rtype: void
	"""
	
	
	# individual simulated
	nb_individuals = float(100)
	
	# For output files
	if not(os.path.isdir(PREFIX)):
		os.mkdir(PREFIX)
	
	# Recording chromosomes
	if CHR != None:
		chr_to_work = CHR.split(':')
		print (chr_to_work)
	
	# recording accession names to draw
	acc_to_draw = RecordAccession(VCF, NAMES, None)
	acc_to_draw2 = RecordAccession(VCF, NAMES2, None)
	acc_to_record = list(set(acc_to_draw + acc_to_draw2))
	print (acc_to_draw)
	print (acc_to_draw2)
	
	# recording accession groups which will be hybridized
	dico_acc_gp = RecordGroup2Hybidise(NAMES2)
	print (dico_acc_gp)
	
	# Recording hybrid cross to perform
	cross_2_do = do_cross(dico_acc_gp, PLOIDY)
	print (cross_2_do)
	
	# recording chromosomes to work with
	dic_chr = Get_chr_size(VCF)
	if CHR != None:
		chr2remove = []
		for n in dic_chr:
			if not(n in chr_to_work):
				chr2remove.append(n)
		for n in chr2remove:
			del dic_chr[n]
	print (dic_chr)
	
	# recording group to work with
	groups = list(dico_acc_gp.keys())
	print (groups)
	
	# recording allele grouping
	print ("Recording allele grouping")
	dico_group = {}
	RecordAlleleGrouping(dico_group, GROUP, dic_chr, groups)
	
	# Recording by accessions grouping information
	print ("Recording accession grouping")
	dico_prop = {}
	RecordGrouping(VCF, acc_to_record, dic_chr, dico_group, dico_prop)
	
	# Creating the hybrids
	print ("Creating hybrids")
	do_dico_hybrid = {}
	dico_hybrid_mean_and_var = {}
	if TYPE == "Simul":
		DoTheHybrid(dico_prop, cross_2_do, dico_acc_gp, do_dico_hybrid, nb_individuals)
	elif TYPE == "Binom":
		CalculateHybridProb(dico_prop, dico_acc_gp, dico_hybrid_mean_and_var, nb_individuals)
	
	# sys.exit()
	
	# Working accession per accession
	print("Calculating probabilities")
	listJobs = []
	for acc in acc_to_draw:
		# CalcProb(acc, dico_prop, PREFIX, groups, do_dico_hybrid, WINDOW, cross_2_do, PLOIDY, nb_individuals, dico_hybrid_mean_and_var, PROPORTION)
		# CalcProbFromBinomial(acc, dico_prop, PREFIX, groups, do_dico_hybrid, WINDOW, cross_2_do, PLOIDY, nb_individuals, dico_hybrid_mean_and_var, PROPORTION)
		listJobs.append([TYPE, acc, dico_prop, PREFIX, groups, do_dico_hybrid, WINDOW, cross_2_do, PLOIDY, nb_individuals, dico_hybrid_mean_and_var, PROPORTION])
	pool = mp.Pool(processes=THREAD)
	results = pool.map(run_mutithread, listJobs)
	for n in results:
		if n != 0:
			sys.stdout.write(str(n)+'\n')
			sys.stdout.flush()
	
	# We need to record data probabilities
	print ("Recording informations for drawing")
	dico2draw = {}
	Record_data(groups, acc_to_draw, dic_chr, PLOIDY, PREFIX, dico2draw)
	
	# Plotting informations
	for acc in acc_to_draw:
		for chr in dic_chr:
			draw_plot(dico2draw, PREFIX, GCOL, sorted(groups), acc, chr, dic_chr[chr], PLOIDY)

def run_mutithread(job):
	
	try:
		if job[0] == "Simul":
			rslt = CalcProbFromSimulated(job[1],job[2],job[3],job[4],job[5],job[6], job[7], job[8], job[9], job[10], job[11])
		elif job[0] == "Binom":
			rslt = CalcProbFromBinomial(job[1],job[2],job[3],job[4],job[5],job[6], job[7], job[8], job[9], job[10], job[11])
		else:
			raise ErrorValue ('Unknown function name: '+str(job[0]))
	except Exception as e:
		print (e)
		rslt = 1
	finally:
		return rslt

def __main__():
	#Parse Command Line
	parser = optparse.OptionParser(usage="python %prog [options]\n\nProgram designed by Guillaume MARTIN : guillaume.martin@cirad.fr")
	# Wrapper options. 
	parser.add_option( '-v',	'--vcf',			dest='vcf',			default=None,			help='The vcf file. [Default: %default]')
	parser.add_option( '-m',	'--mat',			dest='mat',			default=None,			help='Grouping information file. [Default: %default]')
	parser.add_option( '-n',	'--names',			dest='names',		default=None,			help='A one column file containing accession names to treat. [Default: %default]')
	parser.add_option( '-N',	'--namesH',			dest='namesH',		default=None,			help='A one column file containing accession names used to simulate populations. [Default: %default]')
	parser.add_option( '-c',	'--chr',			dest='chr',			default=None, 			help='Chromosome names to exclude from analysis. Each chromosomes should be separated by ":". [Default: %default]')
	parser.add_option( '-w',	'--win',			dest='win',			default='25', 			help='Half window size around a variant site to evaluate the structure at the site. [Default: %default]')
	parser.add_option( '-g',	'--gcol',			dest='gcol',		default=None, 			help='Group color. [Default: %default]')
	parser.add_option( '-P',	'--ploidy',			dest='ploidy',		default='2', 			help='Ploidy level (integer). [Default: %default]')
	parser.add_option( '-t',	'--thread',			dest='thread',		default='1',			help='Number of processors to use (integer), [default: %default]')
	parser.add_option( '-T',	'--type',			dest='type',		default='Binom',		help='Type of estimation performed: "Simul", "Binom". If "Simul", a total of 100 individuals are simulated for '
	'each combinations of haplotype and mean values and sd values are estimated based on these simulation. If "Binom", mean value is calculated as the sum of binomial mean at each point (exact estimator) and '
	'sd value is estimated as sqrt(sum variance at each point). This is not the exact sd but the analysis is a lot more faster! If you do not trust this sd estimation, you can choose to change this estimator '
	' by filling a value between ]0,1] to --prop parameter. In this case, the program will use the maximal_expected_value*prop_argument instead of using the maximal sd observed for all groups for probability '
	'calculation [default: %default]')
	parser.add_option( '-d',	'--prop',			dest='prop',		default='0',			help='Estimator different from sd calculated as mean_value*--prop. Value should be comprised in ]0,1]. A value of '
	'0, means that this parameter is not used. [default: %default]')
	parser.add_option( '-p',	'--prefix',			dest='prefix',		default='WorkOnVcf', 	help='The prefix for output files. [Default: %default]')

	(options, args) = parser.parse_args()
	
	
	# Identify genome blocs in accessions and draw a circos by accessions
	if options.vcf == None:
		sys.exit('Please provide a vcf file to --vcf argument')
	if options.mat == None:
		sys.exit('Please provide a matrix file to --mat argument')
	if options.namesH == None:
		sys.exit('Please provide a namesH file to --namesH argument')
	if options.gcol == None:
		sys.exit('Please provide a gcol file to --gcol argument')
	CalcGroupProp(options.vcf, options.names, options.namesH, options.prefix, options.chr, int(options.win), options.mat, options.gcol, int(options.ploidy), int(options.thread), options.type, float(options.prop))
		
if __name__ == "__main__": __main__()