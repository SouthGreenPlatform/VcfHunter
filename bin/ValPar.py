#!/usr/bin/env python
#
#  Copyright 2014 CIRAD
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

import argparse
import sys
import os
import gzip



#### Section of code from http://python.jpvweb.com/python/mesrecettespython/doku.php?id=combinaisons
def combinliste(seq, k):
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

def combinlisterep(seq, k):
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

#### End of section of code from http://python.jpvweb.com/python/mesrecettespython/doku.php?id=combinaisons

def combinGui(seq):
	seq2 = []
	for n in seq:
		for k in seq:
			seq2.append([n,k])
	
	for i in range(len(seq2)): 
		seq2[i] = sorted(seq2[i])
	
	Set2 = set()
	for n in seq2:
		Set2.add(':'.join(n))
	
	seq2 = []
	for n in Set2:
		seq2.append(n.split(':'))
	
	return seq2

def __main__():
	#Parse Command Line
	parser = argparse.ArgumentParser(
		description="This program look for direct parentage between genotyped individuals.",
		epilog="Program designed by Guillaume MARTIN (guillaume.martin@cirad.fr)")
	# Wrapper options.
	parser.add_argument( '-p',	'--parent',	dest='parent',	default=None,			help='Path to one column file containing parent names')
	parser.add_argument( '-v',	'--vcf',	dest='vcf',		default=None,			help='Path to the vcf file')
	parser.add_argument( '-a',	'--acc',	dest='acc',		default='all',			help='Accession name to calculate parentage. Only one accession allowed')
	parser.add_argument( '-o',	'--output',	dest='output',	default='Valid.tab',	help='Output file name')
	options = parser.parse_args()
	
	dicoPlo = {}
	dicoPlo[2] = [[1,1]]
	dicoPlo[3] = [[1,2]]
	dicoPlo[4] = [[1,3],[2,2]]
	dicoPlo[5] = [[1,4],[2,3]]
	dicoPlo[6] = [[1,5],[2,4],[3,3]]
	
	
	if options.acc == 'all':
		sys.exit('Please provide an accession for this analysis\n')
	else:
		CHILD = options.acc
	
	# Loading parent
	file = open(options.parent)
	MyParent = set()
	for line in file:
		data = line.split()
		if data:
			if data[0] != CHILD:
				MyParent.add(data[0])
	file.close()
	
	# Creating combination of parents
	Combin = combinGui(list(MyParent))
	combinNumber = len(Combin)
	
	# Calculating parent identity based on the vcf
	if options.vcf[-3:] == '.gz':
		file = gzip.open(options.vcf,'rt')
	else:
		file = open(options.vcf)
	
	for line in file:
		data = line.split()
		if data:
			if data[0] == "#CHROM":
				header = data[:]
				
				# Recording parent position 
				dicoParentVcfPos = {}
				for Par in MyParent:
					dicoParentVcfPos[Par] = header.index(Par)
					
				# Recording individuals position
				dicoIndVcfPos = {}
				if options.acc == 'all':
					sys.exit('Please provide an accession for this analysis\n')
				else:
					FormatPos = header.index('FORMAT')
					dicoIndVcfPos[CHILD] = header.index(CHILD)
				
				
				dicoIndOrder = {}
				dicoCrossOrder = {}
			
			# Calculating identity
			elif data[0][0] != "#":
				FORMAT = data[FormatPos].split(':')
				GTPos = FORMAT.index('GT')
				# print(data[0], data[1])
				# Obtaining accession information
				GTacc = data[dicoIndVcfPos[CHILD]].split(':')[GTPos].split('/')
				PLO = len(GTacc)
				GT = ':'.join(sorted(GTacc))
				
				if not (dicoCrossOrder):
					# Initiating crosses tested
					TotalCombineNum = 0
					for i in range(combinNumber):
						for k in range(len(dicoPlo[PLO])):
							dicoCrossOrder[':'.join([Combin[i][0]+'-'+str(dicoPlo[PLO][k][0]), Combin[i][1]+'-'+str(dicoPlo[PLO][k][1])])] = TotalCombineNum
							TotalCombineNum += 1
							if Combin[i][0] != Combin[i][1] and dicoPlo[PLO][k][1] != dicoPlo[PLO][k][0]:
								dicoCrossOrder[':'.join([Combin[i][0]+'-'+str(dicoPlo[PLO][k][1]), Combin[i][1]+'-'+str(dicoPlo[PLO][k][0])])] = TotalCombineNum
								TotalCombineNum += 1
					
					# Initiating for variable calculation
					listInd = list(dicoIndVcfPos.keys())
					for n in range(len(listInd)):
						dicoIndOrder[listInd[n]] = n
					
					List2Print1 = [0]*TotalCombineNum*len(dicoIndOrder)
					List2Print2 = [0]*TotalCombineNum*len(dicoIndOrder)
					accessionNumber = len(dicoIndOrder)
				
				
				
				# Identification of parent genotypes and creating their gametes
				DicoParentGeno = {}
				for Par in MyParent:
					GTPar = list(set(data[dicoParentVcfPos[Par]].split(':')[GTPos].split('/')))
					for plo in dicoPlo[PLO]:
						DicoParentGeno[Par+'-'+str(plo[0])] = combinlisterep(GTPar, plo[0])
						DicoParentGeno[Par+'-'+str(plo[1])] = combinlisterep(GTPar, plo[1])
				
				# Creating potential parental combination
				DicoParentCombin = {}
				for comb in Combin:
					# print('**********************')
					for plo in dicoPlo[PLO]:
						geno = set()
						Parent1 = comb[0]+'-'+str(plo[0])
						Parent2 = comb[1]+'-'+str(plo[1])
						for gam1 in DicoParentGeno[Parent1]:
							for gam2 in DicoParentGeno[Parent2]:
								genotype = ':'.join(sorted(gam1+gam2))
								if not('.' in genotype):
									geno.add(genotype)
						DicoParentCombin[':'.join([Parent1, Parent2])] = geno
						# print(Parent1, Parent2, geno, DicoParentGeno[Parent1], DicoParentGeno[Parent2])
						
						geno = set()
						if comb[0] != comb[1]:
							Parent1 = comb[0]+'-'+str(plo[1])
							Parent2 = comb[1]+'-'+str(plo[0])
							for gam1 in DicoParentGeno[Parent1]:
								for gam2 in DicoParentGeno[Parent2]:
									genotype = ':'.join(sorted(gam1+gam2))
									if not('.' in genotype):
										geno.add(genotype)
							DicoParentCombin[':'.join([Parent1, Parent2])] = geno
							# print(Parent1, Parent2, geno, DicoParentGeno[Parent1], DicoParentGeno[Parent2])
				
				# Calculating identity
				for acc in dicoIndVcfPos:
					GT = ':'.join(sorted(data[dicoIndVcfPos[acc]].split(':')[GTPos].split('/')))
					if not('.' in GT):
						for cross in DicoParentCombin:
							if DicoParentCombin[cross]:
								CrossId = dicoCrossOrder[cross]
								IndId = dicoIndOrder[acc]
								ValuePos = (IndId*TotalCombineNum)+CrossId
								# print(ValuePos)
								if GT in DicoParentCombin[cross]:
									# print(acc, cross, 'yes', GT, DicoParentCombin[cross])
									List2Print1[ValuePos] += 1
									List2Print2[ValuePos] += 1
								else:
									# print(acc, cross, 'no', GT, DicoParentCombin[cross])
									List2Print2[ValuePos] += 1
							else:
								pass
	# print(List2Print1)
	# print(List2Print2)
				# sys.exit()
	file.close()
	
	# Printing results
	outfile = open(options.output,'w')
	Mot2print = ['']*TotalCombineNum
	for i in dicoCrossOrder:
		Mot2print[dicoCrossOrder[i]] = i
	Mot2print = ['Ind\Cross']+Mot2print
	outfile.write('\t'.join(Mot2print))
	outfile.write('\n')
	
	for acc in dicoIndVcfPos:
		Mot2print = [acc]
		IndId = dicoIndOrder[acc]
		for k in range(IndId*TotalCombineNum, (IndId+1)*TotalCombineNum):
			value1 = List2Print1[k]
			value2 = List2Print2[k]
			if value2 != 0:
				Mot2print.append(str(value1/value2))
			else:
				Mot2print.append('NA')
		outfile.write('\t'.join(Mot2print))
		outfile.write('\n')
	outfile.close()

if __name__ == "__main__": __main__()
