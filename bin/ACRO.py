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

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


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



def __main__():
	#Parse Command Line
	parser = argparse.ArgumentParser(
		description="This program calculated the number and proportion of sites in accordance with "
		"parentage between two potential parents and a child or one potential parent and a child. "
		"This program does not work on phased vcf file.",
		epilog="Program designed by Guillaume MARTIN (guillaume.martin@cirad.fr)")
	# Wrapper options.
	parser.add_argument( '-1',	'--parent1',	dest='parent1',	default=None,			help='Potential parent 1 and ploidy of the generated gamete. Parent and ploidy should be separated by ","')
	parser.add_argument( '-2',	'--parent2',	dest='parent2',	default=None,			help='Potential parent 2 and ploidy of the generated gamete. Parent and ploidy should be separated by ","')
	parser.add_argument( '-v',	'--vcf',		dest='vcf',		default=None,			help='Path to the vcf file')
	parser.add_argument( '-a',	'--acc',		dest='acc',		default=None,			help='Accession name to calculate parentage. Only one accession allowed')
	parser.add_argument( '-w', '--window',		dest='window',	default='100',			help='Half window to calculate site proportion in accordance with tested parentage (unit = variant site)')
	parser.add_argument( '-W', '--WINDOW',		dest='WINDOW',	default='100000',		help='The window to calculate site proportion in accordance with tested parentage (unit base pair)')
	parser.add_argument( '-f', '--fasta',		dest='fasta',	default=None,			help='The multifasta reference file')
	parser.add_argument( '-c', '--chr',			dest='chr',		default=None,			help='Chromosome to work with. They should be separated by ":". If omitted, all chromosomes will be used')
	parser.add_argument( '-p',	'--prefix',		dest='prefix',	default='Valid',		help='Prefix for output files')
	options = parser.parse_args()
	
	# recording window size
	if int(options.window) < 1:
		sys.exit("There is a problem, the program could not work with -w argument lower than 1")
	WINSNP = int(options.window)*2 + 1
	HalfWINSNP = int(options.window)
	WIN = int(options.WINDOW)
	
	# recording chromosomes to omit
	if options.chr == None:
		ChrToWorkWith = []
	else:
		ChrToWorkWith = options.chr.split(':')
	
	print("Obtaining chromosome statistics")
	# Getting sequence statistics
	dico_chr = {}
	dico_chr_data = {}
	sequence_dict = SeqIO.index(options.fasta, "fasta")
	ChrSet = set()
	if ChrToWorkWith:
		for n in ChrToWorkWith:
			print("--Loading "+n)
			dico_chr[n] = [0]*len(str(sequence_dict[n].seq))
			dico_chr_data[n] = [0]*len(str(sequence_dict[n].seq))
			ChrSet.add(n)
	else:
		for n in sequence_dict:
			print("--Loading "+n)
			dico_chr[n] = [0]*len(str(sequence_dict[n].seq))
			dico_chr_data[n] = [0]*len(str(sequence_dict[n].seq))
			ChrSet.add(n)
	del sequence_dict
	
	print("Working on the vcf file")
	if options.parent2 != None:
	
		# Loading parents information
		P1name,P1Plo = options.parent1.split(',')
		P2name,P2Plo = options.parent2.split(',')
		P1Plo = int(P1Plo)
		P2Plo = int(P2Plo)
		Acc = options.acc
		
		# Identifications of sites in accordance with the hypothesized parentage and those that are not.
		outfile1 = gzip.open(options.prefix+'_OK.tab.gz','wt')
		outfile2 = gzip.open(options.prefix+'_noOK.tab.gz','wt')
		
		OK = 0
		NoOk = 0
		
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
					dicoVcfPos = {}
					dicoVcfPos[P1name] = header.index(P1name)
					dicoVcfPos[P2name] = header.index(P2name)
					dicoVcfPos[Acc] = header.index(Acc)
					# Recording format position
					FormatPos = header.index('FORMAT')
					# Recording chromosome position
					chrPos = header.index('#CHROM')
					# Recording POS position
					POSPos = header.index('POS')
				
				elif data[0][0] == "#":
					pass
				else:
					chr = data[chrPos]
					pos = data[POSPos]
					FORMAT = data[FormatPos].split(':')
					GTPos = FORMAT.index('GT')
					GTAcc = tuple(data[dicoVcfPos[Acc]].split(':')[GTPos].split('/'))
					GTP1 = data[dicoVcfPos[P1name]].split(':')[GTPos].split('/')
					GTP2 = data[dicoVcfPos[P2name]].split(':')[GTPos].split('/')
					if chr in ChrSet:
						if not('.' in GTAcc) and not('.' in GTP1) and not('.' in GTP2):
							
							IntPos = int(pos)
							
							# Recording position as a studied one
							dico_chr_data[chr][IntPos-1] = 1
							
							# Creating parent gametes
							P1Gametes = combinlisterep(GTP1, P1Plo)
							P2Gametes = combinlisterep(GTP2, P2Plo)
							
							# Creating potential individuals
							dicoInd = set()
							for Gam1 in P1Gametes:
								for Gam2 in P2Gametes:
									Geno = tuple(sorted(Gam1+Gam2))
									dicoInd.add(Geno)
							
							# Testing if individual genotype is possible genotype
							if GTAcc in dicoInd:
								outfile1.write('\t'.join([chr, pos, "1"]))
								outfile1.write('\n')
								OK += 1
							else:
								outfile2.write('\t'.join([chr, pos, "1"]))
								outfile2.write('\n')
								NoOk += 1
								dico_chr[chr][IntPos-1] = 1
		
		outfile1.close()
		outfile2.close()
	
	else:
	
		# Loading parents information
		P1name,P1Plo = options.parent1.split(',')
		P1Plo = int(P1Plo)
		Acc = options.acc
		
		# Identifications of sites in accordance with the hypothesized parentage and those that are not.
		outfile1 = gzip.open(options.prefix+'_OK.tab.gz','wt')
		outfile2 = gzip.open(options.prefix+'_noOK.tab.gz','wt')
		
		OK = 0
		NoOk = 0
		
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
					dicoVcfPos = {}
					dicoVcfPos[P1name] = header.index(P1name)
					dicoVcfPos[Acc] = header.index(Acc)
					# Recording format position
					FormatPos = header.index('FORMAT')
					# Recording chromosome position
					chrPos = header.index('#CHROM')
					# Recording POS position
					POSPos = header.index('POS')
				
				elif data[0][0] == "#":
					pass
				else:
					chr = data[chrPos]
					pos = data[POSPos]
					FORMAT = data[FormatPos].split(':')
					GTPos = FORMAT.index('GT')
					GTAcc = data[dicoVcfPos[Acc]].split(':')[GTPos].split('/')
					GTP1 = data[dicoVcfPos[P1name]].split(':')[GTPos].split('/')
					if chr in ChrSet:
						if not('.' in GTAcc) and not('.' in GTP1):
							
							IntPos = int(pos)
							
							# Recording position as a studied one
							dico_chr_data[chr][IntPos-1] = 1
							
							# Creating potential gametes of parent
							P1Gametes = combinlisterep(GTP1, P1Plo)
							DicP1Gam = set()
							for n in P1Gametes:
								DicP1Gam.add(tuple(sorted(n)))
							
							# Creating potential gametes from child
							ChildGametes = combinliste(GTAcc, P1Plo)
							DicChildGam = set()
							for n in ChildGametes:
								DicChildGam.add(tuple(sorted(n)))
							
							# Testing if parent is a potential donor
							if DicChildGam & DicP1Gam:
								outfile1.write('\t'.join([chr, pos, "1"]))
								outfile1.write('\n')
								OK += 1
							else:
								outfile2.write('\t'.join([chr, pos, "1"]))
								outfile2.write('\n')
								NoOk += 1
								dico_chr[chr][IntPos-1] = 1
		
		outfile1.close()
		outfile2.close()
	
	print("Printing output with nucleotide windows")
	# Calculating statistics along chromosomes with nucleotide windows
	outfile3 = open(options.prefix+'_OK_prop_NUCwin.tab','wt')
	outfile4 = open(options.prefix+'_NoOK_prop_NUCwin.tab','wt')
	listChr = sorted(dico_chr.keys())
	for chr in listChr:
		debut = 0
		fin = WIN
		chr_fin = len(dico_chr[chr])-1
		while debut < chr_fin:
			temp_fin = fin
			if chr_fin < fin:
				temp_fin = chr_fin + 1
			subListNoOk = dico_chr[chr][debut:temp_fin]
			subListPos = dico_chr_data[chr][debut:temp_fin]
			SumPosStudied = sum(subListPos)
			if SumPosStudied > 10:
				Value = (sum(subListNoOk)/float(SumPosStudied))
				outfile3.write('\t'.join([chr, str(debut+1), str(temp_fin), str(1-Value)])+'\n')
				outfile4.write('\t'.join([chr, str(debut+1), str(temp_fin), str(Value)])+'\n')
			debut += WIN
			fin += WIN
			sys.stdout.flush()
	outfile3.close()
	outfile4.close()
	
	print("Printing output with SNP half windows")
	# Calculating statistics along chromosomes with SNP half windows
	outfile3 = open(options.prefix+'_OK_prop_SNPwin.tab','wt')
	outfile4 = open(options.prefix+'_NoOK_prop_SNPwin.tab','wt')
	listChr = sorted(dico_chr.keys())
	for chr in listChr:
		Position = dico_chr_data[chr][:]
		NotOk = dico_chr[chr][:]
		liste_value = []
		liste_pos = []
		First = True
		for i in range(len(Position)):
			if Position[i] == 1:
				liste_pos.append(i+1)
				if NotOk[i] == 1:
					liste_value.append(0)
				else:
					liste_value.append(1)
			if len(liste_value) == WINSNP:
				FinalValue = (sum(liste_value)/float(len(liste_value)))
				FinalPos = liste_pos[HalfWINSNP*-1 - 1]
				if First and HalfWINSNP != 0:
					outfile3.write('\t'.join([chr, str(liste_pos[0]), str(liste_pos[0]), str(FinalValue)])+'\n')
					outfile4.write('\t'.join([chr, str(liste_pos[0]), str(liste_pos[0]), str(1-FinalValue)])+'\n')
					First = False
				outfile3.write('\t'.join([chr, str(FinalPos), str(FinalPos), str(FinalValue)])+'\n')
				outfile4.write('\t'.join([chr, str(FinalPos), str(FinalPos), str(1-FinalValue)])+'\n')
				liste_value = liste_value[HalfWINSNP*-1 - 1:]
				liste_pos = liste_pos[HalfWINSNP*-1 - 1:]
		# Managing last position
		outfile3.write('\t'.join([chr, str(liste_pos[-1]), str(liste_pos[-1]), str(FinalValue)])+'\n')
		outfile4.write('\t'.join([chr, str(liste_pos[-1]), str(liste_pos[-1]), str(1-FinalValue)])+'\n')
	outfile3.close()
	outfile4.close()
		
	print('Total sites inspected:', OK + NoOk)
	print('Sites in accordance with tested origin:', OK, str(OK/float(OK + NoOk)))
	print('Sites not in accordance with tested origin:', NoOk, str(NoOk/float(OK + NoOk)))
	
if __name__ == "__main__": __main__()
