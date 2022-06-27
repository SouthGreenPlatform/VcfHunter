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
import sys
sys.stdout.write('loading modules\n')
import gzip
from operator import itemgetter
import optparse
from os import listdir
from os.path import isfile, join
import re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import math


def __main__():
	#Parse Command Line
	parser = optparse.OptionParser(usage='python %prog [options]\n\nProgram designed by Guillaume MARTIN : guillaume.martin@cirad.fr\n This program convert outputs from PaintArp.py or vcf2linear.1.1.py to inputs for GEMO\n')
	# Wrapper options.
	parser.add_option( '',  '--name',	dest='name',	default=None,			help='The name of the accession. It corresponds to a pattern that will be searched to look for files.')
	parser.add_option( '',  '--dir',	dest='dir',		default=None,			help='The directory that contained the accession.')
	parser.add_option( '',  '--col',	dest='col',		default=None,			help='A color file that will be used to paint accessions')
	parser.add_option( '',  '--size',	dest='size',	default=None,			help='A file containing chromosome size. 2 columns are required: col1 : chromosome name, col2 : chromosome size')
	parser.add_option( '',  '--plo',	dest='plo',		default='2',			help='The ploidy of studied individual. [Default: %default]')
	parser.add_option( '',  '--chro',	dest='chro',	default='all',			help='List of chromosomes to draw separated by ",". If omitted, all chromosomes will be drawn. [Default: %default]')
	parser.add_option( '',  '--prefix',	dest='prefix',	default='ForForGEMO',	help='Prefix for the output files')
	(options, args) = parser.parse_args()

	if options.name == None:
		sys.exit('Please provide an name to --name argument')
	if options.dir == None:
		sys.exit('Please provide an directory to --dir argument')
	if options.size == None:
		sys.exit('Please provide an size file to --size argument')
	
	# Recording color
	DicoFiles = {}
	files = [f for f in listdir(options.dir) if isfile(join(options.dir, f))]
	SetCol = set()
	for f in files:
		fsplit = f.split('_')
		Name = '_'.join(fsplit[:-2])
		chro = fsplit[-2]
		ext = fsplit[-1]
		pattern = re.compile(r'haplo\d.tab')
		if Name == options.name and pattern.search(ext) != None:
			if not (chro in DicoFiles):
				DicoFiles[chro] = []
			DicoFiles[chro].append(f)
	
	for chro in DicoFiles:
		for f in DicoFiles[chro]:
			file = open(join(options.dir, f),'r')
			for line in file:
				data = line.split()
				if data:
					SetCol.add(data[4])
			file.close()
	ListCol = sorted(list(SetCol))
	# print(DicoFiles)
	# print(ListCol)
	
	# Recording potential color
	DicoCol = {}
	IndexColName = {}
	if options.col == None:
		cmap = matplotlib.cm.get_cmap('gist_rainbow')
		for i in range(len(ListCol)):
			if i%2 == 0:
				DicoCol[ListCol[i]] = '#%02x%02x%02x' % cmap((i/2)/(len(ListCol)))[0:3]
				IndexColName[ListCol[i]] = ListCol[i]+'\t#%02x%02x%02x' % cmap((i/2)/(len(ListCol)))[0:3]
			else:
				DicoCol[ListCol[i]] = '#%02x%02x%02x' % cmap((((i-1)/2)+math.ceil(len(ListCol)/2))/(len(ListCol)))[0:3]
				IndexColName[ListCol[i]] = ListCol[i]+'\t#%02x%02x%02x' % cmap((((i-1)/2)+math.ceil(len(ListCol)/2))/(len(ListCol)))[0:3]
	else:
		file = open(options.col)
		# testing the type of color file
		header = file.readline().split()
		if "group" in header and "name" in header and "r" in header and "g" in header and "b" in header:
			rpos = header.index('r')
			gpos = header.index('g')
			bpos = header.index('b')
			grouppos = header.index('group')
			namepos = header.index('name')
			for line in file:
				data = line.split()
				if data:
					DicoCol[data[grouppos]] = '#%02x%02x%02x' % (int(data[rpos]), int(data[gpos]), int(data[bpos]))
					IndexColName[data[grouppos]] = data[namepos]+'\t#%02x%02x%02x' % (int(data[rpos]), int(data[gpos]), int(data[bpos]))
		else:
			data = header
			sur_group = False
			sur_color = False
			while data:
				if data[0] == '[group]':
					sur_group = True
					sur_color = False
				elif data[0] == '[color]':
					sur_group = False
					sur_color = True
				elif sur_color:
					color_list = data[1].split(':')
					dic = {}
					for k in color_list:
						dic[k.split('=')[0]] = float(k.split('=')[1])
					if not('red' in dic):
						sys.exit('There is a problem in the file passed to --col argument: red color proportion is missing')
					if not('green' in dic):
						sys.exit('There is a problem in the file passed to --col argument: green color proportion is missing')
					if not('blue' in dic):
						sys.exit('There is a problem in the file passed to --col argument: blue color proportion is missing')
					if not('alpha' in dic):
						sys.exit('There is a problem in the file passed to --col argument: alpha color proportion is missing')
					DicoCol[data[0]] = '#%02x%02x%02x' % (dic['red'], dic['green'], dic['blue'])
					IndexColName[data[0]] = data[0]+'\t#%02x%02x%02x' % (dic['red'], dic['green'], dic['blue'])
				data = file.readline().split()
		file.close()
		
		for col in ListCol:
			if not (col in DicoCol):
				print(col+' was not found in color file, an arbitrary grey color was attributed')
				DicoCol[col] = '#%02x%02x%02x' % (189, 189, 189)
				IndexColName[col] = col+'\t#%02x%02x%02x' % (189, 189, 189)
	
	# Printing color file
	outfile = open(options.prefix+'_color.tab','w')
	outfile.write("group\tname\thex\n")
	for col in IndexColName:
		outfile.write(col+'\t'+IndexColName[col]+'\n')
	outfile.close()
	
	# Generating files for ideogram
	outfile = open(options.prefix+'_ideoProb.tab','w')
	outfile.write('\t'.join(['chr', 'haplotype', 'start', 'end', 'ancestral_group']))
	outfile.write('\n')
	if options.chro == 'all':
		listChr = sorted(list(DicoFiles.keys()))
	else:
		listChr = options.chro.split(',')
	# dicoChrId = {}
	for ichr in range(len(listChr)):
		chro = listChr[ichr]
		# dicoChrId[chro] = ichr
		for i in range(len(DicoFiles[chro])):
			file = open(join(options.dir, DicoFiles[chro][i]), 'r')
			for line in file:
				data = line.split()
				if data:
					outfile.write(chro+'\t'+str(i)+'\t'+data[2]+'\t'+data[3]+'\t'+data[4]+'\n')
	outfile.close()
	
	# Reformatting chromosome file
	file = open(options.size,'r')
	dicoChr = {}
	for line in file:
		data = line.split()
		if data:
			dicoChr[data[0]] = data[1]
	file.close()
	outfile = open(options.prefix+'_chrom.tab','w')
	outfile.write("\t".join(["chr", "len","centromereInf", "centromereSup", "label"])+"\n")
	for chro in listChr:
		mot = ''
		for i in range(len(DicoFiles[chro])):
			mot = mot + chr(i+65)
		outfile.write(chro+'\t'+dicoChr[chro]+'\t'+str(int(int(dicoChr[chro])/2))+'\t'+str(int(int(dicoChr[chro])/2)+2)+'\t'+mot+'\n')
	outfile.close()
	
	# Reformating curve file
	CurveFile = options.dir+'/'+options.name+'_win_ratio.tab.gz'
	if CurveFile == None:
		pass
	else:
		if CurveFile[-3:] == '.gz':
			file = gzip.open(CurveFile,'rt')
		else:
			file = open(CurveFile, 'r')
		outfile = open(options.prefix+'_curve.tab','w')
		header = file.readline().split()
		CHRPOS = header.index("CHR")
		POSPOS = header.index("POS")
		STAPOS = header.index("Start")
		ENDPOS = header.index("End")
		GRPVAL = header[ENDPOS+1:]
		AfterHeader = True
		for line in file:
			data = line.split()
			if data:
				if AfterHeader:
					outfile.write('\t'.join(['chr','start','end']+header[ENDPOS+1:]))
					outfile.write('\n')
					cchr = data[CHRPOS]
					csta = int(data[STAPOS])
					cpos = int(data[POSPOS])
					cend = int(data[ENDPOS])
					cval = data[ENDPOS+1:]
					ValueFinal = []
					for n in cval:
						if n == "NA":
							ValueFinal.append('0.0')
						else:
							ValueFinal.append(n)
					AfterHeader = False
				elif cchr == data[CHRPOS]:
					outfile.write('\t'.join([cchr,str(csta),str(int((cpos+int(data[POSPOS]))/2.0))]+ValueFinal))
					outfile.write('\n')
					cchr = data[CHRPOS]
					csta = int((cpos+int(data[POSPOS]))/2.0)
					cpos = int(data[POSPOS])
					cend = int(data[ENDPOS])
					cval = data[ENDPOS+1:]
					ValueFinal = []
					for n in cval:
						if n == "NA":
							ValueFinal.append('0.0')
						else:
							ValueFinal.append(n)
				else:
					outfile.write('\t'.join([cchr,str(csta),str(cend)]+ValueFinal))
					outfile.write('\n')
					cchr = data[CHRPOS]
					csta = int(data[STAPOS])
					cpos = int(data[POSPOS])
					cend = int(data[ENDPOS])
					cval = data[ENDPOS+1:]
					ValueFinal = []
					for n in cval:
						if n == "NA":
							ValueFinal.append('0.0')
						else:
							ValueFinal.append(n)
		
		outfile.write('\t'.join([cchr,str(csta),str(cend)]+ValueFinal))
		outfile.write('\n')
		outfile.close()
		file.close()
	
	# Calculating blocks according to GEMO algorithm
	
	file = open(options.prefix+'_curve.tab')
	header = file.readline().split()
	CHRPOS = header.index("chr")
	STAPOS = header.index("start")
	ENDPOS = header.index("end")
	GRPVAL = header[ENDPOS+1:]
	dicoGpPos = {}
	for n in range(len(GRPVAL)):
		dicoGpPos[n] = GRPVAL[n]
	Threshold = 1/int(options.plo)
	Ploidy = int(options.plo)
	LastNonUnHaplo = ['un']*Ploidy
	dicoPierrick = {}
	IDPierrick = 0
	listchr = set()
	chro = ''
	AddTheStart = False
	for line in file:
		data = line.split()
		if data:
			if chro != data[CHRPOS]:# This is the begining of a new chromosome
				if chro == '': # This is the first line
					pass
				else: # This is the junction between two chromosomes so we should add the end of former chromosome
					if end < int(dicoChr[chro]): # The end of the data are before the end of the chromosome
						dicoPierrick[IDPierrick] = {}
						for n in range(Ploidy): # We fill this start region with missing information
							listchr.add(chro)
							dicoPierrick[IDPierrick][n] = {}
							dicoPierrick[IDPierrick][n][0] = chro
							dicoPierrick[IDPierrick][n][1] = n
							dicoPierrick[IDPierrick][n][2] = end
							dicoPierrick[IDPierrick][n][3] = int(dicoChr[chro])
							dicoPierrick[IDPierrick][n][4] = haplo[n]
							# print('Avant  i:'+str(IDPierrick)+' j:'+str(n)+' --- '+str(dicoPierrick[IDPierrick][n][0])+','+str(dicoPierrick[IDPierrick][n][1])+','+str(dicoPierrick[IDPierrick][n][2])+','+str(dicoPierrick[IDPierrick][n][3])+','+str(DicoCol[dicoPierrick[IDPierrick][n][4]]))
						IDPierrick += 1
				if int(data[STAPOS]) > 0: # The start of data for the new chromosome is not the start of the chromosome
					AddTheStart = True
			
			# Standard treatment of line
			chro = data[CHRPOS]
			start = int(data[STAPOS])
			end = int(data[ENDPOS])
			haplo = []
			Values = data[ENDPOS+1:]
			for Gp in range(len(GRPVAL)):
				for n in range(min(int(float(Values[Gp])/Threshold),Ploidy)):
					haplo.append(dicoGpPos[Gp])
			# print(chro, start, end, haplo)
			## At this step  we have the haplotype deduced from curves in a list we should now decide if the number of deduced haplotypes is ok or not
			if len(haplo) <= Ploidy: # The number of deduced haplotypes is compatible with the ploidy. 
				for i in range(Ploidy-len(haplo)): # Adding unknown haplotype if necessary according to ploidy
					haplo.append("un")
					# haplo.insert(0,"un")
			else : # The number of deduced haplotypes is not compatible with the ploidy. Haplotypes in this regions are this converted to unknown
				haplo = ["un"]*Ploidy
			# print('toto', chro, start ,haplo)
			
			if AddTheStart: # extending the start region of the chromosome with the following region
				dicoPierrick[IDPierrick] = {}
				for n in range(Ploidy): # We fill this start region with missing information
					listchr.add(chro)
					dicoPierrick[IDPierrick][n] = {}
					dicoPierrick[IDPierrick][n][0] = chro
					dicoPierrick[IDPierrick][n][1] = n
					dicoPierrick[IDPierrick][n][2] = 0
					dicoPierrick[IDPierrick][n][3] = int(data[STAPOS])
					dicoPierrick[IDPierrick][n][4] = haplo[n]
					# print('Avant  i:'+str(IDPierrick)+' j:'+str(n)+' --- '+str(dicoPierrick[IDPierrick][n][0])+','+str(dicoPierrick[IDPierrick][n][1])+','+str(dicoPierrick[IDPierrick][n][2])+','+str(dicoPierrick[IDPierrick][n][3])+','+str(DicoCol[dicoPierrick[IDPierrick][n][4]]))
				IDPierrick += 1
				AddTheStart = False
				
			
			dicoPierrick[IDPierrick] = {}
			for n in range(Ploidy):
				listchr.add(chro)
				dicoPierrick[IDPierrick][n] = {}
				dicoPierrick[IDPierrick][n][0] = chro
				dicoPierrick[IDPierrick][n][1] = n
				dicoPierrick[IDPierrick][n][2] = start
				dicoPierrick[IDPierrick][n][3] = end
				dicoPierrick[IDPierrick][n][4] = haplo[n]
				# print('Avant  i:'+str(IDPierrick)+' j:'+str(n)+' --- '+str(dicoPierrick[IDPierrick][n][0])+','+str(dicoPierrick[IDPierrick][n][1])+','+str(dicoPierrick[IDPierrick][n][2])+','+str(dicoPierrick[IDPierrick][n][3])+','+str(DicoCol[dicoPierrick[IDPierrick][n][4]]))
			
			IDPierrick += 1
	
	# This is for the last line
	if end < int(dicoChr[chro]): # The end of the data are before the end of the chromosome
		dicoPierrick[IDPierrick] = {}
		for n in range(Ploidy): # We fill this start region with missing information
			listchr.add(chro)
			dicoPierrick[IDPierrick][n] = {}
			dicoPierrick[IDPierrick][n][0] = chro
			dicoPierrick[IDPierrick][n][1] = n
			dicoPierrick[IDPierrick][n][2] = end
			dicoPierrick[IDPierrick][n][3] = int(dicoChr[chro])
			dicoPierrick[IDPierrick][n][4] = haplo[n]
			# print('Avant  i:'+str(IDPierrick)+' j:'+str(n)+' --- '+str(dicoPierrick[IDPierrick][n][0])+','+str(dicoPierrick[IDPierrick][n][1])+','+str(dicoPierrick[IDPierrick][n][2])+','+str(dicoPierrick[IDPierrick][n][3])+','+str(DicoCol[dicoPierrick[IDPierrick][n][4]]))
		IDPierrick += 1
	
	file.close()
	
	## Partie de Pierrick recodee (GroupByColor)
	dicoGroupHaplo = {}
	for i in dicoPierrick:
		for j in dicoPierrick[i]:
			group = dicoPierrick[i][j][4]
			haplo = dicoPierrick[i][j][1]
			chro = dicoPierrick[i][j][0]
			start = dicoPierrick[i][j][2]
			if not(group in dicoGroupHaplo and group != 'un'):
				dicoGroupHaplo[group] = haplo
			
			elif haplo != dicoGroupHaplo[group]:# le groupe n'est pas sur le bon haplotype
				GroupInBlock = set()
				for k in range(Ploidy):
					GroupInBlock.add(dicoPierrick[i][k][4])
				if len(GroupInBlock) > 1:
					for k in range(Ploidy):
						if dicoPierrick[i][k][1] == dicoGroupHaplo[group]:
							dicoPierrick[i][k][1] = haplo
					dicoPierrick[i][j][1] = dicoGroupHaplo[group]
			else:
				pass
				# print(group in dicoGroupHaplo, group)
			GroupInBlock = set()
	
	## Ordering by chromosomes, haplotype and position
	ListToSort = []
	for i in dicoPierrick:
		for j in dicoPierrick[i]:
			chro = dicoPierrick[i][j][0]
			haplo = dicoPierrick[i][j][1]
			start = dicoPierrick[i][j][2]
			end = dicoPierrick[i][j][3]
			group = dicoPierrick[i][j][4]
			ListToSort.append((chro, haplo, start, end, group))
	ListToSort = sorted(ListToSort, key=itemgetter(0,1,2))
	
	
	
	## Block for output file header
	outfile = open(options.prefix+'_ideo.tab','w')
	outfile.write('\t'.join(['chr', 'haplotype', 'start', 'end', 'ancestral_group']))
	outfile.write('\n')
	## end of block for output file header now printing output
	for n in ListToSort:
		outfile.write('\t'.join(list(map(str,n))))
		outfile.write('\n')
	outfile.close()
	file.close()
	
if __name__ == "__main__": __main__()	
