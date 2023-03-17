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
import os
import sys
import gzip
import math
import optparse
import itertools
import numpy as np
from operator import itemgetter

def NormDens(n,p,q,x):
	mu = n*p
	si = math.sqrt(n*p*q)
	v = (x-mu)/si
	val = math.exp(-0.5*v*v)/(si*math.sqrt(2*math.pi))
	return val

def __main__():
	#Parse Command Line
	parser = optparse.OptionParser(usage="python %prog [options]\n\nProgram designed by Guillaume MARTIN (guillaume.martin@cirad.fr)"
	"\n\nThis program use allele specific ratios calculated in one accession to"
	"\ncharacterize the ancestry along chromosomes of studied accession.")
	# Wrapper options.
	parser.add_option( '-r', '--ratio',		dest='ratio',		default=None,		help='The ratio file. Tabulated file with 6 columns with headers. Col1: chr, col2: pos, col3: allele, col4: obs_ratio, col5: exp_ratio, col6: grp')
	parser.add_option( '-c', '--color',		dest='color',		default=None,		help='Color file name. Tabulated file with 5 columns with header. Col1: group, col2: name, col3: r, col4: g, col5: b')
	parser.add_option( '-p', '--ploidy',	dest='ploidy',		default='2',		help='Ploidy of studied accession')
	parser.add_option( '-w', '--win',		dest='win',			default='300',		help='Half window size of markers')
	parser.add_option( '-O', '--overlap',	dest='overlap',		default=None,		help='Overlap between two windows (number of SNP positions). Of not filled, the overlap between two windows of size n will be of n-1. ')
	parser.add_option( '-N', '--noise',		dest='noise',		default='0.05',		help='Maximal mean read proportion threshold in which absence of haplotype probability is put to 1')
	parser.add_option( '-n', '--threshold',	dest='threshold',	default='0.05',		help='Minimal mean read proportion threshold relative to the expected proportion in which haplotype probability is put to 1')
	parser.add_option( '-s', '--size',		dest='size',		default=None,		help='A file containing chromosome size. 2 columns are required: col1 : chromosome name, col2 : chromosome size')
	parser.add_option( '-m', '--MinAll',	dest='MinAll',		default='10',		help='Minimal allele number of an origin to keep this region for prediction for the origin')
	parser.add_option( '-a', '--acc',		dest='acc',			default='Unknown',	help='Accession name')
	parser.add_option( '-o', '--out',		dest='out',			default='Test',		help='Prefix for output files')
	(options, args) = parser.parse_args()
	
	
	
	MinAll = int(options.MinAll)
	pNOISE = float(options.noise)
	qNoise = 1-pNOISE
	THRESHOLD = float(options.threshold)
	PLO = int(options.ploidy)
	HalfWIN = int(options.win)
	TotW = HalfWIN*2+1
	
	if options.overlap == None:
		PosToRmNum = TotW-(HalfWIN*2)
	else:
		PosToRmNum = TotW-int(options.overlap)
	
	UnSet = set()
	UnSet.add('un')
	
	# Loading chromosome size information
	DiChSi = {}
	file = open(options.size)
	for line in file:
		data = line.split()
		if data:
			DiChSi[data[0]] = data[1]
	file.close()
	
	# Preparing ploidy calculation
	DiPlo = {}
	for plo in range(1, PLO+1):
		p = plo/PLO
		q = 1-p
		DiPlo[plo] = (p-THRESHOLD,q+THRESHOLD)
		
	# Recording color information
	DiCol = {}
	file  = open(options.color)
	header = file.readline().split()
	for line in file:
		data = line.split()
		if data:
			DiCol[data[header.index('group')]] = [data[header.index('name')], int(data[header.index('r')])/255, int(data[header.index('g')])/255, int(data[header.index('b')])/255]
	file.close()
	
	# Calculating possible combination
	FinalComb = []
	FinalComb.append(tuple())
	for i in range(1, PLO+1):
		FinalComb += list(itertools.combinations_with_replacement(DiCol.keys(), i))
	
	# Validating file
	if options.ratio[-3:] == '.gz':
		file = gzip.open(options.ratio,'rt')
	else:
		file = open(options.ratio)
	header = file.readline().split()
	chrP = header.index('chr')
	PosP = header.index('pos')
	AllP = header.index('allele')
	ObsP = header.index('obs_ratio')
	ExpP = header.index('exp_ratio')
	GrpP = header.index('grp')
	fistLine = file.readline().split()
	PrecChr = fistLine[chrP]
	PrecPos = int(fistLine[PosP])
	ChrDone = set()
	ChrDone.add(PrecChr)
	for line in file:
		data = line.split()
		Chr = data[chrP]
		Pos = int(data[PosP])
		if Chr == PrecChr:
			if PrecPos > Pos:
				sys.exit('Oups, the file is not chromosome and  position sorted. The program could not work.\n')
		elif Chr in ChrDone:
			sys.exit('Oups, the file is not chromosome and  position sorted. The program could not work.\n')
	LastLine = data[:]
	sys.stdout.write('File format seems OK. continue analysis...\n')
	
	# Estimating probabilities
	dico = {}
	dicoGP = {}
	for gp in DiCol:
		dicoGP[gp] = {}
		dicoGP[gp]['Obs'] = []
		dicoGP[gp]['Exp'] = []
	dico['Pos'] = []
	dico['gp'] = []
	
	outratio = gzip.open(options.out+'_win_ratio.tab.gz', 'wt')
	listGP = list(sorted(list(dicoGP.keys())))
	outratio.write('\t'.join(["CHR", "POS", "Start","End"]+listGP))
	outratio.write('\n')
	
	if options.ratio[-3:] == '.gz':
		file = gzip.open(options.ratio,'rt')
	else:
		file = open(options.ratio)
	
	header = file.readline().split()
	chrP = header.index('chr')
	PosP = header.index('pos')
	AllP = header.index('allele')
	ObsP = header.index('obs_ratio')
	ExpP = header.index('exp_ratio')
	GrpP = header.index('grp')
	CHR = ''
	Start = True
	i = 0
	Counter = 0
	LastEndBlock = 0                                                                                                          # for management of last block and compatibility with GEMO
	for line in file:
		i += 1
		data = line.split()
		Chr = data[chrP]
		Pos = int(data[PosP])
		Obs = float(data[ObsP])
		Exp = float(data[ExpP])
		Grp = data[GrpP]
		if not(Pos) in dico['Pos'] or LastLine == data:                                                                       # adding "or LastLine == data"
			NbPos = len(set(dico['Pos']))
			if NbPos == TotW:
				
				PotentialOK = False
				if LastLine == data:                                                                                          # This as been deplaced before if CHR == Chr and Counter >= PosToRmNum
					PotentialOK = True                                                                                        # This as been deplaced before if CHR == Chr and Counter >= PosToRmNum
					# Adding new value to list                                                                                # This as been deplaced before if CHR == Chr and Counter >= PosToRmNum
					dico['Pos'].append(Pos)                                                                                   # This as been deplaced before if CHR == Chr and Counter >= PosToRmNum
					dico['gp'].append(Grp)                                                                                    # This as been deplaced before if CHR == Chr and Counter >= PosToRmNum
					for gp in DiCol:                                                                                          # This as been deplaced before if CHR == Chr and Counter >= PosToRmNum
						if gp == Grp:                                                                                         # This as been deplaced before if CHR == Chr and Counter >= PosToRmNum
							dicoGP[gp]['Obs'].append(Obs)                                                                     # This as been deplaced before if CHR == Chr and Counter >= PosToRmNum
							dicoGP[gp]['Exp'].append(Exp)                                                                     # This as been deplaced before if CHR == Chr and Counter >= PosToRmNum
						else:                                                                                                 # This as been deplaced before if CHR == Chr and Counter >= PosToRmNum
							dicoGP[gp]['Obs'].append(0)                                                                       # This as been deplaced before if CHR == Chr and Counter >= PosToRmNum
							dicoGP[gp]['Exp'].append(0)                                                                       # This as been deplaced before if CHR == Chr and Counter >= PosToRmNum
					print('This is the last line, if this message appear more than once, this means that there is a problem') # This as been deplaced before if CHR == Chr and Counter >= PosToRmNum
				elif CHR == Chr and Counter >= PosToRmNum:
					PotentialOK = True
				elif CHR != Chr:
					PotentialOK = True
				else:
					pass
				
				if PotentialOK:
					Counter = 0
					diCal = {}
					NoEnAll = set()
					MedianPos = int(np.median(np.array(list(set(dico['Pos'])))))
					LIST2Print = [CHR, MedianPos, min(dico['Pos']), max(dico['Pos'])]
					
					for gp in listGP:
						###### A piece of verification
						if len(dico['Pos']) != len(dico['gp']):
							sys.exit('Oups there is a bug or a problem in files')
						if len(dico['Pos']) != len(dicoGP[gp]['Obs']):
							sys.exit('Oups there is a bug or a problem in files')
						if len(dico['Pos']) != len(dicoGP[gp]['Exp']):
							sys.exit('Oups there is a bug or a problem in files')
						###### Verification completed
						
						
						NbGp = dico['gp'].count(gp)
						# Managing insufficient specific markers
						if NbGp < MinAll:
							NoEnAll.add(gp)
							LIST2Print.append("NA")
						else:
							NormRatio = min(1, sum(dicoGP[gp]['Obs'])/sum(dicoGP[gp]['Exp']))
							LIST2Print.append(NormRatio)
							diCal[gp] = {}
							# Absence of haplotype probability calculation
							if NormRatio <= pNOISE:
								diCal[gp][0] = 1
							else:
								diCal[gp][0] = NormDens(NbGp,pNOISE,qNoise,NbGp*NormRatio)/NormDens(NbGp,pNOISE,qNoise,NbGp*pNOISE)
							###################################
							
							# Haplotype dose probability calculation
							for plo in DiPlo:
								if NormRatio >= DiPlo[plo][0]:
									diCal[gp][plo] = 1
								else:
									diCal[gp][plo] = NormDens(NbGp,DiPlo[plo][0],DiPlo[plo][1],NbGp*NormRatio)/NormDens(NbGp,DiPlo[plo][0],DiPlo[plo][1],NbGp*DiPlo[plo][0])
							###################################
					if CHR == Chr and LastLine != data:                                                                       # for management of last block and compatibility with GEMO
						outratio.write('\t'.join(map(str, LIST2Print)))                                                       # for management of last block and compatibility with GEMO
						outratio.write('\n')                                                                                  # for management of last block and compatibility with GEMO
						LastEndBlock = LIST2Print[3]+1                                                                        # for management of last block and compatibility with GEMO
					else:                                                                                                     # for management of last block and compatibility with GEMO
						LIST2Print[2] = LastEndBlock                                                                          # for management of last block and compatibility with GEMO
						outratio.write('\t'.join(map(str, LIST2Print)))                                                       # for management of last block and compatibility with GEMO
						outratio.write('\n')                                                                                  # for management of last block and compatibility with GEMO
						print("In the win_ratio.tab.gz file, the last windows of chromosome ",CHR," has been adjusted so that it is compatible with GEMO.") # for management of last block and compatibility with GEMO
					
					if len(NoEnAll) > 0:
						print("For region",  CHR, MedianPos, min(dico['Pos']), max(dico['Pos']), 'there are not enough alleles for origin(s)', NoEnAll)
					
					# Calculating probabilities
					ListHaploProb = []
					listProb = []
					for comb in FinalComb:
						OKForProba = True
						for GPInComb in comb:
							if GPInComb in NoEnAll:
								OKForProba = False
						if OKForProba:
							Prob = 1
							for gp in DiCol:
								if not(gp in NoEnAll):
									NbGpComb = comb.count(gp)
									Prob = Prob*diCal[gp][NbGpComb]
							ListHaploProb.append((comb, Prob))
							listProb.append(Prob)
						else:
							ListHaploProb.append((comb, 0))
							listProb.append(0)
					
					ListHaploProb = sorted(ListHaploProb, key=itemgetter(1))
					
					BestCouple, BestVal = ListHaploProb[-1]
					
					# Managing case where several haplotypes compositions have best probability values
					if listProb.count(BestVal) > 1:
						SeveralBest = ListHaploProb[listProb.count(BestVal)*(-1):]
						BestCouple = ()
						BESTCOUPLE = []
						OK = True
						
						# Looking if there is a uniq combination with a maximum of haplotypes attributed.
						BestHaplo = []
						for best in range(len(SeveralBest)):
							BESTCOUPLE.append((SeveralBest[best][0], len(SeveralBest[best][0])))
						ListBestHaploNum = sorted(BESTCOUPLE, key=itemgetter(1))
						BestCouple, BestValNum = ListBestHaploNum[-1]
						
						# Now testing if all other best haplotypes are included in final best haplotype, if it is not the case all haplotypes are converted to missing
						for couple in ListBestHaploNum:
							BESTCOUPLENUM = list(BestCouple)
							for CPL in couple[0]:
								if CPL in BESTCOUPLENUM:
									BESTCOUPLENUM.remove(CPL)
								else:
									OK = False
						if OK:
							pass
						else:
							print('NO', BestVal, ListHaploProb[listProb.count(BestVal)*(-1):], listProb.count(BestVal))
							BestCouple = tuple(['un']*PLO)
					
					###########################
					
					if Start:
						StartPos = min(dico['Pos'])
						# Preparing haplotype construction
						DicoHaplo = {}
						ListHaplo = [[1]*PLO, [StartPos-1]*PLO, ['un']*PLO]
						for plo in range(PLO):
							outfile = open(options.out+'_'+CHR+'_haplo'+str(plo+1)+'.tab', 'w')
							outfile.close()
							DicoHaplo[plo] = {}
							for gp in DiCol:
								DicoHaplo[plo][gp] = 0
							DicoHaplo[plo]['un'] = 0
						###################################
					
					# Attributing origin to haplotype
					ListTooAttribute = list(BestCouple)+['un']*(PLO-len(BestCouple))
					if len(ListTooAttribute) != PLO:
						sys.exit('Bug1\n')
					HaploTaken = set()
					ToRemove = []
					# Attributing origin already present
					for n in ListTooAttribute:
						if n != 'un':
							Possiblehaplo = []
							for plo in range(PLO):
								if not(plo in HaploTaken):
									if ListHaplo[2][plo] == n:
										Possiblehaplo.append((plo, DicoHaplo[plo][n]))
							# Selecting best position
							if Possiblehaplo:
								Possiblehaplo = sorted(Possiblehaplo, key=itemgetter(1))
								HaploTaken.add(Possiblehaplo[-1][0])
								DicoHaplo[Possiblehaplo[-1][0]][n] += 1
								ToRemove.append(n)
								ListHaplo[1][Possiblehaplo[-1][0]] = MedianPos
					# Removing attributed origin
					for n in ToRemove:
						ListTooAttribute.remove(n)
					
					# Attributing other origin
					while set(ListTooAttribute) - UnSet:
						Possiblehaplo = []
						for n in ListTooAttribute:
							for plo in range(PLO):
								if not(plo in HaploTaken):
									Possiblehaplo.append((plo, DicoHaplo[plo][n], n))
						Possiblehaplo = sorted(Possiblehaplo, key=itemgetter(1))
						plo = Possiblehaplo[-1][0]
						gp = Possiblehaplo[-1][2]
						HaploTaken.add(plo)
						ListTooAttribute.remove(gp)
						DicoHaplo[plo][gp] += 1
						outfile = open(options.out+'_'+CHR+'_haplo'+str(plo+1)+'.tab', 'a')
						outfile.write('\t'.join([options.acc, CHR, str(ListHaplo[0][plo]), str(ListHaplo[1][plo]), ListHaplo[2][plo]]))	
						outfile.write('\n')	
						outfile.close()
						if Start:
							ListHaplo[0][plo] = ListHaplo[1][plo]+1
						else:
							ListHaplo[0][plo] = MedianPos
						ListHaplo[2][plo] = gp
						ListHaplo[1][plo] = MedianPos
					
					# Attributing "un" origin
					if ListTooAttribute:
						# print(ListTooAttribute)
						for plo in range(PLO):
							if not(plo in HaploTaken):
								if ListHaplo[2][plo] == "un":
									ListHaplo[1][plo] = MedianPos
								else:
									outfile = open(options.out+'_'+CHR+'_haplo'+str(plo+1)+'.tab', 'a')
									outfile.write('\t'.join([options.acc, CHR, str(ListHaplo[0][plo]), str(ListHaplo[1][plo]), ListHaplo[2][plo]]))	
									outfile.write('\n')	
									outfile.close()
									ListHaplo[1][plo] = MedianPos
									ListHaplo[0][plo] = MedianPos
									ListHaplo[2][plo] = "un"
								HaploTaken.add(plo)
					
					Start = False
					
				if CHR != Chr:
					print('End of' , CHR, 'processing. Now working on', Chr)
					Start = True
					EndPos = max(dico['Pos'])
					for plo in range(PLO):
						outfile = open(options.out+'_'+CHR+'_haplo'+str(plo+1)+'.tab', 'a')
						outfile.write('\t'.join([options.acc, CHR, str(ListHaplo[0][plo]), str(EndPos), ListHaplo[2][plo]]))	
						outfile.write('\n')
						outfile.write('\t'.join([options.acc, CHR, str(EndPos+1), str(DiChSi[CHR]), 'un']))	
						outfile.write('\n')	
						outfile.close()
					dico['Pos'] = []
					dico['gp'] = []
					for gp in DiCol:
						dicoGP[gp] = {}
						dicoGP[gp]['Obs'] = []
						dicoGP[gp]['Exp'] = []
				
				elif LastLine == data:
					print('End of' , CHR, 'processing. Writing its last line')
					Start = True
					EndPos = max(dico['Pos'])
					for plo in range(PLO):
						outfile = open(options.out+'_'+CHR+'_haplo'+str(plo+1)+'.tab', 'a')
						outfile.write('\t'.join([options.acc, CHR, str(ListHaplo[0][plo]), str(EndPos), ListHaplo[2][plo]]))	
						outfile.write('\n')
						outfile.write('\t'.join([options.acc, CHR, str(EndPos+1), str(DiChSi[CHR]), 'un']))	
						outfile.write('\n')	
						outfile.close()
					dico['Pos'] = []
					dico['gp'] = []
					for gp in DiCol:
						dicoGP[gp] = {}
						dicoGP[gp]['Obs'] = []
						dicoGP[gp]['Exp'] = []
					
				
				# Removing first position of the list #################
				if dico['Pos']:
					minPos = min(dico['Pos'])
					while minPos in dico['Pos']:
						PosToRm = dico['Pos'].index(minPos)
						for lts in dicoGP:
							del dicoGP[lts]['Obs'][PosToRm]
							del dicoGP[lts]['Exp'][PosToRm]
						del dico['Pos'][PosToRm]
						del dico['gp'][PosToRm]
			Counter += 1
		# Adding new value to list
		dico['Pos'].append(Pos)
		dico['gp'].append(Grp)
		for gp in DiCol:
			if gp == Grp:
				dicoGP[gp]['Obs'].append(Obs)
				dicoGP[gp]['Exp'].append(Exp)
			else:
				dicoGP[gp]['Obs'].append(0)
				dicoGP[gp]['Exp'].append(0)
		# Recording chromosome name
		CHR = Chr
	
	file.close()
	
	# Printing last line of last chromosome
	# Start = True
	# EndPos = max(dico['Pos'])
	# for plo in range(PLO):
		# outfile = open(options.out+'_'+CHR+'_haplo'+str(plo+1)+'.tab', 'a')
		# outfile.write('\t'.join([options.acc, CHR, str(ListHaplo[0][plo]), str(EndPos), ListHaplo[2][plo]]))	
		# outfile.write('\n')
		# outfile.write('\t'.join([options.acc, CHR, str(EndPos+1), str(DiChSi[CHR]), 'un']))	
		# outfile.write('\n')	
		# outfile.close()
	
	
if __name__ == "__main__": __main__()
