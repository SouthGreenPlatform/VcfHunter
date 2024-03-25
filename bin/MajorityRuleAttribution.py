#!/usr/bin/env python
#
#  Copyright 2024 CIRAD
#  Guillaume Martin
#
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
# -*- coding: utf-8 -*-

import optparse
import gzip

def __main__():
	#Parse Command Line
	parser = optparse.OptionParser()
	# Wrapper options. 
	parser.add_option( '-f', '--files',			dest='files',		default=None,		help='List of density files separated by ",". These files should be formated as follows: ID1:File1,ID2:File2,...')
	parser.add_option( '-c', '--color',			dest='color',		default=None,		help='A color file in tabulated format (in RGB). Col1: RefxName, Col2: Red, Col3: Green, Col4: Blue')
	parser.add_option( '-m', '--mini_total',	dest='mini_total',	default=10,			help='A minimal value (excluded) for the sum of all observed values in the windows for all origin to attribute an origin. [default: 10]')
	parser.add_option( '-M', '--maxi_total',	dest='maxi_total',	default=None,		help='A maximal value (excluded) for the sum of all observed values in the windows for all origin to attribute an origin. No default value, if this value is omitted, the filter is not applicated.')
	parser.add_option( '-t', '--threshold',		dest='threshold',	default=0.6,		help='Threshold value representing the minimal proportion (excluded) of the best origin relative to other origins at the sudied windows. [default: 0.6]')
	parser.add_option( '-o', '--outfile',		dest='outfile',		default='Stacked',	help='Output file name')
	(options, args) = parser.parse_args()
	
	# Obtaining variables
	MINI=float(options.mini_total)
	THRESHOLD=float(options.threshold)
	if options.maxi_total != None:
		MAXI=float(options.maxi_total)
	else:
		MAXI=False
	
	# Opening files and identification of origins
	FilesList = options.files.split(',')
	dicoFile = {}
	OUTFILE = open(options.outfile, 'w')
	for File in FilesList:
		ID, filename = File.split(":")
		dicoFile[ID] = gzip.open(filename, 'rt')
	Origin=sorted(dicoFile.keys())
	
	# Obtaining colors if any
	DicoCol = {}
	if options.color == None:
		cmap = matplotlib.cm.get_cmap('gist_rainbow')
		for i in range(len(Origin)):
			if i%2 == 0:
				DicoCol[Origin[i]] = cmap((i/2)/(len(Origin)))
			else:
				DicoCol[Origin[i]] = cmap((((i-1)/2)+math.ceil(len(Origin)/2))/(len(Origin)))
	else:
		file = open(options.color)
		for line in file:
			data = line.split()
			if data:
				DicoCol[data[0]] = ','.join((data[1], data[2], data[3]))
		file.close()
		for ref in Origin:
			if not (ref in DicoCol):
				print(ref+' was not found in color file, an arbitrary color was attributed')
				DicoCol[ref] = ','.join(("238", "16", "16"))
	
	# Attributing an origin
	i=0
	for line in dicoFile[Origin[0]]:
		chr1,start1,end1,value1=line.split()
		value1 = float(value1)
		total = float(value1)
		listValue=[float(value1)]
		for ori in Origin[1:]:
			chr2,start2,end2,value2=dicoFile[ori].readline().split()
			if chr1!=chr2 or start1!=start2 or end1!=end2:
				sys.exit("There is a bug in the files")
			value2 = float(value2)
			total += value2
			listValue.append(value2)
		
		if total > MINI:
			if MAXI:
				if total < MAXI:
					MaxVal=max(listValue)
					if listValue.count(MaxVal) == 1:
						if MaxVal/total > THRESHOLD:
							i += 1
							BestOri=Origin[listValue.index(MaxVal)]
							OUTFILE.write("Mark"+str(i)+"\t"+chr1+"\t"+start1+"\t"+end1+"\t"+BestOri+"\t"+DicoCol[BestOri]+"\n")
			else:
				MaxVal=max(listValue)
				if listValue.count(MaxVal) == 1:
					if MaxVal/total > THRESHOLD:
						i += 1
						BestOri=Origin[listValue.index(MaxVal)]
						OUTFILE.write("Mark"+str(i)+"\t"+chr1+"\t"+start1+"\t"+end1+"\t"+BestOri+"\t"+DicoCol[BestOri]+"\n")
	
	for ori in Origin:
		dicoFile[ori].close()
	OUTFILE.close()
	
if __name__ == "__main__": __main__()


