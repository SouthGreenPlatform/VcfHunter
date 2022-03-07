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
import optparse
import os

sys.stdout.write('modules loaded\n')

def __main__():
	#Parse Command Line
	parser = optparse.OptionParser(usage="python %prog [options]\n\nProgram designed by Guillaume MARTIN : guillaume.martin@cirad.fr")
	# Wrapper options. 
	parser.add_option( '',	'--tree',	dest='tree',	default=None,			help='The tree file obtained from PhyMl. [Default: %default]')
	parser.add_option( '',	'--fasta',	dest='fasta',	default=None,			help='The corresponding fasta file. [Default: %default]')
	parser.add_option( '',	'--color',	dest='color',	default=None,			help='A color file with column 1 = accession name, column 2 = color in hexadecimal. [Default: %default]')
	parser.add_option( '',	'--layout',	dest='layout',	default='RECTILINEAR',	help='The layout of the tree. Possible values: "RECTILINEAR", "RADIAL" or "POLAR" [Default: %default]')
	parser.add_option( '',	'--fsize',	dest='fsize',	default='12',			help='Labels font size. [Default: %default]')
		
	(options, args) = parser.parse_args()
	
	# Filtering vcf file
	if options.tree == None:
		sys.exit('Please provide a tree file to --tree argument')
	if options.fasta == None:
		sys.exit('Please provide a fasta file to --fasta argument')
	if options.color == None:
		sys.exit('Please provide a color file to --color argument')
	
	# recoding color information
	DicoColor = {}
	file = open(options.color)
	for line in file:
		data = line.split()
		if data:
			if len(data) >= 2:
				DicoColor[data[0]] = data[1]
	file.close()
	
	# recoding accessions in the tree based on the fasta file
	DicoAcc = set()
	file = open(options.fasta)
	for line in file:
		data = line.split()
		if data:
			if data[0][0] == '>':
				DicoAcc.add(data[0][1:])
	file.close()
	
	# Printing New tree
	outfile = open(options.tree+'-'+options.layout+'.tree', 'w')
	outfile.write("#NEXUS\nbegin taxa;\n\tdimensions ntax="+str(len(DicoAcc))+";\n\ttaxlabels\n")
	for acc in DicoAcc:
		if acc in DicoColor:
			outfile.write('\t'+acc+'[&!color='+DicoColor[acc]+']\n')
		else:
			outfile.write('\t'+acc+'\n')
	outfile.write(";\nend;\n\nbegin trees;\n\ttree tree_1 = ")
	
	file = open(options.tree)
	lineToprint = ""
	for line in file:
		data = line.split()
		if data:
			lineToprint += line
	file.close()
	outfile.write(lineToprint+"\nend;\n")
	
	outfile.write('begin figtree;\n')
	outfile.write('\tset appearance.backgroundColorAttribute="Default";\n')
	outfile.write('\tset appearance.backgroundColour=#ffffff;\n')
	outfile.write('\tset appearance.branchColorAttribute="User selection";\n')
	outfile.write('\tset appearance.branchColorGradient=false;\n')
	outfile.write('\tset appearance.branchLineWidth=1.0;\n')
	outfile.write('\tset appearance.branchMinLineWidth=0.0;\n')
	outfile.write('\tset appearance.branchWidthAttribute="Fixed";\n')
	outfile.write('\tset appearance.foregroundColour=#000000;\n')
	outfile.write('\tset appearance.hilightingGradient=false;\n')
	outfile.write('\tset appearance.selectionColour=#2d3680;\n')
	outfile.write('\tset branchLabels.colorAttribute="User selection";\n')
	outfile.write('\tset branchLabels.displayAttribute="Branch times";\n')
	outfile.write('\tset branchLabels.fontName="Agency FB";\n')
	outfile.write('\tset branchLabels.fontSize=8;\n')
	outfile.write('\tset branchLabels.fontStyle=0;\n')
	outfile.write('\tset branchLabels.isShown=false;\n')
	outfile.write('\tset branchLabels.significantDigits=4;\n')
	outfile.write('\tset layout.expansion=0;\n')
	outfile.write('\tset layout.layoutType="'+options.layout+'";\n')
	outfile.write('\tset layout.zoom=447;\n')
	outfile.write('\tset legend.attribute="label";\n')
	outfile.write('\tset legend.fontSize=10.0;\n')
	outfile.write('\tset legend.isShown=false;\n')
	outfile.write('\tset legend.significantDigits=4;\n')
	outfile.write('\tset nodeBars.barWidth=4.0;\n')
	outfile.write('\tset nodeBars.displayAttribute=null;\n')
	outfile.write('\tset nodeBars.isShown=false;\n')
	outfile.write('\tset nodeLabels.colorAttribute="User selection";\n')
	outfile.write('\tset nodeLabels.displayAttribute="Node ages";\n')
	outfile.write('\tset nodeLabels.fontName="Agency FB";\n')
	outfile.write('\tset nodeLabels.fontSize=8;\n')
	outfile.write('\tset nodeLabels.fontStyle=0;\n')
	outfile.write('\tset nodeLabels.isShown=false;\n')
	outfile.write('\tset nodeLabels.significantDigits=4;\n')
	outfile.write('\tset nodeShapeExternal.colourAttribute="User selection";\n')
	outfile.write('\tset nodeShapeExternal.isShown=false;\n')
	outfile.write('\tset nodeShapeExternal.minSize=10.0;\n')
	outfile.write('\tset nodeShapeExternal.scaleType=Width;\n')
	outfile.write('\tset nodeShapeExternal.shapeType=Circle;\n')
	outfile.write('\tset nodeShapeExternal.size=4.0;\n')
	outfile.write('\tset nodeShapeExternal.sizeAttribute="Fixed";\n')
	outfile.write('\tset nodeShapeInternal.colourAttribute="User selection";\n')
	outfile.write('\tset nodeShapeInternal.isShown=false;\n')
	outfile.write('\tset nodeShapeInternal.minSize=10.0;\n')
	outfile.write('\tset nodeShapeInternal.scaleType=Width;\n')
	outfile.write('\tset nodeShapeInternal.shapeType=Circle;\n')
	outfile.write('\tset nodeShapeInternal.size=4.0;\n')
	outfile.write('\tset nodeShapeInternal.sizeAttribute="Fixed";\n')
	outfile.write('\tset polarLayout.alignTipLabels=false;\n')
	outfile.write('\tset polarLayout.angularRange=0;\n')
	outfile.write('\tset polarLayout.rootAngle=0;\n')
	outfile.write('\tset polarLayout.rootLength=100;\n')
	outfile.write('\tset polarLayout.showRoot=true;\n')
	outfile.write('\tset radialLayout.spread=0.0;\n')
	outfile.write('\tset rectilinearLayout.alignTipLabels=false;\n')
	outfile.write('\tset rectilinearLayout.curvature=0;\n')
	outfile.write('\tset rectilinearLayout.rootLength=100;\n')
	outfile.write('\tset scale.offsetAge=0.0;\n')
	outfile.write('\tset scale.rootAge=1.0;\n')
	outfile.write('\tset scale.scaleFactor=1.0;\n')
	outfile.write('\tset scale.scaleRoot=false;\n')
	outfile.write('\tset scaleAxis.automaticScale=true;\n')
	outfile.write('\tset scaleAxis.fontSize=8.0;\n')
	outfile.write('\tset scaleAxis.isShown=false;\n')
	outfile.write('\tset scaleAxis.lineWidth=1.0;\n')
	outfile.write('\tset scaleAxis.majorTicks=1.0;\n')
	outfile.write('\tset scaleAxis.minorTicks=0.5;\n')
	outfile.write('\tset scaleAxis.origin=0.0;\n')
	outfile.write('\tset scaleAxis.reverseAxis=false;\n')
	outfile.write('\tset scaleAxis.showGrid=true;\n')
	outfile.write('\tset scaleBar.automaticScale=true;\n')
	outfile.write('\tset scaleBar.fontSize=10.0;\n')
	outfile.write('\tset scaleBar.isShown=true;\n')
	outfile.write('\tset scaleBar.lineWidth=1.0;\n')
	outfile.write('\tset scaleBar.scaleRange=0.0;\n')
	outfile.write('\tset tipLabels.colorAttribute="User selection";\n')
	outfile.write('\tset tipLabels.displayAttribute="Names";\n')
	outfile.write('\tset tipLabels.fontName="Calibri";\n')
	outfile.write('\tset tipLabels.fontSize='+options.fsize+';\n')
	outfile.write('\tset tipLabels.fontStyle=0;\n')
	outfile.write('\tset tipLabels.isShown=true;\n')
	outfile.write('\tset tipLabels.significantDigits=4;\n')
	outfile.write('\tset trees.order=false;\n')
	outfile.write('\tset trees.orderType="increasing";\n')
	outfile.write('\tset trees.rooting=false;\n')
	outfile.write('\tset trees.rootingType="User Selection";\n')
	outfile.write('\tset trees.transform=false;\n')
	outfile.write('\tset trees.transformType="cladogram";\n')
	outfile.write('end;\n')
	
	outfile.close()
		
if __name__ == "__main__": __main__()
