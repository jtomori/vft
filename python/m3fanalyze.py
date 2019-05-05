#!/usr/bin/env python

"""
-----------------------------------------------------------------------------
This source file has been developed within the scope of the
Technical Director course at Filmakademie Baden-Wuerttemberg.
http://technicaldirector.de
 
Written by Juraj Tomori.
Copyright (c) 2019 Animationsinstitut of Filmakademie Baden-Wuerttemberg
-----------------------------------------------------------------------------
"""

import os
import fnmatch

def countInstructions(filePath):
	with open(filePath) as file:
		lines = file.readlines()
		lines = [n.replace(" ","") for n in lines] # kill spaces for finding exact line

		# find start and ending lines of code section
		start = lines.index("[CODE]\n")
		end = lines.index("[END]\n")
		length = end-start-1

		# create a list with lines containing instructions
		code = []

		for i in range(length):
			code.append( lines[start + i + 1].replace("\n", "") )

		# concat the list and divide by 2 (hexadecimal instrcutions)
		return len( "".join(code) )/2

# path to the formulas folder
path = "/home/jtomori/Downloads/mb3d-master/M3Formulas"

# find files only
files = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]

# keep only .m3f files, but not JIT
files = fnmatch.filter(files, "*.m3f")
files = [n for n in files if not fnmatch.fnmatch(n.lower(), "*jit*")]

# a list containing amount of instructions per file
instructions = []

for file in files:
	instructions.append( countInstructions( os.path.join(path, file) ) )

# generate an index list
sortedIndices = sorted(range(len(instructions)), key=lambda k: instructions[k])

# print from longest to shortest code formulas
for index in sortedIndices[::-1]:
	print files[index] + " (" + str(instructions[index]) + " instructions)"