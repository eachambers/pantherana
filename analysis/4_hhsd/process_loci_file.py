# Import libraries
import pandas as pd
import itertools
import csv
import sys
import os
import glob

# Set up args
outfileformat = sys.argv[1] # path and format of output files e.g., "loci/forreri_0.75miss_ldp_locus#.txt"
locuslist = sys.argv[2] # one line per locus number
inputfile = sys.argv[3]

splittingtxt = '//'

# Open and read in relevant files
locikeep = open(locuslist, "r") # open list of loci
locuslines = locikeep.readlines() # read each line of loci to keep so you can index later
file = open(inputfile) # open .loci file from iPyrad
lines = file.read().split(splittingtxt) # split .loci file based on splitting text

# Define function
# The retrieveloc arg will be the indexed file containing all retained loci. The lines arg specifies the pre-split .loci file you read in already.
def output(retrieveloc, lines):
    filename = outfileformat.replace('#', str(retrieveloc.rstrip('\n')))
    fout = open(filename, "w")
    fout.write(lines[int(retrieveloc)])
    fout.close
    
# Run function
# For each indexed kept locus line, retrieve the split locus from the .loci file
for i in range(0, len(locuslines)):
    output(locuslines[i], lines)
