#!/usr/bin/env python3

import re, csv, sys

# Usually a .fasta or .fastq file
# - Scaffolds.fasta
inFile = sys.argv[1]
# Nodes that you want to grab from your file
# - 1,2,3
putative_nodes = sys.argv[2]
# The name you want for your output file
# - putative.fasta
outFile = sys.argv[3]

nodes = []
contigs = []
sequence = []
temp = ""

# This loop is creating the nodes and contigs array
# After each node if found a 'None' value is entered in the contigs to separate
# - the sequences
with open(inFile) as f:
  for l, line in enumerate(f):
    if re.search(r"\>", line):
      nodes.append(line.strip())
      contigs.append(None)
    else:
      contigs.append(line.strip())

# This adds a 'None' value at the end so that the sequences can be parsed
# - correctly
contigs.append(None)

# This loop merges the correct contigs into sequences
for i in range(len(contigs)):
  if contigs[i] != None:
    temp = temp + contigs[i]
  elif contigs[i] == None:
    sequence.append(temp)
    temp = ""
sequence.pop(0)

# Cleaning up temporary lists
del contigs
del temp

## Compair second input to list of available node

keep = []

for i in range(len(nodes)):
  temp = str(re.findall(r"\>NODE_(\d+)_", str(nodes[i])))
  temp = temp[2:len(temp)-2]
  if temp in putative_nodes:
    keep.append(i)

del temp

file = open(outFile, "w")

for i in range(len(keep)):
  file.write(nodes[keep[i]])
  file.write("\n")
  file.write(sequence[keep[i]])
  file.write("\n")

file.close()