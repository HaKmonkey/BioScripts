#!/usr/bin/env python3

# will be using positional information from sam file to inform positional
# - changes in the gff data. This data will then be output as a new gff file,
# - that will hopefull work... fingers crossed
## ONLY MODIFY GFF DATA

import csv, sys, re

inFile = []
samData = []
gffData = []
header = 0

# Might not need some of these either. Will edit after I figure out how to
# - propperly edit the gff information
samNodes = []
samStart = []
samSeq = []
samLength = []

# Both of the input files
# Index 0 is the `.sam` file from bwa
# Index 1 is the `.gff` file from prokka
[inFile.append(sys.argv[i]) for i in range(1,3)]

# Storing the `.sam` file data into a list
with open(inFile[0], 'r') as sam:
    samreader = csv.reader(sam, delimiter = '\t')
    # Skip the first 2 lines in the sam file
    [next(samreader) for _ in range(2)]
    for row in samreader:
        samData.append(row)

# Storing the `.gff` file data into a list
with open(inFile[1], 'r') as gff:
    gffreader = csv.reader(gff, delimiter = '\t')
    for row in gffreader:
        gffData.append(row)

# Finding where the sequence information starts and remove it because we don't
# - need it for the changes we want to make
for i in range(len(gffData)):
    if gffData[i] == ['##FASTA']:
        fasta = i
        break

del gffData[fasta:len(gffData)]

# Finding and removing the header to the `.gff` file because we don't need it
for i in range(len(gffData)):
    if re.search(r"\#", str(gffData[i])):
        header += 1

del gffData[0:header]

# Need position 1 [index 0] (Node name), position 4 [index 3] (left most base),
# - and position 10 [index 9] (sequence, for length and whatnot)
for i in range(len(samData)):
    temp = str(samData[i])
    temp = temp.split(" ")
    node = temp[0]
    start = temp[3]
    start = start[1:len(start)-2]
    seq = temp[9]
    seq = seq[1:len(seq)-2]
    node = node[2:len(node)-2]
    samNodes.append(node)
    samStart.append(start)
    samSeq.append(seq)
    samLength.append(len(seq))

# Need position 4 [index 3] (feature start), position 5 [index 4] (feature end)
for i in range(len(gffData)):
    temp = gffData[i]
    # if temp = samNode
    # - change temp[3]
    # - change temp[4]
    # - gffData[i] = temp

print(samNodes[3])
print(samStart[3])
print(samLength[3])
print(samNodes[4])
print(samStart[4])
print(samLength[5])
