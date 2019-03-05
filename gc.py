#!/usr/bin/env python3

import re, csv, sys

inFile = sys.argv[1]
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
for i in range(0, len(contigs)):
    if contigs[i] != None:
        temp = temp + contigs[i]
    elif contigs[i] == None:
        sequence.append(temp)
        temp = ""
sequence.pop(0)

g = [None] * len(sequence)
c = [None] * len(sequence)
t = [None] * len(sequence)
a = [None] * len(sequence)
gc = [None] * len(sequence)

# This loop calculates the number of different bases for each sequence and then
# - calculates the gc percent
for i in range(len(sequence)):
    g[i] = sequence[i].count('G')
    c[i] = sequence[i].count('C')
    t[i] = sequence[i].count('T')
    a[i] = sequence[i].count('A')
    gc[i] = ((g[i] + c[i])/len(sequence[i])) * 100

# This block creates the csv file. You can change the delimeter to a ',' if you
# - would like, but it only accepts something that is one character in length
# You can then take the csv file and import it into Excel or Google Sheets
# Google Sheets will automatically recognize the delimeters, but you may have to
# - specify in Excel
with open('gc_content.csv', 'w') as csvfile:
    fieldnames = ['Node ID', 'GC %', 'Length']
    writer = csv.DictWriter(csvfile, delimiter = ' ', fieldnames = fieldnames)
    writer.writeheader()
    for i in range(len(nodes)):
        writer.writerow({'Node ID': str(nodes[i]) ,'GC %': str(gc[i]),
        'Length': len(sequence[i])})
