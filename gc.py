#!/usr/bin/env python3

import re, csv, sys

inFile = sys.argv[1]
nodes = []
contigs = []
sequence = []
temp = ""

with open(inFile) as f:
    for l, line in enumerate(f):
        if re.search(r"\>", line):
            nodes.append(line.strip())
            contigs.append(None)
        else:
            contigs.append(line.strip())

contigs.append(None)

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

for i in range(len(sequence)):
    g[i] = sequence[i].count('G')
    c[i] = sequence[i].count('C')
    t[i] = sequence[i].count('T')
    a[i] = sequence[i].count('A')
    gc[i] = ((g[i] + c[i])/len(sequence[i])) * 100

with open('gc_content.csv', 'w') as csvfile:
    fieldnames = ['Node ID', 'GC %', 'Length']
    writer = csv.DictWriter(csvfile, delimiter = ' ', fieldnames = fieldnames)
    writer.writeheader()
    for i in range(len(nodes)):
        writer.writerow({'Node ID': str(nodes[i]) ,'GC %': str(gc[i]), 'Length': len(sequence[i])})
