#!/usr/bin/env python3

##################
# Import & Input #
##################

import re
import sys
import csv
import numpy as np
from astropy.table import Table, Column

# blastn output file
inFile = sys.argv[1]

#############
# Functions #
#############

# Function that splits the match percentages into different groups
def splitz(seq, sep):
  group = []
  for num in seq:
    if num != sep:
      group.append(num)
    elif group:
      yield group
      group = []

#################
# Parsing Blast #
#################

blast = []

with open(inFile, 'r') as b:
    blastreader = csv.reader(b, delimiter = '\t')
    # Removing the header
    [next(blastreader) for k in range(11)]
    for row in blastreader:
      # Removing the sequences and 'Score' lines
      # - (No length of match %)
      # re.match(r'\[\' Score =', str(row)) or 
      if re.match(r'\[\'Query\s+[0-9]+', str(row)) or re.match(r'\[\'Sbjct\s+[0-9]+', str(row)) or re.match(r'\[\'\s+\|+|\s+', str(row)) or re.match(r'\[\]', str(row)) or re.match(r'\[\' Strand=', str(row)) or re.match(r'\[\'\s+Score', str(row)) or re.match(r'\[\'Length=' ,str(row)) or re.match(r'\[\'Lambda' ,str(row)) or re.match(r'\[\'Gapped' ,str(row)) or re.match(r'\[\'\s+[0-9]+' ,str(row)):
        continue
      elif re.match(r'\[\'Effective search space used\: [0-9]+\'\]', str(row)) and re.match(r'\[\'\*\*\*\*\* No hits found \*\*\*\*\*\'\]', str(blast[-1])):
        del blast[-2:]
      # Removing the entries with no hits
      else:  
        blast.append(row)

b.close()

########################
# Finding Query Chunks #
########################

start_chunk = []
end_chunk = []

for i in range(len(blast)):
  if re.match(r'\[\'Query=', str(blast[i])):
    start_chunk.append(i)
  if re.match(r'\[\'Effective', str(blast[i])):
    end_chunk.append(i)

########################
# Parsing Query Chunks #
########################

file = open("parsed_blast.txt","w")

for i in range(len(start_chunk)):

  #########################################################
  # Getting the Query and Significant Alignment Sequences #
  #########################################################

  query = ""
  start_sig = None
  end_sig = None

  for j in range(start_chunk[i], end_chunk[i]+1):
    if re.match(r'\[\'Query= ', str(blast[j])):
      query = str(blast[j])
      query = query[2:len(query)-2]
    if re.match(r'\[\'Sequences producing significant alignments', str(blast[j])):
      start_sig = j
    if re.match(r'\[\'\> ', str(blast[j])):
      end_sig = j
      break
  
  #############################
  # Finding Hits Per Sequence #
  #############################
  
  seq = [] # For sequences column
  hits = []
  x = False

  for j in range(start_chunk[i], end_chunk[i]+1):
    if re.match(r'\[\'\> ', str(blast[j])):
      temp_s = str(blast[j])
      temp_s = temp_s[4:len(temp_s)-2]
      seq.append(temp_s)
      x = True
      continue
    while x:
      if re.match(r'\[\'\> ', str(blast[j])) or j == end_chunk[i]:
        x = False
        hits.append(" ")
        j += 1
      else:
        hits.append(str(blast[j]))
        j += 1
  
  #########################
  # Pulling Hits and Gaps #
  #########################
  
  match = [] 
  gaps = []
  length = []

  for j in range(len(hits)):
    temp_m = str(re.findall(r'Identities = [0-9]+/[0-9]+ \(([0-9]+)\%\)', hits[j]))
    temp_m = temp_m[2:len(temp_m)-2]
    match.append(temp_m)
    temp_g = str(re.findall(r'Gaps = [0-9]+/[0-9]+ \(([0-9]+)\%\)', hits[j]))
    temp_g = temp_g[2:len(temp_g)-2]
    gaps.append(temp_g)
    temp_l = str(re.findall(r'Identities = [0-9]+/([0-9]+) \([0-9]+\%\)', hits[j]))
    temp_l = temp_l[2:len(temp_l)-2]
    length.append(temp_l) 

  ############
  # Grouping #
  ############

  # match and gap % grouped by which sequence they hit
  matchp = list(splitz(match, '')) # for avg match % column
  gapp = list(splitz(gaps, '')) # for avg gap % column
  length = list(splitz(length, '')) # for avg hit length column
  num_hits = [] # for the hits column

  # gets the mean for each group
  for j in range(len(matchp)):
    matchp[j] = [int(x) for x in matchp[j]]
    num_hits.append(len(matchp[j]))
    matchp[j] = np.around(np.mean(matchp[j]), decimals = 3)
    gapp[j] = [int(x) for x in gapp[j]]
    gapp[j] = np.around(np.mean(gapp[j]), decimals = 3)
    length[j] = [int(x) for x in length[j]]
    length[j] = np.around(np.mean(length[j]), decimals = 3)

  ####################
  # Printing a Table #
  ####################

  header1 = "Sequence" + " " * (len(max(seq, key=len))-8)
  header2 = "Hits" + " " * (len(max(str(num_hits), key=len))-5)
  header3 = "Avg. Hit Length" + " " * (len(max(str(length), key=len))-15)
  header4 = "Avg. Match %" + " " * (len(max(str(matchp), key=len))-12)
  header5 = "Avg. Gaps %" + " " * (len(max(str(gapp), key=len))-11)

  sep1 = '-' * len(header1)
  sep2 = '-' * len(header2)
  sep3 = '-' * len(header3)
  sep4 = '-' * len(header4)
  sep5 = '-' * len(header5)

  file.write(query)
  file.write("\n\n")

  file.write("|\t")
  file.write(header1)
  file.write("\t|\t")
  file.write(header2)
  file.write("\t|\t")
  file.write(header3)
  file.write("\t|\t")
  file.write(header4)
  file.write("\t|\t")
  file.write(header5)
  file.write("\t|")
  file.write("\n")

  file.write("|\t")
  file.write(sep1)
  file.write("\t|\t")
  file.write(sep2)
  file.write("\t|\t")
  file.write(sep3)
  file.write("\t|\t")
  file.write(sep4)
  file.write("\t|\t")
  file.write(sep5)
  file.write("\t|")
  file.write("\n")

  for j in range(len(seq)):  
    file.write("|\t")
    file.write(str(seq[j]))
    file.write(" " * (len(header1) - len(str(seq[j]))))
    file.write("\t|\t")
    file.write(str(num_hits[j]))
    file.write(" " * (len(header2) - len(str(num_hits[j]))))
    file.write("\t|\t")
    file.write(str(length[j]))
    file.write(" " * (len(header3) - len(str(length[j]))))
    file.write("\t|\t")
    file.write(str(matchp[j]))
    file.write(" " * (len(header4) - len(str(matchp[j]))))
    file.write("\t|\t")
    file.write(str(gapp[j]))
    file.write(" " * (len(header5) - len(str(gapp[j]))))
    file.write("\t|")
    file.write("\n")

  file.write("\n\n\n")    
 
file.close()