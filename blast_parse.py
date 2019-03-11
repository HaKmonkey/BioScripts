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

'''
query = ""
start_sig = None
end_sig = None

for i in range(len(start_chunk)):
  for j in range(start_chunk[i], end_chunk[i] + 1):
    if re.match(r'Query= ', blast[j]):
      query = str(blast[j])
    if re.match(r'\[\'Sequences producing significant alignments', blast[j]):
      start_sig = j
    if re.match(r'\[\'\> ', blast[j]):
      end_sig = j
'''

##############################################################
# Getting the Query and Significant Alignment Sequences TEST #
##############################################################

query = ""
start_sig = None
end_sig = None

for j in range(start_chunk[0], end_chunk[0]+1):
  if re.match(r'\[\'Query= ', str(blast[j])):
      query = str(blast[j])
  if re.match(r'\[\'Sequences producing significant alignments', str(blast[j])):
    start_sig = j
  if re.match(r'\[\'\> ', str(blast[j])):
    end_sig = j
    break

##################################
# Finding Hits Per Sequence TEST #
##################################

seq = [] # For sequences column
hits = []
x = False

for k in range(start_chunk[0], end_chunk[0]+1):
  if re.match(r'\[\'\> ', str(blast[k])):
    seq.append(str(blast[k]))
    x = True
    continue
  while x:
    if re.match(r'\[\'\> ', str(blast[k])) or k == end_chunk[0]:
      x = False
      hits.append(" ")
      k += 1
    else:
      hits.append(str(blast[k]))
      k += 1

##############################
# Pulling Hits and Gaps TEST #
##############################

match = [] 
gaps = []
length = []

for i in range(len(hits)):
  temp_m = str(re.findall(r'Identities = [0-9]+/[0-9]+ \(([0-9]+)\%\)', hits[i]))
  temp_m = temp_m[2:len(temp_m)-2]
  match.append(temp_m)
  temp_g = str(re.findall(r'Gaps = [0-9]+/[0-9]+ \(([0-9]+)\%\)', hits[i]))
  temp_g = temp_g[2:len(temp_g)-2]
  gaps.append(temp_g)
  temp_l = str(re.findall(r'Identities = [0-9]+/([0-9]+) \([0-9]+\%\)', hits[i]))
  temp_l = temp_l[2:len(temp_l)-2]
  length.append(temp_l)


#################
# Grouping TEST #
#################

# match and gap % grouped by which sequence they hit
matchp = list(splitz(match, '')) # for avg match % column
gapp = list(splitz(gaps, '')) # for avg gap % column
length = list(splitz(length, '')) # for avg hit length column
num_hits = [] # for the hits column

# gets the mean for each group
for i in range(len(matchp)):
  matchp[i] = [int(x) for x in matchp[i]]
  num_hits.append(len(matchp[i]))
  matchp[i] = np.mean(matchp[i])
  gapp[i] = [int(x) for x in gapp[i]]
  gapp[i] = np.mean(gapp[i])
  length[i] = [int(x) for x in length[i]]
  length[i] = np.mean(length[i])

#########################
# Printing a Table TEST #
#########################

t = Table()
#t['Sequence'] = seq
t['Hits'] = num_hits
t['Avg. Hit Length'] = length
t['Avg. Match %'] = matchp
t['Avg. Gaps %'] = gapp
#t['Best']

print(seq)
print(num_hits)
print(length)
print(matchp)
print(gapp)
print(" ")

print(t)

# Table
## Query
## Sequence(s) that hit
### If more than one sequence, list which is better (higher bit score)
## Number of hits for each sequence
## Average % match and gaps for each sequence
## Average Length of hits


# NODE_1
#
# | Sequence    | Hits | Avg Hit Length | Avg Match % | Avg Gap % | Best |
# | ----------- | ---- | -------------- | ----------- | --------- | ---- |
# | Pseudomonas | 135  | 40000          | 95%         | 0.5%      | X    |
# | Annandia    | 5    | 60             | 90%         | 1%        |      |


'''
file = open("parsed_blast.txt","w")

# writing everything but the footer
for i in range(len(blast)-6):
  temp = str(blast[i])
  blast[i] = temp[2:len(temp)-2]
  file.write(str(blast[i]))
  file.write("\n")

file.close()
'''
