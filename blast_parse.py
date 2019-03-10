#!/usr/bin/env python3

import re
import sys
import csv
from astropy.table import Table, Column

# blastn output file
inFile = sys.argv[1]

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

# START
# Query= NODE_1_length_937914_cov_4.753038 
# END
# Effective search space used: 2035474639056

start = []
end = []

for i in range(len(blast)):
  if re.match(r'\[\'Query=', str(blast[i])):
    start.append(i)
  if re.match(r'\[\'Effective', str(blast[i])):
    end.append(i)
      
print(start)
print(len(start))
print(end)
print(len(end))

# Table
## Query
## Sequence(s) that hit
### If more than one sequence, list which is better (higher bit score)
## Number of hits for each sequence
## Average % match and gaps for each sequence
## Average Length of hits


# | Query Sequence | Database Sequence | Hits | Avg Hit Length | Avg Match % | Avg Gap % | Best |
# | -------------- | ----------------- | ---- | -------------- | ----------- | --------- | ---- |
# | NODE_1         | Pseudomonas       | 135  | 40000          | 95%         | 0.5%      | X    |
# | NODE_1         | Annandia          | 5    | 60             | 90%         | 1%        |      |



file = open("parsed_blast.txt","w")

# writing everything but the footer
for i in range(len(blast)-6):
  temp = str(blast[i])
  blast[i] = temp[2:len(temp)-2]
  file.write(str(blast[i]))
  file.write("\n")

file.close()
