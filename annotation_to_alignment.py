#!/usr/bin/env python3

##################
# Import & Input #
##################

import csv, sys, re

# A list of the input files
# The fist input should be the `.sam file
# The second input should be the `.gff` file
inFile = [sys.argv[1], sys.argv[2]]

# Warnings
if inFile[0].endswith('.sam') == False:
  print("The first input needs to be a .sam file \n")
  print("EX: python3 annotation_to_alignment.py alignment.sam prokka.gff")
  sys.exit()
elif inFile[1].endswith('.gff') == False:
  print("The second input needs to be a .gff file \n")
  print("EX: python3 annotation_to_alignment.py alignment.sam prokka.gff")
  sys.exit()
elif len(inFile) != 2:
  print("You must input a .sam file followed by a .gff file \n")
  print("EX: python3 annotation_to_alignment.py alignment.sam prokka.gff")
  sys.exit()

#############
# Functions #
#############

# Function to reverse a string
def reverse(string):
  string = string[::-1]
  return string

####################
# Parsing Inputs 1 #
####################

# This is the `.sam` file
samData = []

with open(inFile[0], 'r') as sam:
  samreader = csv.reader(sam, delimiter = '\t')
  # Skip the first 2 lines in the sam file
  [next(samreader) for i in range(2)]
  for row in samreader:
    samData.append(row)
sam.close()

#####################
# Grabbing SAM Data #
#####################

# This is the `.sam` file information
samNodes = []
samFlags = []
samStarts = []
samSeqs = []
samCigars = []

# Need [index 0] Node name, [index 3] left most base,[index 5] the cigar string
# - , and [index 9] sequence
for i in range(len(samData)):
  temp = str(samData[i])
  temp = temp.split(" ")
  node = temp[0]
  node = node[2:len(node)-2]
  flag = temp[1]
  flag = flag[1:len(flag)-2]
  start = temp[3]
  start = start[1:len(start)-2]
  start = int(start)
  seq = temp[9]
  seq = seq[1:len(seq)-2]
  cigar = temp[5]
  cigar = cigar[1:len(cigar)-2]
  samNodes.append(node)
  samFlags.append(flag)
  samStarts.append(start)
  samSeqs.append(seq)
  samCigars.append(cigar)

# Cleaning up temporary lists
del temp
del node
del flag
del start
del seq
del cigar

#######################
# Modifying Sequences #
#######################

# Variable to contain the modified sequences
mod_seq = []

# Variable to contain end position in the alignment
samEnds = []

# Now we need to adjust the sequences that remain to match their cigar string so
# - that everything will match nicely with the alignment
# Nodes that have the 16 flag (0x10) are mapped in reverse and must therefore
# - be flipped before they are modified, again, so they match the alignment
for i in range(len(samSeqs)):
  smoke = re.findall(r'([\d]+)([A-Z])', samCigars[i])
  base_num = 1
  if samFlags[i] == '16':
    temp_seq = list(reverse(samSeqs[i]))
  else:
    temp_seq = list(samSeqs[i])
  for j in range(len(smoke)):
    count = 0
    op = smoke[j]
    if op[1] == 'M' or op[1] == 'X' or op[1] == '=':
      # iterate through the sequence without modifying the base
      while count != int(op[0]):
        base_num += 1
        count += 1
    elif op[1] == 'D':
      # add an N to the read to inflate its length
      while count != int(op[0]):
        temp_seq.insert(base_num-1, '*')
        base_num += 1
        count += 1
    elif op[1] == 'N':
      while count != int(op[0]):
        temp_seq.insert(base_num-1, '-')
        base_num += 1
        count += 1
    elif op[1] == 'S' or op[1] == 'I':
      # Remove the bases from a range
      while count != int(op[0]):
        del temp_seq[base_num-1]
        count += 1
    elif op[1] == 'H':
      while count != int(op[0]):
        count += 1
    else:
      # I don't yet know how padding affects the reads for this...
      print('P')
  samEnds.append(samStarts[i] - 1 + len(temp_seq))
  if samFlags[i] == '16':
    mod_seq.append(''.join(reverse(temp_seq)))
  else:
    mod_seq.append(''.join(temp_seq))

# Cleaning up temporary lists
del smoke
del base_num
del temp_seq
del count
del op

####################
# Parsing Inputs 2 #
####################

# This is the `.gff` file
gffData = []

# Storing the `.gff` file data into a list
with open(inFile[1], 'r') as gff:
  gffreader = csv.reader(gff, delimiter = '\t')
  for row in gffreader:
    gffData.append(row)
gff.close()

#####################
# Cleaning GFF File #
#####################

# Finding where the sequence information starts and remove it because we don't
# - need it for the changes we want to make
for i in range(len(gffData)):
  if gffData[i] == ['##FASTA']:
    fasta = i
    break
if type(fasta) == int:
    del gffData[fasta:len(gffData)]

# Variable for filtering the `.gff` information
header = 0

# Finding and removing the header to the `.gff` file because we don't need it
for i in range(len(gffData)):
  if re.search(r"\#", str(gffData[i])):
    header += 1
del gffData[0:header]

# Cleaning up temporary lists
if type(fasta) == int:
  del fasta
del header

#####################
# Grabbing GFF Data #
#####################

# This is the `.gff` file information
gffNodes = []
gffTypes = []
gffStarts = []
gffEnds = []
gffScores = []
gffStrands = []
gffPhases = []
gffGenes = []

# Need [index 0] Node name, [index 3] start of the feature, and [index 9] end
# - of the feature
for i in range(len(gffData)):
  temp = str(gffData[i])
  temp = temp.split(" ")
  node = temp[0]
  node = node[2:len(node)-2]
  ftype = temp[2]
  ftype = ftype[1:len(ftype)-2]
  start = temp[3]
  start = start[1:len(start)-2]
  end = temp[4]
  end = end[1:len(end)-2]
  score = temp[5]
  score = score[1:len(score)-2]
  strand = temp[6]
  strand = strand[1:len(strand)-2]
  phase = temp[7]
  phase = phase[1:len(phase)-2]
  gene = str(re.findall(r'(Name=[A-Za-z0-9]+);' ,temp[8]))
  gene = gene[2:len(gene)-2]
  gffNodes.append(node)
  gffTypes.append(ftype)
  gffStarts.append(start)
  gffEnds.append(end)
  gffScores.append(score)
  gffStrands.append(strand)
  gffPhases.append(phase)
  gffGenes.append(gene)

# Cleaning up temporary lists
del temp
del node
del ftype
del start
del end
del score
del strand
del phase
del gene

######################
# Modifying Features #
######################

# Variable for grouping the gene features by NODE_ID
groups = []

# Sorting feature indices by which node they match to
for i in range(len(samNodes)):
  temp = []
  for j in range(len(gffNodes)):
    if gffNodes[j] == samNodes[i]:
      temp.append(j)
  groups.append(temp)

# Variable for storing new position info for all features
info = []

# Adjusting annotations
for i in range(len(groups)):
  # need to parse through samCigar[i]
  smoke = re.findall(r'([\d]+)([A-Z])', samCigars[i])
  indices = groups[i]
  for j in range(len(indices)):
    index = int(indices[j])
    start =  int(gffStarts[index])
    end = int(gffEnds[index])
    base_num = 1
    temp = []
    for k in range(len(smoke)):
      count = 0
      op = smoke[k]
      # This won't affect any of the features because it is simply
      # - 'reading' through the sequence
      if op[1] == 'M' or op[1] == 'X' or op[1] == '=':
        while count != int(op[0]):
          base_num += 1
          count += 1
      elif op[1] == 'D' or op[1] == 'N':
        while count != int(op[0]):
          base_num += 1
          count += 1
          # checking where the edit falls in the feature
          if start < base_num <= end:
            end += 1
          elif base_num <= start:
            start += 1
            end += 1
      elif op[1] == 'S' or op[1] == 'I' or op[1] == 'H':
        while count != int(op[0]):
          count += 1
          if base_num <= start:
            start -= 1
            end -= 1
          elif start < base_num <= end:
            end -= 1
      else:
        print('P')
    temp.append(index)
    temp.append(start)
    temp.append(end)
    info.append(temp)

# Cleaning up temporary lists
del temp
del smoke
del base_num
del indices
del index
del start
del end
del count
del op

############################
# Adjusting SAM Node Names #
############################

# Variable for storing original node number
id = []

# Pulling node number for each index
for i in range(len(samNodes)):
  temp = re.findall(r'NODE_(\d+)_', samNodes[i])
  temp = str(temp)
  id.append(temp)

# Variable for finding chimeric alignments
count = 1
bool = False

# Changing the names of nodes to fit the NODE_[\d]+_p[\d]+
# - NODE_1_p1
# - NODE_2_p1
# - NODE_2_p2
for i in range(len(samNodes)):
  if bool == True and id[i] == id[i-1]:
    count += 1
    samNodes[i] = ('%s%s%s%d') %("NODE_", str(id[i])[2:len(id[i])-2], "_p", count)
  else:
    count = 1
    samNodes[i] = ('%s%s%s%d') %("NODE_", str(id[i])[2:len(id[i])-2], "_p", count)
    bool = True

# Cleaning up temporary lists
del temp
del id
del count
del bool

###############################
# Parsing Input 1 as Original #
###############################

# This is the `.sam` file unedited
originalSam = []

# Storing the `.sam` file data into a list
with open(inFile[0], 'r') as sam:
  samreader = csv.reader(sam, delimiter = '\t')
  for row in samreader:
    originalSam.append(row)
        
sam.close()

###############################
# Modify Original SAM Headers #
###############################

samString = ""
temp = originalSam[0]

for i in range(len(temp)):
  samString += temp[i]
  samString += "\t"
originalSam[0] = samString

samString = ""
temp = originalSam[1]

for i in range(len(temp)):
  samString += temp[i]
  samString += "\t"
originalSam[1] = samString

# Cleaning up temporary lists
del samString
del temp

################################
# Modify Original SAM Node IDs #
################################

# Changing node names in `.sam` file to match `.gff` file
for i in range(2,len(originalSam)):
  temp = str(originalSam[i])
  temp = temp.split(' ')
  temp[0] = samNodes[i-2]
  for j in range(1,len(temp)):
    spot = temp [j]
    spot = spot[1:len(spot)-2]
    temp[j] = spot
  temp = '\t'.join(temp)
  originalSam[i] = temp

# Cleaning up temporary lists
del temp
del spot

##########################################
# Adjusting GFF Node Names & Feature IDs #
##########################################

# Variables for storing the 'new' sets of nodes
newNodes = []
used_index = []

for i in range(len(info)):
  temp = info[i]
  index = temp[0]
  name = re.findall(r'NODE_(\d+)_', gffNodes[index])
  name = str(name)
  used_index.append(index)
  count = used_index.count(index)
  newNodes.append(('%s%s%s%d') %("NODE_", name[2:len(name)-2], "_p", count))

# Cleaning up temporary lists
del temp
del index
del name
del count
del used_index

####################################
# Duplicate Necessary Feature Info #
####################################

# Variable for storing 'new' types
newTypes = []
newScores = []
newStrands = []
newPhases = []
newGenes = []
newStarts = []
newEnds = []
newID = []

for i in range(len(info)):
  temp = info[i]
  index = temp[0]
  start = temp[1]
  end = temp[2]
  newTypes.append(gffTypes[index])
  newScores.append(gffScores[index])
  newStrands.append(gffStrands[index])
  newPhases.append(gffPhases[index])
  newGenes.append(gffGenes[index])
  newStarts.append(start)
  newEnds.append(end)
  if gffGenes[index] == '':
    newID.append("ID=MPLGGMNL_" + str(i))
  else:
    newID.append("ID=MPLGGMNL_" + str(i) + ";")

# Cleaning up temporary lists
del temp
del index

######################
# Write New SAM File #
######################

# Making a new `.sam` file with updated names
file = open("a2a.sam","w")

# Writing the lines to the `.sam` file
for i in range(len(originalSam)):
  file.write(originalSam[i])
  file.write("\n")

# Closing the new `.sam` file
file.close()

#######################################
# Storing Index of Features to Remove #
#######################################

# Variable that stores index of features that don't map to alignment, so they
# - can be removed before writing to file
emp = []

for i in range(len(info)):
  temp = info[i]
  start = temp[1]
  end = temp[2]
  if start == 0 or end == 0 or start == end:
    emp.append(i)

# Cleaning up temporary lists
del temp
del start
del end

######################
# Write New GFF File #
######################

# Making a new `.gff` file
file = open("a2a.gff", "w")

# Writing the header to the `.gff` file
file.write("##gff-version 3 \n")
'''
for i in range(len(samNodes)):
  file.write("##sequence-region\t" + samNodes[i] + "\t1\t" + str(len(samSeqs[i])) + "\n")

print(len(info))
'''

## SOME FEATURES ARE WAY OUTSIDE THE BOUNDS OF THE ALIGNMENT
# - Check `a2a.gff`
# - Check Adjusting Features
## newNodes has the correct nodes and duplicates...
## info has the correct indexes and duplicates as well
## seeing 40 NODE_22_p1 instead of 20...

#for i in range(len(newNodes)):
  #print(newNodes[i])

#print(newNodes[127])

#print(len(samSeqs))
#print(len(newNodes))
#print(newNodes[146])
#print(newNodes[166])

# Adjusting starts and ends of gene features to match alignment
for i in range(len(samSeqs)):
  #print(samNodes[i])
  for j in range(len(newNodes)):
    temp = info[j]
    index = temp[0]
    if newNodes[j] == samNodes[i]:
      #print(samNodes[i])
      #print(j)
      ## NODE_22_p1 give 40 instead of 20
      ## NODE_112_p1 gives 33 instead of 25
      ## NODE_312_p1 gives 30 instead of 16
      ## NODE_1108 gives 19 instead of 12
      ## NODE_1297_ gives 81 instead of 11
      ## NODE_2333_ gives 6 instead of 7
      ## NODE_3154_, NODE_3811_, NODE_6448_, NODE_6540_, NODE_13192_, NODE_22013_, NODE_43464_ not showing hits
      #print(samStarts[i])
      newStarts[j] = int(newStarts[j]) + samStarts[i]
      #print(newStarts[index])
      newEnds[j] = int(newEnds[j]) + samStarts[i]

#for i in range(len(newStarts)):
  #print(newStarts[i])

# Writing the features to the `.gff` files
for i in range(len(newNodes)):
  if (i in emp) == False:
    file.write(newNodes[i] + "\ta2a\t" + newTypes[i] + "\t" + str(newStarts[i]) + "\t" + str(newEnds[i]) + "\t" + str(newScores[i]) + "\t" + newStrands[i] + "\t" + str(newPhases[i]) + "\t" + newID[i] + newGenes[i] + "\n")
'''
# Writing the ##FASTA section to the file
file.write("##FASTA \n")

# Wtiring node names and sequences to the file
for i in range(len(samNodes)):
  file.write(">" + samNodes[i] + "\n")
  temp = mod_seq[i]
  for j in range(0,len(temp), 60):
    file.write(temp[j:j+60] + "\n")

# Cleaning up temporary lists
del temp
'''
# Closing the new `.gff` file
file.close()
