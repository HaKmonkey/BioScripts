#!/usr/bin/env python3

import csv, sys, re

#############
# Functions #
#############

# Function to reverse a string
def reverse(string):
    string = string[::-1]
    return string

#########################
# Variable Declarations #
#########################

# A list of the input files
# The fist input should be the `.sam file
# The second input should be the `.gff` file
inFile = [sys.argv[1], sys.argv[2]]

# This is the `.sam` file information
samData = []
# Might not need some of these either. Will edit after I figure out how to
# - propperly edit the gff information
samNodes = []
samFlags = []
samStarts = []
samSeqs = []
samCigars = []

# This is the `.gff` file information
gffData = []
gffNodes = []
gffStarts = []
gffEnds = []

# Variable for filtering the `.gff` information
header = 0

# Variable for grouping the gene features by NODE_ID
groups = []
info_groups = []
info = []

# Variable to contain the modified sequences
mod_seq = []

##################
# Parsing Inputs #
##################

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

##################
# Cleaning Files #
##################

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
head = gffData[0:header] # just in case I need this for making the gff file
del gffData[0:header]

#################
# Grabbing Data #
#################

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
    seq = temp[9]
    seq = seq[1:len(seq)-2]
    cigar = temp[5]
    cigar = cigar[1:len(cigar)-2]
    samNodes.append(node)
    samFlags.append(flag)
    samStarts.append(start)
    samSeqs.append(seq)
    samCigars.append(cigar)

# Need [index 0] Node name, [index 3] start of the feature, and [index 9] end
# - of the feature
for i in range(len(gffData)):
    temp = str(gffData[i])
    temp = temp.split(" ")
    node = temp[0]
    node = node[2:len(node)-2]
    start = temp[3]
    start = start[1:len(start)-2]
    end = temp[4]
    end = end[1:len(end)-2]
    gffNodes.append(node)
    gffStarts.append(start)
    gffEnds.append(end)

##################################
# Modifying Sequences & Features #
##################################

# Sorting feature indices by which node they match to
for i in range(len(samNodes)):
    temp = []
    print(samNodes[i])
    for j in range(len(gffNodes)):
        if gffNodes[j] == samNodes[i]:
            temp.append(j)
    groups.append(temp)

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
    if samFlags[i] == '16':
        mod_seq.append(''.join(reverse(temp_seq)))
    else:
        mod_seq.append(''.join(temp_seq))

# Adjusting annotations
for i in range(len(groups)):
    # need to parse through samCigar[i]
    smoke = re.findall(r'([\d]+)([A-Z])', samCigars[i])
    base_num = 1
    indices = groups[i]
    for j in range(len(groups[i])):
        index = int(indices[j])
        start =  int(gffStarts[index])
        end = int(gffEnds[index])
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
                    # checking where the edit falls in the feature
                    if start < base_num <= end:
                        end += 1
                    elif base_num <= start:
                        start += 1
                        end += 1
                    base_num += 1
                    count += 1
            elif op[1] == 'S' or op[1] == 'I' or op[1] == 'H':
                while count != int(op[0]):
                    if base_num <= start:
                        start -= 1
                        end -= 1
                    elif start < base_num <= end:
                        end -= 1
                    count += 1
            else:
                print('P')
            temp = []
            temp.append(index)
            temp.append(start)
            temp.append(end)
        info.append(temp)
