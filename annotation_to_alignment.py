#!/usr/bin/env python3

import csv, sys, re

# The fist input should be the `.sam file
# The second input should be the `.gff` file
inFile = [sys.argv[1], sys.argv[2]]
samData = []
gffData = []
supAlign = []
header = 0

# Might not need some of these either. Will edit after I figure out how to
# - propperly edit the gff information
samNodes = []
samFlags = []
samStart = []
samSeq = []
samCigar = []

# Function to reverse a string
def reverse(string):
    string = string[::-1]
    return string

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
    flag = temp[1]
    flag = flag[1:len(flag)-2]
    start = temp[3]
    start = start[1:len(start)-2]
    seq = temp[9]
    seq = seq[1:len(seq)-2]
    node = node[2:len(node)-2]
    cigar = temp[5]
    cigar = cigar[1:len(cigar)-2]
    samNodes.append(node)
    samFlags.append(flag)
    samStart.append(start)
    samSeq.append(seq)
    samCigar.append(cigar)

# Filter `.sam` data
# `.sam` enteries that have the 2048 flag (0x800) are supplimentary alignments
# - and are clipped when viewing the alignments
# So these entries need to be removed so the annotations can be ligned up with
# - the parts of the nodes that are present in the actual alignment
for i in range(len(samNodes)):
    if samFlags[i] == '2048':
        supAlign.append(i)

# will probably just modify lists in place once everything is working
# this way there won't be so many lists taking up memory :)
node = []
start = []
seq = []
flag = []
cigar = []

for i in range(len(samNodes)):
    if  (i in supAlign) == False:
        node.append(samNodes[i])
        start.append(samStart[i])
        seq.append(samSeq[i])
        flag.append(samFlags[i])
        cigar.append(samCigar[i])

# Now we need to adjust the sequences that remain to match their cigar string so
# - that everything will match nicely with the alignment
# Nodes that have the 16 flag (0x10) are mapped in reverse and must therefore
# - be flipped before they are modified, again, so they match the alignment

# Case1 ( M/X/=):
# - start at the specified mapping position, set counter to 1
# - Add 1 to both the counts of the bases from that position and the counter.
# - Move to the next position.
# - Repeat this process till counter is the same as the number associated with
# - the operator.
## SIMPLY COUNT THE BASES PRESENT

# Case2 (N/D):
# - Move the specified mapping position by the number associated with
# - the operator.
## ADD SPACE INTO READ

# Case3 (I/S/H/P):
# - Do nothing
## P ADDS SPACE TO READ ? (DON'T KNOW YET, NO READS HAVE HAD PADDING!)
## S/H CLIP FROM READ (SHORTEN IT)
## I REMOVES BASES FROM READ (looks that way in tablet and artemis)

mod_seq = []
mod_seq_test = []

for i in range(len(seq)):
    smoke = re.findall(r'([\d]+)([A-Z])', cigar[i])
    base_num = 1
    if flag[i] == '16':
        temp_seq = list(reverse(seq[i]))
        print('-')
    else:
        temp_seq = list(seq[i])
        print('+')
    print(node[i])
    print("start len seq: " + str(len(temp_seq)))
    print("start base num: " + str(base_num))
    for j in range(len(smoke)):
        count = 0
        op = smoke[j]
        if op[1] == 'M' or op[1] == 'X' or op[1] == '=':
            # iterate through the sequence without modifying the base
            while count != int(op[0]):
                base_num += 1
                count += 1
            #print("iter")
            #print("current base num: " + str(base_num))
        elif op[1] == 'D':
            # add an N to the read to inflate its length
            while count != int(op[0]):
                temp_seq.insert(base_num-1, '*')
                base_num += 1
                count += 1
            #print("add deletion")
            #print("current base num: " + str(base_num))
        elif op[1] == 'N':
            while count != int(op[0]):
                temp_seq.insert(base_num-1, '-')
                base_num += 1
                count += 1
        elif op[1] == 'S' or op[1] == 'H' or op[1] == 'I':
            # Remove the bases from a range
            while count != int(op[0]):
                del temp_seq[base_num-1]
                count +=1
            #print("remove")
            #print("current base num: " + str(base_num))
        else:
            # I don't yet know how padding affects the reads for this...
            print("P")
    mod_seq.append(str(temp_seq))
    print("end len seq: " + str(len(temp_seq)))
    print("end base num: " + str(base_num))
    print(" ")

print(cigar[23])

##### TEST LOOP WITH SEQUENCE INFLATION

for i in range(len(seq)):
    smoke = re.findall(r'([\d]+)([A-Z])', cigar[i])
    base_num = 1
    if flag[i] == '16':
        temp_seq = list(reverse(seq[i]))
        print('-')
    else:
        temp_seq = list(seq[i])
        print('+')
    print(node[i] + " INFLATION TEST")
    print("start len seq: " + str(len(temp_seq)))
    print("start base num: " + str(base_num))
    for j in range(len(smoke)):
        count = 0
        op = smoke[j]
        if op[1] == 'M' or op[1] == 'X' or op[1] == '=':
            # iterate through the sequence without modifying the base
            while count != int(op[0]):
                base_num += 1
                count += 1
            #print("iter")
            #print("current base num: " + str(base_num))
        elif op[1] == 'D':
            # add an N to the read to inflate its length
            while count != int(op[0]):
                temp_seq.insert(base_num-1, '*')
                base_num += 1
                count += 1
            #print("add deletion")
            #print("current base num: " + str(base_num))
        elif op[1] == 'N':
            while count != int(op[0]):
                temp_seq.insert(base_num-1, '-')
                base_num += 1
                count += 1
        elif op[1] == 'S' or op[1] == 'H' or op[1] == 'I':
            # Remove the bases from a range
            while count != int(op[0]):
                temp_seq.insert(base_num-1, '$')
                base_num += 1
                count += 1
            #print("remove")
            #print("current base num: " + str(base_num))
        else:
            # I don't yet know how padding affects the reads for this...
            print("P")
    mod_seq_test.append(str(temp_seq))
    print("end len seq: " + str(len(temp_seq)))
    print("end base num: " + str(base_num))
    print(" ")

print(cigar[23])
print(len(mod_seq[23]))
print(mod_seq[23])
print(" ")
print(cigar[23])
print(len(mod_seq_test[23]))
print(mod_seq_test[23])

# Need position 4 [index 3] (feature start), position 5 [index 4] (feature end)
for i in range(len(gffData)):
    temp = gffData[i]
    # if temp = samNode
    # - change temp[3]
    # - change temp[4]
    # - gffData[i] = temp

#for i in range(len(node)):
    #print(node[i])
    #print(start[i])
    #print(cigar[i])
    #print(flag[i])
    #print(" ")
