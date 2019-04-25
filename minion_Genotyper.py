import sys
from Bio import SeqIO

readFile = sys.argv[1]
minimumReadLength = int(sys.argv[2])

print "Filtering sequences longer than",minimumReadLength,"nucleotides"
numSeq= 0
sequences = {}
for seq_record in SeqIO.parse(readFile,"fastq"):
    numSeq+=1
    if numSeq%100000 == 0:
        print numSeq,"reads have been tested...."
    if len(str(seq_record.seq)) > minimumReadLength:
        if not str(seq_record.id) in sequences:
            sequences[str(seq_record.id)] = str(seq_record.seq)





