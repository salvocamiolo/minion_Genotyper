import sys
import os
from Bio import SeqIO

detectionTresholdPercentage = 0.002 #To be diminished to increse sensitivity


orderedHyperLoci = ["rl5a","rl6","rl12","rl13","ul1","ul9","ul11","ul20","ul73","ul74","ul120","ul139","ul146"]

kmerDB = "./kmerDB/mainDB_seqs_filtered.txt"

readFile = sys.argv[1]
minimumReadLength = int(sys.argv[2])
outputPrefix = sys.argv[3]
print "Starting elaboration....\n"
print "Filtering sequences longer than",minimumReadLength,"nucleotides\n"
numSeq= 0
foundSequences = 0
sequences = {}
for seq_record in SeqIO.parse(readFile,"fastq"):
    numSeq+=1
    if numSeq%100000 == 0:
        print numSeq,"reads have been tested...."
    if len(str(seq_record.seq)) > minimumReadLength:
        if not str(seq_record.id) in sequences:
            foundSequences +=1
            sequences[str(seq_record.id)] = str(seq_record.seq)


print "\n",foundSequences,"reads longer than",minimumReadLength,"nucleotides were found...."


detectionTreshold= foundSequences *detectionTresholdPercentage

dbfile = open(kmerDB)
dbfile.readline()

genotypeKmers = {}
geneGenotypes = {}


outfile = open(outputPrefix+"_kmersCount.txt","w")
outfile.write("Read\tLocus\tGenotype\tNumHits\tReadLength\n")

print "kmer database aquisition...."
while True:
    line = dbfile.readline().rstrip()
    if not line:
        break
    fields = line.split("\t")
    if not fields[0] in geneGenotypes:
        geneGenotypes[fields[0]] = {}
    if not fields[1] in geneGenotypes[fields[0]]:
        geneGenotypes[fields[0]][fields[1]] = []
    for item in fields[2].split(","):
        if not item == '':
            geneGenotypes[fields[0]][fields[1]].append(item)
print "Done!"

print "Searching kmers in reads...."


num = 0
for read in sequences:
    num+=1
    if num%10000 == 0:
        print "Analyzed",num,"sequences"
    query = sequences[read]
    print "Analyzing read",read,"( length:",len(sequences[read]),")"
    for locus in orderedHyperLoci:
        for genotype in geneGenotypes[locus]:
            numKmers = len(geneGenotypes[locus])
            numHits = 0
            for item in geneGenotypes[locus][genotype]:
                if item in query:
                    numHits+=1
            if (float(numHits)/float(numKmers)) > 0:
                print "Found genotype",genotype,"for gene",locus,"with a number of hits",numHits
                outfile.write(read+"\t"+locus+"\t"+genotype+"\t"+str(numHits)+"\t"+str(len(query))+"\n")

outfile.close()

print "mker search finished...."

print "Assemblying the hypervariable genes genotypes"



infile = open(outputPrefix+"_kmersCount.txt")
outfile1 = open(outputPrefix+"_readsCodes.txt","w")
outfile1.write("Read\tGenotypeCode\n")

outfile2 = open(outputPrefix+"_statistics.txt","w")


orderedHyperLoci = ["rl5a","rl6","rl12","rl13","ul1","ul9","ul11","ul20","ul73","ul74","ul120","ul146","ul139"]
genotypeCode = {
    "G1" : "A",
"G10" : "S",
"G1A" : "D",
"G1B" : "F",
"G1C" : "G",
"G2" : "H",
"G2A" : "K",
"G2B" : "L",
"G3" : "Q",
"G3A" : "W",
"G3B" : "E",
"G4" : "R",
"G4A" : "T",
"G4B" : "Y",
"G4D" : "I",
"G5" : "P",
"G6" : "C",
"G7" : "V",
"G8" : "N",
"G9" : "M",
"G11" : "J",
"G4C" : "O",
"G12" : "Z",
"G13" : "X",
"G14" : "B"

}

aa = ["A","S","D","F","G","H","K","L","Q","W","E","R","T","Y","I","P","C","V","N","M","J","O","Z","X","B"]



genotypesCombination = {}
infile.readline() #Reading header

while True:
    line = infile.readline().rstrip()
    if not line:
        break
    fields = line.split("\t")
    if not fields[0] in genotypesCombination:
        genotypesCombination[fields[0]] = {}
    if not fields[1] in genotypesCombination[fields[0]]:
        genotypesCombination[fields[0]][fields[1]] = (genotypeCode[fields[2]], int(fields[3]))
    else:

        if int(fields[3]) > genotypesCombination[fields[0]][fields[1]][1]:
            genotypesCombination[fields[0]][fields[1]] = (genotypeCode[fields[2]], int(fields[3]))
    

#Reconstruct the genotypes succession 
genotypeStringList = []
genotypeStringList_discarded = []
for read in genotypesCombination:
    genotypeString = ""
    genotypesKmers = ""

    for locus in orderedHyperLoci:
        if locus in genotypesCombination[read]:
            genotypeString += genotypesCombination[read][locus][0]
            genotypesKmers += str(genotypesCombination[read][locus][1])+"_"
        else:
            genotypeString += "-"
            genotypesKmers += "-"

    discarded = False
    for a in range(1,len(genotypeString)-1):
        if genotypeString[a] == "-" and not genotypeString[a-1]=="-" and not genotypeString[a+1]=="-":
            #print "Discarding",genotypeString,read
            genotypeStringList_discarded.append(genotypeString)

            discarded = True
    if discarded == False:
        genotypeStringList.append(genotypeString)
        
    print genotypeString,read
    outfile1.write(read+"\t"+genotypeString+"\n")
    #print genotypesKmers
    

#Calculate the frequency of each aminoacid (genotype) per position
outfile2.write("Percentage of each code per position\n")
for a in range(13):
    foundGenotypes = {}
    print "Analyzing gene",orderedHyperLoci[a]
    for b in genotypeStringList:
        for c in aa:
            if c in b[a]:
                if not c in foundGenotypes:
                    foundGenotypes[c] = 0
                foundGenotypes[c] +=1 
    for item in foundGenotypes:
        outfile2.write(item+"\t"+str(foundGenotypes[item])+"\n")


#Perform the graph reconstruction
#Start with the first trimer. 
a=0
trimerCount = {}
for sequence in genotypeStringList:
    if not "-" in sequence[a:a+3]:
        if not sequence[a:a+3] in trimerCount:
            trimerCount[sequence[a:a+3]] = 0
        trimerCount[sequence[a:a+3]] += 1
trimerToElong = set()
for item in trimerCount:
    if trimerCount[item] > detectionTreshold:
        trimerToElong.add(item)

#print "Set to elong"
#print trimerToElong

#Keep on elonging all the trimers to reconstruct the final seuqnece of genotypes
for a in range(1,11):
    trimerCount = {}
    sequencesToElong = set()
    for trimer in trimerToElong:
        overlap = trimer[-2:]
        for sequence in genotypeStringList:
            if not "-"  in sequence[a:a+3]:
                if sequence[a:a+2] == overlap:
                    if not sequence[a:a+3] in trimerCount:
                        trimerCount[sequence[a:a+3]] = 0
                    trimerCount[sequence[a:a+3]] += 1

        for item in trimerCount:
            if trimerCount[item]>detectionTreshold and item[:2]==overlap:
                sequencesToElong.add(trimer+item[2])
    
    trimerToElong = sequencesToElong


#Getting the deposited sequences codes
depCodeFile = open("depositedSequences_codes.txt")
depCodes = {}
while True:
    line = depCodeFile.readline().rstrip()
    if not line:
        break
    fields = line.split("\t")
    if not fields[1] in depCodes:
        depCodes[fields[1]] = fields[0]
depCodeFile.close()
 

print "\n\n"
print "These are the genotypes combinations for the analyzed strains"
outfile2.write("\nFound strains\n")

for item in sequencesToElong: #Sequence to elong now contain the found final strains
    print item
    outfile2.write(item+"\t")
    if item in depCodes:
        outfile2.write(depCodes[item]+"\n")
    else:
        outfile2.write("No match found\n")
    

    
outfile2.write("\nCode trimers frequency in reads\n")
for item in sequencesToElong:
    outfile2.write(item+"\n")
    #print "Analyzing sequence",item
    for a in range(11):
        count = 0
        for sequence in genotypeStringList:
            if item[a:a+3] in sequence[a:a+3]:
                count +=1
        outfile2.write(item[a:a+3]+"\t"+str(count)+"\n")
    outfile2.write("\n")


outfile2.close()
outfile.close()

