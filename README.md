# minion_Genotyper

## Main algorithm
For 13 HCMV hypervariable genes, this tool search genotype-specific motifs in all the sequeneced reads whose length exceed a user-defined threshold (see section below "how to run the software"). Genes are ordered according to their position along the HCMV genome. For each read, if a genotype is detected by a specific motif, a letter is used to code it within a 13-letter string, hereafter referred to as a genotype string. If none of the genotype-specific motifs are found in the read, a hyphen symbol is inserted in the  genotype string. As a way of example, if the first two hypervariable genes are not found in the read, then 4 genotypes are detected for the following 4 hypervariable genes and no genotype is detected for the remaining ones the reported genotype string would look like the following:

-ACGA-------

A read resulting in genotyping string in which a hyphen symbol separates two coded genotypes, is considered unreliable and therefore discarded. Similarly, if at any give position of the genotyping string, a three consecutive letters are lowly represented in all the other read genotyping strings, the read is discarded (see section "output files")
After all the genotyping strings are computed, they are aligned and complete haplotypes are reconstructed by exploting the overlap between consecutive trimers in a fashion recalling the Overlap Layout Consensus (OLC) method, as shown in the following example:

ACCF---------    Genotype string 1
-CCFAA-------    Genotype string 2
--CFAAGFGG---    Genotype string 3
-------FGGADA    Genotype string 4

ACCFAAGFGGADA    Final haplotype


The found haplotype are compared to those of 244 deposited HCMV genes. 


## Folder content
minion_Genotyper.py: 		This is the main program

kmerDB folder:			This is the folder where the kmers that are used for the
genotyping are kept 

depositedSequences_codes.txt	This is a file were the genotyping codes for all the 244 
					deposited genomes are reported


example_kmersCount.txt		This an example of how the kmer count will be reported for 
					each read

example_readsCodes.txt		This is an example of how the reads are coded 

example_statistics.txt		This is an example of the main output file (see below)

codeTable.xls				This is the table that links the used one letter code with the					original genotype codes.

## How to run the software
The software runs with the following command within the minion_Genotyper folder:

python minion_Genotyper.py  inputFileName minimumReadLength outPrefix

inputFileName: 		fastq file name with its complete path
minimumReadLength;	Minimum length for a read to be included in the downstream
				analyses (so far 50000 proved to be a good number)
outPrefix:			this is a name that will constitute the prefix of all the output files (e.g. 
				the outPrefix example will generate the output files that are now 
				present in the folder

 
## Output file
The file _statistics.txt is the main output folder. 

In the first section you will find for each hypervariable gene the code of the genotype that has been found and how many times this occurred. 

In the second section you will find the reconstructed strains that are reported as the genotypes string code. If the string code has been previously observed in any of the 244 deposited HCMV genomes, then the name of such genome will be reported close to the coded strain. If not, a “No match found” will be reported. 

In the third section some “reliability” statistics are reported. Each found strain code is divided in trimers and all the code is scanned by using a sliding window (window size = 3, window step = 1). The software will count how many times that trimer has been called in the reads. Strain codes featuring some of the trimers with low count as compared to the other observed values are probably artefact and possibly need to be removed.  As a way of example, in the file example_statistics.txt all strain codes have at least one trimer that is present twice or less in the reads, with the exception of the last two (which in this case are also the only ones to get a match in the 244 deposited genomes as seen in the section before). The number of unreliable strain codes should be very low and it is not in this case just because I used a very low number of reads to test the software. You can use this section to copy and paste the trimer counts for each strain code into Excel and calculate for each trimer position the frequency and then average all of them to obtain the frequency of each strain in a mixed strain infection.
The sensitivity of the software can be tuned by changing the following line in the main software

detectionTresholdPercentage = 0.002

If you low down that 0.002 the sensitivity will increase, but of course in this case you may risk to take in some low frequent / possibly artifactual strains.


