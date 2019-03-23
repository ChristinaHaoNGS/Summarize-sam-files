# Summarize-sam-files
Summarizes sam files from amplicon-seq to detect single nucleotide variations (SNV)

 ####APPLICATION####
This script is useful for detecting sequence variations from amplicon-seq data. An example application is as follows: PCR amplification using loci specific primers is performed to determine the mutations (SNV and indels) present in those regions. The PCR product is sequenced, and Bowtie is used to map reads to expected amplicons provided in the FASTA file. This script summarizes the SAM file, reporting the number of times each mapped and unmapped reads occur. 

####EXAMPLE OUTPUT####

Example summary of mapped reads: shown below are top three most abundant reads mapping to “ALK-3” amplicon defined in the Fasta file. The first sequence occurs 1965 times, and it maps exactly to the reference sequence. The second most abundant read occurs 52 times, and it has a mismatch at position 22, where a G in the reference is replaced by an A. The length of the sequence that is analyzed is user defined (options -s and -e). As default, the script analyzes 60 bases from positions 10 to 70. 

ALK-3 CATGGAAGCCCTGATCATCAGGTAAAGCCACAGAGAGACACCCTCACCCCAACTCCCCTC 1965	

ALK-3	CATGGAAGCCCTGATCATCAGATAAAGCCACAGAGAGACACCCTCACCCCAACTCCCCTC	52 22:G>A

ALK-3	CATGGAAGCCCTGATCATCAGGTGAAGCCACAGAGAGACACCCTCACCCCAACTCCCCTC 45 24:A>G

Example summary of unmapped reads: the script summarizes reads that do not map to expected sequences, as shown below. Only reads that occur above a user specified threshold (option -u) is reported. The default -u value is 1000.

CCTTCTGCACCTTGCTGCCTGAACAGCTCCCTCATCGGCTTTCCCCTGTGGTCCCTGGGG	23433
AGCTCATCAAACGGAGGGCTGCTTTTTCCTGATCTCTTTGAGCACCGTGTTTGGATTAGC	23714
GAGCCAATATTGTCTTTGTGTTCCCGGACATAGTCCAGGAGGCAGCCGAAGGGCATGAGC	26434

####REQUIRED INPUT FILES####

1)	SAM file: output from Bowtie or other mapping algorithm
2)	Fasta file of reference sequences 

####HOW TO RUN####

Place fullseq_exam.py, sam file, fasta file in the same directory. To run type: 
python fullseq_exam.py -f your_file.sam -a your_file.fa 

Arguments

-f, --samfile	[required] sam file
-a, --fastafile	[required] fasta file
-s, --start	[optional] Position of read to begin analysis. Default: 10 
-e, --end	[optional] Position of read to end analysis. Default: 70
-c, --summarycutoff	[optional] In the summary of mapped reads, report only reads with greater than -c reads. Default: 10
-u, --unmappedcutoff	[optional] In the summary of unmapped reads, report only reads with greater than -u reads. Default: 100

