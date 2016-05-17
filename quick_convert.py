# this script is used to convert fastq files to fasta files 

from Bio import SeqIO
import sys 

# grabbing the file and the name 
seq_file = sys.argv[1]
labels = seq_file.split(".")

# converting the file from fastq to fasta
SeqIO.convert(seq_file,labels[1],labels[0]+".fasta","fasta")
