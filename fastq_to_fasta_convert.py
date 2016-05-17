# this script is used to convert fastq files to fasta files 
# then to rename the fasta ID with the sample ID from the lab

from Bio import SeqIO
import sys 

# grabbing the file and the name 
seq_file = sys.argv[1]
labels = seq_file.split(".")

# converting the file from fastq to fasta
SeqIO.convert(seq_file,"fastq",labels[0]+".fasta","fasta")

# taking the converted file and then changing the fasta header
handle = open(labels[0]+".fasta","a")

for seq_record in SeqIO.parse(handle,"fasta"):
    old_header = seq_record.id
    new_header = labels[0]
    seq_record.id = new_header + "_" + old_header # renaming the pseudogene with
                                                  # the lab id and the referance 
                                                  # used
    seq_record.description = "" # this strips the old header out
    SeqIO.write(seq_record, handle,"fasta")

handle.close()


 
