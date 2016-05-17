from Bio.SeqUtils import GC
from Bio import SeqIO
import sys

fastq_file = sys.argv[1]
GC_cutoff = sys.argv[2]

low_gc = open("low_gc.fastq","a")
high_gc = open("high_gc.fastq","a")

for rec in SeqIO.parse(fastq_file,"fastq"):
    if GC(rec.seq) <= float(GC_cutoff):
        SeqIO.write(rec,low_gc,"fastq")
    else:
        SeqIO.write(rec,high_gc,"fastq")

low_gc.close()
high_gc.close()

