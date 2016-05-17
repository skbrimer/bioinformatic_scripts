from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio import SeqIO
import argparse 

## this makes a dict of the samtools depth coverage input file of all the 0's
coverageDict = {} ## first the empty dictionary

## this loops over the input depth information and appends the dictionary with 
## a key,value as genome position, depth of coverage
coverage = open("/home/sbrimer/Desktop/Columbia ST/1819_1/coverage.txt","r")
for line in coverage:
    coverages = line.strip().split('\t') #added strip, values also had newline character 
    coverageDict[int(coverages[1])] = int(coverages[2])
coverage.close()

## making a set of only the index of the missing
missing = {index for index,value in coverageDict.items() if value == 0}

def filter_coverage(index,value):
        if index in missing:
            return "N"
        else:
            return value 

## this should read in the sequence file as an oblect, change it to a mutable
## object, then anywhere the coverageDict value = 0. Replace that nt with "N"
append_handle = open("NewFile.fasta","a")
for seq_record in SeqIO.parse("Salmonella enterica subsp.fasta","fasta"):
    newSeq= "".join(filter_coverage(index,value)for index,value in enumerate(seq_record.seq, start= 1))
    SeqIO.write(SeqRecord(Seq(newSeq),id="something",description="something_else"),append_handle,"fasta")
append_handle.close()
