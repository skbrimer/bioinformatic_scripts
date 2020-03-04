import sys
import csv
import pandas as pd

## CLI passed file of interest, tab deliminated

BLAST_File = sys.argv[0]

## read in the file

df1 = pd.read_csv(BLAST_File, sep="\t", header = None)

## create a unique list of the nodes/contigs

nodeList = {node for node in list(df1[0])}

## A list to store just the max for each node

topHits = []
   
for node in nodeList:
    contig = df1[0] == node            #looks up each node from the nodelist
    hits = list(df1[contig].max())     #pulls the max value for that node 
                                       #and writes it as a list
    topHits.append(hits)
    
with open("topHits.csv", "w+") as f:
    writer = csv.writer(f)
    writer.writerows(topHits)

