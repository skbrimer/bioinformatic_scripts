''' file conversion for from genbank to fasta and to make a bedfile as well. 
    Code modified from Brant Faircloth's "bed_from_genbank.py" you can find the 
    original code here https://gist.github.com/brantfaircloth/893580 '''

from Bio import SeqIO
import sys

# Lets gets some files 
gb_file = sys.argv[1]
lables = gb_file.split(".")

# Producing the bed file, first make a new writeable file and creat header
# Added the line break \n for a cleaner look with a text editor, sorry my OCD. 

out_bedfile = open(lables[0]+".bed", 'w')
bed_header = "track name="+ '"'+lables[0]+'"'+" description="+'"'+lables[0] +" genes"+'"'+" itemRgb=On"+"\n"
out_bedfile.write(bed_header)

# loop will work with multi genbank files. 

for record in SeqIO.parse(open(gb_file,"rU"),"genbank"):
    for feature in record.features:
        if feature.type == 'CDS':
            start = int(feature.location.start)
            stop = int(feature.location.end)
            try:
                name = feature.qualifiers['gene'][0]
                #name = feature.qualifiers['product'][0]
            except:
                #some features only have locus tags
                name = feature.qualifiers['locus_tag'][0]
            if feature.strand < 0:
                strand = "-"
            else:
                strand = "+"
            bed_line = record.id +"\t{0}\t{1}\t{2}\t500\t{3}\t{0}\t{1}\t50,205,50\n".format(start, stop, name,strand)
            out_bedfile.write(bed_line)

out_bedfile.close()

#convert to fasta - so easy (thank you biopython)
SeqIO.convert(gb_file,"genbank",lables[0]+".fasta","fasta")







