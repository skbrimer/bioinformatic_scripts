#!/bin/bash
###############################################################################
#####
#####  A Pipeline for training
#####  $0 <fastq directory> for convenience,
#####  Fastq file extension is assumed to be 'fastq'
#####
###############################################################################


###############################################################################
#####
#####  Set up locations
#####  Unnecessary if in $PATH, to test use 'which'
#####  Check write permission
#####  
###############################################################################

FQ_DIR=$1
REF=/home/byoo/proj/ceva/ref/AM933172.fa
NUM_THREAD=10
FASTQC=/illumina/software/FastQC/fastqc
TMAP=/illumina/software/iontorrent/TMAP/tmap
SAMSTATS=/illumina/software/ea-utils.1.1.2-484/sam-stats
SAMSTAT=/illumina/software/samstat/src/samstat
FREEBAYES=/home/byoo/sw/freebayes/bin/freebayes
SNPEFF=/home/byoo/sw/snpEff
TMAP_OUT_DIR=${FQ_DIR}/../tmap
FastQC_OUT_DIR=${FQ_DIR}/../FastQC
FREEBAYES_OUT_DIR=${FQ_DIR}/../freebayes
SNPEFF_OUT_DIR=${FQ_DIR}/../snpEff

: '
###############################################################################
#####
#####  Running Fastq file quality assessment using FastQC
#####  Java should be in PATH or you need to specify --java
#####  
###############################################################################

echo "Running Fastq file quality assessment using FastQC"
mkdir $FastQC_OUT_DIR
nohup /usr/bin/time --verbose $FASTQC -o $FastQC_OUT_DIR -t $NUM_THREAD ${FQ_DIR}/*.fastq &> $FastQC_OUT_DIR/FastQC.log &






###############################################################################
#####
#####  Mapping Sequences using TMAP
#####  Genome is assumed to be prepared and indexed
#####  
###############################################################################

echo "Mapping Sequences using TMAP"
GLOBAL_OPTION="-f ${REF} -i fastq -n $NUM_THREAD"
mkdir $TMAP_OUT_DIR
find ${FQ_DIR} -name "*.fastq" | while read FILENAME;
do
	BASENAME=`basename $FILENAME`
	/usr/bin/time --verbose $TMAP mapall -r $FILENAME -s ${TMAP_OUT_DIR}/${BASENAME}.sam $GLOBAL_OPTION -v stage1 map1 map2 map3 &> ${TMAP_OUT_DIR}/${BASENAME}.log
done



###############################################################################
#####
#####  Mapping Quality Assessment using Sam-stats
#####  
###############################################################################

echo "Mapping Quality Assessment using Sam-stats"
find ${TMAP_OUT_DIR} -name "*.sam" | while read FILENAME;
do
	BASENAME=`basename $FILENAME`
	nohup /usr/bin/time --verbose $SAMSTATS -D ${TMAP_OUT_DIR}/${BASENAME} &> ${TMAP_OUT_DIR}/${BASENAME}.stats &
done



###############################################################################
#####
#####  Mapping Quality Assessment using Samstat
#####  
###############################################################################

echo "Mapping Quality Assessment using Samstat"
nohup /usr/bin/time --verbose $SAMSTAT ${TMAP_OUT_DIR}/*.sam &> ${TMAP_OUT_DIR}/samstat.log &


###############################################################################
#####
#####  Variant Calling using freebayes
#####  
###############################################################################

echo "Variant Calling using freebayes"
mkdir $FREEBAYES_OUT_DIR
find ${TMAP_OUT_DIR} -name "*.sam" | while read FILENAME;
do
	BASENAME=`basename $FILENAME`
	samtools view -bS $FILENAME > ${TMAP_OUT_DIR}/${BASENAME%.*}.bam
	samtools sort ${TMAP_OUT_DIR}/${BASENAME%.*}.bam ${TMAP_OUT_DIR}/${BASENAME%.*}.sorted
	samtools index ${TMAP_OUT_DIR}/${BASENAME%.*}.sorted.bam
	/usr/bin/time --verbose $FREEBAYES -p 1 -f ${REF} -v ${FREEBAYES_OUT_DIR}/${BASENAME%.*}.vcf ${TMAP_OUT_DIR}/${BASENAME%.*}.sorted.bam &> ${FREEBAYES_OUT_DIR}/${BASENAME%.*}.freebayes.log
done

sed -i -e 's/AM933172\t/NC_011294\t/g' $FREEBAYES_OUT_DIR/*.vcf
'

###############################################################################
#####
#####  Variant Annotation using Snpeff
#####  java -jar ./snpEff.jar download Salmonella_enterica_serovar_Enteritidis_P125109_uid59247
#####  
###############################################################################

echo "Variant Calling using freebayes"
mkdir $SNPEFF_OUT_DIR
find ${FREEBAYES_OUT_DIR} -name "*.vcf" | while read FILENAME;
do
	BASENAME=`basename $FILENAME`
	/usr/bin/time --verbose java -Xmx4g -jar $SNPEFF/snpEff.jar eff -c $SNPEFF/snpEff.config -v Salmonella_enterica_serovar_Enteritidis_P125109_uid59247 -s "snpEff_summary_${BASENAME}.html" ${FILENAME} 1> ${SNPEFF_OUT_DIR}/${BASENAME%.*}.eff.vcf 2> ${SNPEFF_OUT_DIR}/${BASENAME}.log
done
SNPEFF_SUM_OUT_DIR=${FQ_DIR}/../snpEff_summary
mkdir $SNPEFF_SUM_OUT_DIR
mv ./snpEff_summary_* $SNPEFF_SUM_OUT_DIR
