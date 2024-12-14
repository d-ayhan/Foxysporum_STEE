# load software
module load bwa/0.7.15
module load picard/2.0.1
module load samtools/1.4.1
module load bedtools/2.26.0

set -e

sample="Fol4287"
fastqdir="reads/"
picard="/picard.jar"
ref="/Fol4287_GCA_003315725_genomic.fna"
RG="@RG\tID:"$sample"\tPL:illumina\tLB:LIB-"$sample"\tSM:"$sample

bwa index $ref
for A in H3K27me3_1 H3K27me3_2 H3K4me2_1 H3K4me2_2
do

# map reads to reference
date
echo "bwa mem"
bwa mem -t 8 $ref $fastqdir$A".fastq" > $A".sam"

# Clean and fix SAM
date
echo "Picard CleanSam"
java -Xmx15g -jar $picard CleanSam I=$A".sam" O=$A".clean.sam"

rm $A".sam"

# SAM>BAM
date
echo "Samtools View"
samtools view -b -T $ref -o $A".bam" $A".clean.sam"

rm $A".clean.sam"

samtools sort -o $A".sorted.bam" $A".bam"
rm $A."bam"

# Mark duplicate reads
date
echo "Picard MarkDuplicates"
java -Xmx9g -jar $picard MarkDuplicates I=$A".sorted.bam" O=$A".sorted.dedup.bam" M=$A".sorted.dedup.metrics" REMOVE_DUPLICATES=true

# index BAM
date
echo "Samtools Index"
samtools index $A".sorted.dedup.bam"

date
echo "flagstat "$A".sorted.dedup.bam"
samtools flagstat $A".sorted.dedup.bam"

# get genome coverage
date
echo "genomecov"
bedtools genomecov -d -ibam $A".sorted.dedup.bam" -g $ref".fai" > $A".genomecov.txt"
done