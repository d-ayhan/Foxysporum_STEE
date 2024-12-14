###################### 
map() {
local sample=$1
A=${sample:0:2}
RG="@RG\tID:"$sample"\tPL:illumina\tLB:LIB-"$sample"\tSM:"$sample

# map reads to reference
date
echo "bwa mem"
bwa mem -t 8 -M -a -v 1 -R $RG $ref $fastqdir$sample"_R1_001.fastq" $fastqdir$sample"_R2_001.fastq" > $A".sam"

# Clean and fix SAM
date
echo "Picard CleanSam"
java -Xmx9g -jar $picard CleanSam I=$A".sam" O=$A".clean.sam"

rm $A".sam"

# SAM>BAM
date
echo "Samtools View"
samtools view -b -T $ref -o $A".bam" $A".clean.sam"

rm $A".clean.sam"

# Fixes paired reads and sorts
date
echo "Picard FixMateInformation"
java -Xmx9g -jar $picard FixMateInformation SO=coordinate I=$A".bam" O=$A".fixmate.bam"

rm $A."bam"

# Mark duplicate reads
date
echo "Picard MarkDuplicates"
java -Xmx9g -jar $picard MarkDuplicates I=$A".fixmate.bam" O="bam/"$A".fixmate.dedup.bam" M="dedup_metrics/"$A".fixmate.dedup.metrics"

rm $A".fixmate.bam"

# index BAM
date
echo "Samtools Index"
samtools index "bam/"$A".fixmate.dedup.bam"

date
echo "flagstat "$A".fixmate.dedup.bam"
samtools flagstat "bam/"$A".fixmate.dedup.bam" > "flagstat/"$A".txt"

# Get coverage statistics 
date
echo "Picard CollectRawWgsMetrics"
java -Xmx18g -jar $picard CollectRawWgsMetrics I="bam/"$A".fixmate.dedup.bam" O="wgs_metrics/"$A".raw_wgs_metrics" R=$ref

# Plot insert size histogram
date
echo "Picard CollectInsertSizeMetrics"
java -jar $picard CollectInsertSizeMetrics I="bam/"$A".fixmate.dedup.bam" O="insert_metrics/"$A".insertMetrics.txt" H="insert_metrics/
"$A".insertHistogram.pdf" M=0.5 ASSUME_SORTED=false

date
echo "FreeBayes"
$freebayes -f $ref -C 3 -p 1 "bam/"$A".fixmate.dedup.bam" > "fb_vcf/"$A".fb.vcf"

date
echo "GRIDSS"
java -Xmx31g -jar $gridss OUTPUT="SV/"$A".SV.vcf" ASSEMBLY="SV/"$A".SV.assembly.bam" INPUT="bam/"$A".fixmate.dedup.bam" \
REFERENCE_SEQUENCE=$ref TMP_DIR=/home/da87a/TEMP/ WORKING_DIR=/home/da87a/TEMP/

# get genome coverage
date
echo "genomecov"
bedtools genomecov -d -ibam "bam/"$A".fixmate.dedup.bam" -g $ref".fai" > "genomecov/"$A".genomecov.txt"
}
#############################
bqrecal() {
local A=$1
local gatk="bin/GATK/3.5/GenomeAnalysisTK.jar"

echo "$A is initialized."

# Calibrate base quality scores
java -Xmx10g -jar $gatk -T BaseRecalibrator -nct 10 \
-R $ref \
-I "bam/"$A".fixmate.dedup.bam" \
--knownSites "filtered_variants.vcf" \
-o $A".fixmate.dedup.recal.table"

# Re-create BAM file with calibrated bases
java -Xmx10g -jar $gatk -T PrintReads \
-R $ref \
-I "bam/"$A".fixmate.dedup.bam" -nct 10 \
-BQSR $A".fixmate.dedup.recal.table" \
-o "bam/"$A".fixmate.dedup.recal.bam"  -filterNoBases

# Post re-calibration
java -Xmx10g -jar $gatk -T BaseRecalibrator \
-R $ref \
-I "bam/"$A".fixmate.dedup.bam" \
--knownSites "filtered_variants.vcf" \
-BQSR $A".fixmate.dedup.recal.table" \
-o $A".fixmate.dedup.recal_after.table" -nct 10

# Plot calibration results
java -Xmx10g -jar $gatk -T AnalyzeCovariates \
-R $ref -before $A".fixmate.dedup.recal.table" \
-after $A".fixmate.dedup.recal_after.table" \
-plots $A".recalQC.pdf"

echo "$A is done."
}
###################

set -e
# load software
module load bwa/0.7.15
module load picard/2.0.1
module load R/3.2.2
module load samtools/1.4.1
module load bedtools/2.26.0
module load gnuplot/4.6.5
picard="/share/pkg/picard/2.0.1/picard.jar"
freebayes="bin/freebayes"
gatk="bin/gatk-4.1.4.1/gatk"
gridss="bin/gridss-1.4.1-jar-with-dependencies.jar"

ref="Fol4287_GCA_003315725_genomic.fna"
fastqdir="illuminareads/"

###################
# map read and get stats
for run in "P1-10_S1" "P2-10_S2" "P3-10_S3" "P4-10_S4" "P5-10_S5"; do map "$run" & done
wait
for run in "Y1-10_S6" "Y2-10_S7" "Y3-10_S8" "Y4-10_S9" "Y5-10_S10"; do map "$run" & done
wait
for run in "M1-10_S12" "M2-10_S13" "M3-10_S14" "M4-10_S15" "M5-10_S16"; do map "$run" & done
wait
for run in "P2-1" "P2-2" "P2-3" "P2-4" "P2-5" "P2-6" "P2-7" "P2-8" "P2-9" "Y3-1" "Y3-5" "Y3-9" "M4-1" "M4-5" "M4-9"; do foo "$run" & done
wait

###################
# recalibrate base quality scores
for run in P1 P2 P3 P4 P5 Y1 Y2 Y3 ; do bqrecal "$run" & done
wait
for run in Y4 Y5 M1 M2 M3 M4 M5 ; do bqrecal "$run" & done
wait

###################
# SNV calling
A=".fixmate.dedup.recal.bam"

$gatk --java-options -Xmx45g Mutect2 -R $ref \
-I "bam/WT"$A \
-I "bam/P1"$A \
-I "bam/P2"$A \
-I "bam/P3"$A \
-I "bam/P4"$A \
-I "bam/P5"$A \
-normal "WT_S11" \
-O somatic.P.vcf.gz

$gatk --java-options -Xmx45g Mutect2 -R $ref \
-I "bam/WT"$A \
-I "bam/Y1"$A \
-I "bam/Y2"$A \
-I "bam/Y3"$A \
-I "bam/Y4"$A \
-I "bam/Y5"$A \
-normal "WT_S11" \
-O somatic.Y.vcf.gz

$gatk --java-options -Xmx45g Mutect2 -R $ref \
-I "bam/WT"$A \
-I "bam/M1"$A \
-I "bam/M2"$A \
-I "bam/M3"$A \
-I "bam/M4"$A \
-I "bam/M5"$A \
-normal "WT_S11" \
-O somatic.M.vcf.gz

$gatk --java-options -Xmx45g FilterMutectCalls -R $ref -V somatic.P.vcf.gz -O somatic.P_filtered.vcf.gz
$gatk --java-options -Xmx45g FilterMutectCalls -R $ref -V somatic.Y.vcf.gz -O somatic.Y_filtered.vcf.gz
$gatk --java-options -Xmx45g FilterMutectCalls -R $ref -V somatic.M.vcf.gz -O somatic.M_filtered.vcf.gz


A=".fixmate.dedup.bam"

$gatk --java-options -Xmx45g Mutect2 -R $ref \
-I $WT \
-I "bam/P2-1"$A \
-I "bam/P2-2"$A \
-I "bam/P2-3"$A \
-I "bam/P2-4"$A \
-I "bam/P2-5"$A \
-I "bam/P2-6"$A \
-I "bam/P2-7"$A \
-I "bam/P2-8"$A \
-I "bam/P2-9"$A \
-I "bam/P2.fixmate.dedup.recal.bam" \
-normal "WT_S11" \
-O somatic.P2_inter.vcf.gz

$gatk --java-options -Xmx45g Mutect2 -R $ref 
-I $WT \
-I "bam/Y3-1"$A \
-I "bam/Y3-5"$A \
-I "bam/Y3-9"$A \
-I "bam/Y3.fixmate.dedup.recal.bam" \
-normal "WT_S11" \
-O somatic.Y3_inter.vcf.gz

$gatk --java-options -Xmx45g Mutect2 -R $ref \
-I $WT \
-I "bam/M4-1"$A \
-I "bam/M4-5"$A \
-I "bam/M4-9"$A \
-I "bam/M4.fixmate.dedup.recal.bam" \
-normal "WT_S11" \
-O somatic.M4_inter.vcf.gz

$gatk --java-options -Xmx45g FilterMutectCalls -R $ref -V somatic.P2_inter.vcf.gz -O somatic.P2_inter.filtered.vcf.gz
$gatk --java-options -Xmx45g FilterMutectCalls -R $ref -V somatic.Y3_inter.vcf.gz -O somatic.Y3_inter.filtered.vcf.gz
$gatk --java-options -Xmx45g FilterMutectCalls -R $ref -V somatic.M4_inter.vcf.gz -O somatic.M4_inter.filtered.vcf.gz
