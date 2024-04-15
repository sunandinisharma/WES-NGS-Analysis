#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=168:00:00
#SBATCH --mem=40GB
#SBATCH --job-name=hg38WESPipeline
#SBATCH --licenses=common
#SBATCH --workdir=/lustre/work/javeediqbal/shared/MCL_Sunandini/MCL_WES

######## Sunandini Sharma 01302020 #############
######## Weiwei Zhang 01302020 ########
######## AlyssaBouska 022020 Mod ########
######## Start Variables ########

module load compiler/gcc/4.9 java/1.8 R/3.3 bwa/0.7 samtools/1.3 picard/2.9 gatk4/4.1 varscan/2.4
r1=${1}
r2=${2}
r3=${3}
r4=${4}
pfx=${5}

ref="/common/javeediqbal/abouska/reference_GRCh38_tophat/Homo_sapiens/NCBI/GRCh38/Sequence/BWAIndex/version0.6.0/genome.fa"
indels1="/common/javeediqbal/abouska/Sequencing_Pipeline/GATK.resources/GRCh38/Mills_and_1000G_gold_standard.indels.hg38.vcf"
indels2="/common/javeediqbal/abouska/Sequencing_Pipeline/GATK.resources/GRCh38/Homo_sapiens_assembly38.known_indels.vcf"
dbsnp="/work/javeediqbal/abouska/Sequencing_Pipeline/GATK.resources/GRCh38/dbsnp_146.hg38.vcf"
BAIT="/common/javeediqbal/abouska/Sequencing_Pipeline/SureSelectAllExonV7_S31285117_hs_hg38/SureSelectAllExonV7_hg38_S31285117_Covered.interval.list"
TARGET="/common/javeediqbal/abouska/Sequencing_Pipeline/SureSelectAllExonV7_S31285117_hs_hg38/SureSelectAllExonV7_hg38_S31285117_Padded.interval.list"
REGIONS="/common/javeediqbal/abouska/Sequencing_Pipeline/SureSelectAllExonV7_S31285117_hs_hg38/SureSelectAllExonV7_hg38_S31285117_Regions.interval.list"
MERGEDPROBES="/common/javeediqbal/abouska/Sequencing_Pipeline/SureSelectAllExonV7_S31285117_hs_hg38/SureSelectAllExonV7_hg38_S31285117_MergedProbes.interval.list"


############info on Capture Design Files from Agilent#################################
#regions=Targeted exon intervals, curated and targeted by Agilent Technologies
#covered=Merged probes and sequences with 95% homology or above
#Merged probes and sequences with 95% homology or above extended 100 bp at each side
#######################################################################################

annovar="/common/javeediqbal/abouska/annovar/"
hsMetrics="/common/javeediqbal/abouska/Sequencing_Pipeline/CalculateHsMetrics.jar"


memoryi="512m"
memorym="40g"
threads="16"

######## End Variables ########

######## Alignment ########

mkdir $pfx $pfx/Alignment $pfx/Variants $pfx/Annotation $pfx/metrics $pfx/metrics/CollectHsMetrics $pfx/metrics/CollectInsertSizeMetrics $pfx/metrics/CalculateHsMetrics $pfx/metrics/samtoolsstats 

echo "[`date`] Job starting"

###for reads 70bp-1Mbp start###
echo "[`date`] Read aligning for $pfx using BWA mem"
bwa mem -t $threads $ref $r1 $r2 > $pfx/Alignment/$pfx.L1.sam &&
bwa mem -t $threads $ref $r3 $r4> $pfx/Alignment/$pfx.L2.sam &&
###for reads 70bp-1Mbp end###

###for reads <70bp start###
#echo "[`date`] Read aligning for $pfx using BWA aln"
#$bwa aln -t $threads $ref $r1 > $pfx/Alignment/$pfx.r1sai
#$bwa aln -t $threads $ref $r2 > $pfx/Alignment/$pfx.r2sai
#echo "[`date`] Pair-end mapping for $pfx using BWA sampe"
#$bwa sampe $ref $pfx/Alignment/$pfx.r1sai $pfx/Alignment/$pfx.r2sai $r1 $r2 -f $pfx/Alignment/$pfx.sam
###for reads <70bp end###

echo "[`date`] Transforming SAM file to BAM file using Samtools view"
samtools view -Sb $pfx/Alignment/$pfx.L1.sam > $pfx/Alignment/$pfx.L1.bam &&
samtools view -Sb $pfx/Alignment/$pfx.L2.sam > $pfx/Alignment/$pfx.L2.bam &&
rm $pfx/Alignment/$pfx.L1.sam
rm $pfx/Alignment/$pfx.L2.sam


echo "[`date`] Sorting the BAM file using Samtools sort"
samtools sort $pfx/Alignment/$pfx.L1.bam -o $pfx/Alignment/$pfx.L1.sort.bam &&
samtools sort $pfx/Alignment/$pfx.L2.bam -o $pfx/Alignment/$pfx.L2.sort.bam &&

echo "[`date`] Indexing the sorted BAM file using Samtools index"
samtools index $pfx/Alignment/$pfx.L1.sort.bam &&
samtools index $pfx/Alignment/$pfx.L2.sort.bam &&

###Add or Replace read group for each Lane start###
echo "[`date`] Fixing Read Group in the bam file using PicardTools AddOrReplaceReadGroups.jar for $pfx."
picard AddOrReplaceReadGroups I=$pfx/Alignment/$pfx.L1.sort.bam O=$pfx/Alignment/$pfx.L1.fixed_RG.bam SO=coordinate RGID=$pfx.L1 RGLB=$pfx RGPL=illumina RGPU=$pfx RGSM=$pfx VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true &&

picard AddOrReplaceReadGroups I=$pfx/Alignment/$pfx.L2.sort.bam O=$pfx/Alignment/$pfx.L2.fixed_RG.bam SO=coordinate RGID=$pfx.L2 RGLB=$pfx RGPL=illumina RGPU=$pfx RGSM=$pfx VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true &&


echo "[`date`] Merging the BAM file using Samtools merge"
samtools merge $pfx/Alignment/$pfx.all.bam $pfx/Alignment/$pfx.L1.fixed_RG.bam $pfx/Alignment/$pfx.L2.fixed_RG.bam

echo "[`date`] Sorting the BAM file using Samtools sort"
samtools sort $pfx/Alignment/$pfx.all.bam -o $pfx/Alignment/$pfx.all.sort.bam
echo "[`date`] Indexing the sorted BAM file using Samtools index"
samtools index $pfx/Alignment/$pfx.all.sort.bam

######## PostAlignment ########

#java -Xmx${memory}g -jar $gatk --help
echo "PostAlignment For $pfx"
echo "[`date`] Job starting"

###Mark duplicates start###
echo "[`date`] Marking Duplicates using PicardTools MarkDuplicates.jar for $pfx."
picard MarkDuplicates INPUT=$pfx/Alignment/$pfx.all.sort.bam OUTPUT=$pfx/Alignment/$pfx.dup_rmv.bam M=$pfx/Alignment/$pfx.metric VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true && 
#SORTING_COLLECTION_SIZE_RATIO=0.125, default is 0.25, if run out of memory then reduce the number
#REMOVE_DUPLICATES=false is default so do not need to add in

echo "[`date`] Sorting the bam file using Samtools sort for $pfx."
samtools sort $pfx/Alignment/$pfx.dup_rmv.bam -o $pfx/Alignment/$pfx.dup_rmv.sort.bam &&

echo "[`date`] Indexing $pfx.recalibrated.bam using samtools index for $pfx."
samtools index $pfx/Alignment/$pfx.dup_rmv.sort.bam &&


###Recalibrate Bases start###
echo "[`date`] Create recal_data.1.table using BaseRecalibrator for $pfx." 
gatk --java-options "-Xms${memoryi} -Xmx${memorym}" BaseRecalibrator -R $ref -I $pfx/Alignment/$pfx.dup_rmv.sort.bam --known-sites $dbsnp --known-sites $indels1 --known-sites $indels2 -O $pfx/Alignment/$pfx.recal_data.1.table &&

echo "[`date`] Recalibrate $pfx.indel_realigned.bam using PrintReads for $pfx."
gatk --java-options "-Xms${memoryi} -Xmx${memorym}" ApplyBQSR -I $pfx/Alignment/$pfx.dup_rmv.sort.bam -bqsr $pfx/Alignment/$pfx.recal_data.1.table -O $pfx/Alignment/$pfx.recalibrated.final.bam &&

##BQSR the second pass##
echo "[`date`] Create recal_data.2.table using BaseRecalibrator for $pfx."
gatk --java-options "-Xms${memoryi} -Xmx${memorym}" BaseRecalibrator -R $ref -I $pfx/Alignment/$pfx.recalibrated.final.bam --known-sites $dbsnp --known-sites $indels1 --known-sites $indels2 -O $pfx/Alignment/$pfx.recal_data.2.table &&

echo "[`date`] Recalibrate $pfx.recalibrated.bam using ApplyBQSR for $pfx."
gatk --java-options "-Xms${memoryi} -Xmx${memorym}" ApplyBQSR -I $pfx/Alignment/$pfx.recalibrated.final.bam -bqsr $pfx/Alignment/$pfx.recal_data.2.table -O $pfx/Alignment/$pfx.recalibrated.pass2.bam &&

##Compare before and after outputs##
gatk --java-options "-Xms${memoryi} -Xmx${memorym}" AnalyzeCovariates -before $pfx/Alignment/$pfx.recal_data.1.table -after $pfx/Alignment/$pfx.recal_data.2.table -plots $pfx/Alignment/$pfx.BQSR.pdf &&

#######Submit Mutect #################################################################

sbatch -o ${pfx}/err.out.files/${pfx}.mutect2.out -e ${pfx}/err.out.files/${pfx}.mutect2.err GATK4.MUTECT2.single.sample.FLAGGEDdbSNPfiltered.GDC.PoN_latestAnnovar.092320.sh $pfx/Alignment/$pfx.recalibrated.final.bam ${pfx}

################################################################

#######Submit CW #################################################################

norm=TGEN_PBMC_13_combined
sbatch -o ${pfx}/err.out.files/${pfx}.submitCW.out -e ${pfx}/err.out.files/${pfx}.submitCW.err submit.CW.sh ${pfx} ${norm}

################################################################



####### metrics #######

bam="$pfx/Alignment/$pfx.recalibrated.final.bam"

#picard CollectHsMetrics
picard -Xms4096m -Xmx8g CollectHsMetrics I=$bam O=$pfx/metrics/CollectHsMetrics/$pfx.metrics R=$ref BAIT_INTERVALS=$MERGEDPROBES TARGET_INTERVALS=$REGIONS PER_TARGET_COVERAGE=$pfx/metrics/CollectHsMetrics/$pfx.metrics.per.target.tsv

picard -Xmx4g CollectHsMetrics I=$bam O=$pfx/metrics/CollectHsMetrics/$pfx.REGIONS.metrics R=$ref BAIT_INTERVALS=$REGIONS TARGET_INTERVALS=$REGIONS PER_TARGET_COVERAGE=$pfx/metrics/CollectHsMetrics/$pfx.REGIONS.metrics.per.target.tsv


#picard CollectInsertSizeMetrics
picard -Xms4096m -Xmx30000m CollectInsertSizeMetrics I=$bam O=$pfx/metrics/CollectInsertSizeMetrics/$pfx.insert_size_metrics.txt H=$pfx/metrics/CollectInsertSizeMetrics/$pfx.insert_size_histogram.pdf M=0.5

#CalculateHsMetrics
java -jar $hsMetrics BI=$TARGET TI=$TARGET I=$bam O=$pfx/metrics/CalculateHsMetrics/$pfx.hsmetrics.tsv R=$ref PER_TARGET_COVERAGE=$pfx/metrics/CalculateHsMetrics/$pfx.pertarget_coverage.tsv VALIDATION_STRINGENCY=SILENT

#samtoolsstats
samtools stats $bam -r $ref > $pfx/metrics/samtoolsstats/$pfx.samtoolsstats.txt
########################################################

######## Variants ########
### the -q 1 and -B are the setting that the GDC uses, so we decided to add them
samtools mpileup --ignore-RG --ff 0xF00 -f $ref -q 1 -B $pfx/Alignment/$pfx.recalibrated.final.bam -o $pfx/Variants/$pfx.recalibrated.final.bam.mpileup &&

varscan mpileup2snp $pfx/Variants/$pfx.recalibrated.final.bam.mpileup --mincoverage 8 --min-reads2 2 --min-avg-qual 15 --min-var-freq 0.01 --min-freq-for-hom 0.75 --p-value 0.99 --output-vcf 1 > $pfx/Variants/$pfx.varscan.snp.vcf &&

varscan mpileup2indel $pfx/Variants/$pfx.recalibrated.final.bam.mpileup --mincoverage 8 --min-reads2 2 --min-avg-qual 15 --min-var-freq 0.01 --min-freq-for-hom 0.75 --p-value 0.99 --output-vcf 1 > $pfx/Variants/$pfx.varscan.indel.vcf &&

######## annotation ########

$annovar/table_annovar.pl $pfx/Variants/$pfx.varscan.snp.vcf $annovar/humandb/ -buildver hg38 -out $pfx/Annotation/$pfx.varscan.snp.vcf -remove -protocol refGene,cytoBand,genomicSuperDups,exac03,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_eas,1000g2015aug_eur,avsnp150,avsnp150_NonFlagged,dbnsfp33a,ljb26_all,cosmic83_coding,cosmic83_noncoding,cosmic83_coding.hem.lym,clinvar_20190305,dbnsfp31a_interpro,esp6500siv2_all,gnomad211_exome,gnomad211_genome,gatk4_mutect2_4136_pon.with.avsnp150_NonFlagged -operation g,r,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f -nastring . -vcfinput &&

$annovar/table_annovar.pl $pfx/Variants/$pfx.varscan.indel.vcf $annovar/humandb/ -buildver hg38 -out $pfx/Annotation/$pfx.varscan.indel.vcf -remove -protocol refGene,cytoBand,genomicSuperDups,exac03,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_eas,1000g2015aug_eur,avsnp150,avsnp150_NonFlagged,dbnsfp33a,ljb26_all,cosmic83_coding,cosmic83_noncoding,cosmic83_coding.hem.lym,clinvar_20190305,dbnsfp31a_interpro,esp6500siv2_all,gnomad211_exome,gnomad211_genome,gatk4_mutect2_4136_pon.with.avsnp150_NonFlagged -operation g,r,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f -nastring . -vcfinput &&

####### combine snp and indel #######

tail -n+2 $pfx/Annotation/$pfx.varscan.indel.vcf.hg38_multianno.txt|cat $pfx/Annotation/$pfx.varscan.snp.vcf.hg38_multianno.txt - > $pfx/Annotation/$pfx.varscan.snp.indel.vcf.hg38_multianno.txt &&

################uGT######################

gatk="/lustre/work/javeediqbal/abouska/Sequencing_Pipeline/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar"
bam="$pfx/Alignment/$pfx.recalibrated.final.bam"

java -Xmx15g -jar $gatk -T UnifiedGenotyper -R $ref -I $bam -o $pfx/Variants/$pfx.uGT.hg38.vcf --genotype_likelihoods_model BOTH


######## actual used memory ########
cat /cgroup/memory/slurm/uid_${UID}/job_${SLURM_JOBID}/memory.max_usage_in_bytes
