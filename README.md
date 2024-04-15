# Next Generation Sequencing (NGS) data processing, quality control, annotation and variant calling
This script genetaes .vcf files from fastq files. This comprehensive pipeline is used for Whole Exome Sequencing (WES) analysis, particularly mapping reads to the human reference genome (hg38), identifying variants, and annotating these variants. Here’s a step-by-step explanation of each segment of the script and why these steps are critical in a genomic sequencing workflow:

1. Environment Setup and Variable Initialization
The script starts by setting up a job on a high-performance computing cluster using SLURM directives. These directives specify the job resources such as the number of nodes, memory, and execution time. Variables are initialized with paths to reference files and software modules like bwa, samtools, picard, and gatk are loaded. These tools are essential for sequence alignment, data manipulation, and variant calling.

2. Alignment
BWA-MEM for Alignment: Reads (r1, r2, r3, r4) are aligned to the reference genome using bwa mem, which is suitable for 70bp-1Mbp reads. This step is crucial as it maps DNA sequences to the human genome, allowing for the identification of the genomic coordinates of each read.
SAM to BAM Conversion: The alignment output in SAM format is converted to BAM for more efficient storage and handling.
Sorting and Indexing: BAM files are sorted and indexed using samtools. Sorting is necessary for many downstream processes like duplicate marking and variant calling.
3. Post-Alignment Processing
Read Group Addition: Picard is used to add read group information, which is important for distinguishing data among samples and analysis runs.
Merging BAM Files: Alignments from different lanes are merged to handle them as a single dataset.
Marking Duplicates: Duplicate reads are marked using Picard MarkDuplicates. This step is crucial to ensure that duplicates arising from sequencing or library preparation do not bias variant calling.
Base Quality Score Recalibration (BQSR): Performed using GATK. This step adjusts the base quality scores based on empirical data, improving the accuracy of variant calls.
4. Variant Calling and Recalibration
Variant Calling with GATK’s UnifiedGenotyper or MuTect2: Variants are identified from the sequence data. These tools use statistical models to call SNPs and indels accurately.
Variant Recalibration: Further recalibrates the variants to adjust the quality scores based on various covariates and known variant sites.
5. Annotation
Using ANNOVAR, variants are annotated with information from various databases. This annotation is critical for understanding the biological impact of each variant, such as whether they are known disease-causing mutations or are located in important genomic regions like exons, introns, or regulatory elements.

6. Metrics and Quality Control
Throughout the pipeline, various metrics and quality control steps are implemented:

HsMetrics: Provides metrics specific to hybrid selection and target enrichment efficiency.
Insert Size Metrics: Helps in identifying library preparation artifacts.
Samtools stats: Offers comprehensive statistics about the sequencing run, helping to diagnose problems in sequencing and alignment.
7. Additional Processes
Variant Filtering and Combination: Filters and combines SNP and indel data into a single file for easier downstream analysis.
Memory Usage Monitoring: Tracks the peak memory usage of the job, which is vital for optimizing the computational resources for future runs.
Importance of the Steps
Each step in the pipeline is designed to refine and improve the quality of the final data. Accurate alignment, careful handling of duplicates, quality recalibration, and robust variant calling are essential to ensure that the final variant calls are reliable and accurate. These steps are necessary to derive meaningful biological insights, which can be critical for research and clinical applications such as disease diagnosis, personalized medicine, and understanding genetic diversity.
