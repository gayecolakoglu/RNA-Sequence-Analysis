# RNA-Sequence-Analysis

## Abstract
In this project, I am going to try to characterize natural mutations found in HCoV-229E variants. The
data that have been used in that project can be obtained from [here](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA656274&o=acc_s%3Aa). This data provided by monitoring the cells that were infected with 0.01 multiplicity of infection (MOI) of HCoV-229E virus.

Then total RNA was extracted at 0, 6, 24, 48, and 72 post-HCoV-229E infection. At each time point,
three biological repeats were tested.

Samples in the dataset were downloaded in fastq format. These fastq format files obtained were summarized and visualized with the FastQC tool. After summarizing and visualizing all these fastq format files with the FastQC tool, Multiqc has been applied to summarize the fastqc results.

Six samples have been chosen from fifteen samples of High throughput sequencing (3 normal, 3 diseases). Considering that this will be the time when the disease is most advanced, the cell sample that was kept for 72 hours was selected.

We can obtain detailed information about the project and data with the following accession numbers and URL links.
- PMC: [PMC7906355](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7906355/)
- BioProject: [PRJNA656274](https://www.ncbi.nlm.nih.gov//bioproject/PRJNA656274) (same PRJN for normal/controls)
- GEO: [GSE155986](https://www.ncbi.nlm.nih.gov//geo/query/acc.cgi?acc=GSE155986)
- SRA: [SRP276927](https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP276927)

## Pre-processing
In the coding part, I am going to follow the repository of Griffith lab for [Environment](https://github.com/griffithlab/rnaseq_tutorial/wiki/Environment) and [Installation](https://github.com/griffithlab/rnaseq_tutorial/wiki/Installation).

Getting fastq files from .sra files which are inside SRR files
```
	fastq-dump SRR12424243 --split-files
	fastq-dump SRR12424244 --split-files
	fastq-dump SRR12424245 --split-files
	fastq-dump SRR12424255 --split-files
	fastq-dump SRR12424256 --split-files
	fastq-dump SRR12424257 --split-files
```

We got 1 fastq file for each SRR because our samples are single strand if else we will get 2 fastq file for each SRR

To examine fastq files get FastQC report for each sample
```
	fastqc *.fastq
```

Generate a single summary report across all samples which is MultiQC report
```
	 multiqc . 
```

Since we do not have any adapter contamination we do not have to perform felxbar trim but here trim code(for two strands) just to know
```
	$RNA_HOME/student_tools/flexbar-3.4.0-linux/flexbar --adapter-min-overlap 7
	--adapter-trim-end RIGHT --adapters $RNA_REFS_DIR/illumina_multiplex.fa
	--pre-trim-left 13 --max-uncalled 300 --min-read-length 25 --threads 8 --zip-output
	GZ --reads $RNA_DATA_DIR2/SRR12424243/SRR12424243_1.fastq --reads2
	$RNA_DATA_DIR2/SRR12424243/SRR12424243_2.fastq --target $RNA_DATA_TRIM/SRR12424243
 ```
 
 ## RNA-seq Pipeline
 ### Alignment
 In this section we map the reads in our FASTQ files to a reference genome. The output of this step will be a SAM/BAM file for each data set.

Since we have single strand we are going to use -u option(not -1 and -2):
```
hisat2 -p 16 --rg-id=SRR12424243 --rg SM:NORMAL --rg LB:SRR12424243-NORMAL --summary-file SRR12424243.txt --rg PL:ILLUMINA -x $RNA_REF_INDEX --dta --rna-strandness RF -u $RNA_DATA_DIR2/SRR12424243.fastq.gz -S ./SRR12424243.sam

hisat2 -p 16 --rg-id=SRR12424244 --rg SM:NORMAL --rg LB:SRR12424244-NORMAL --summary-file SRR12424244.txt --rg PL:ILLUMINA -x $RNA_REF_INDEX --dta --rna-strandness RF -u $RNA_DATA_DIR2/SRR12424244.fastq.gz -S ./SRR12424244.sam

hisat2 -p 16 --rg-id=SRR12424245 --rg SM:NORMAL --rg LB:SRR12424245-NORMAL --summary-file SRR12424245.txt --rg PL:ILLUMINA -x $RNA_REF_INDEX --dta --rna-strandness RF -u $RNA_DATA_DIR2/SRR12424245.fastq.gz -S ./SRR12424245.sam

hisat2 -p 16 --rg-id=SRR12424255 --rg SM:INFECTED --rg LB:SRR12424255-INFECTED --summary-file SRR12424255.txt --rg PL:ILLUMINA -x $RNA_REF_INDEX --dta --rna-strandness RF -u $RNA_DATA_DIR2/SRR12424255.fastq.gz -S ./SRR12424255.sam

hisat2 -p 16 --rg-id=SRR12424256 --rg SM:INFECTED --rg LB:SRR12424256-INFECTED --summary-file SRR12424256.txt --rg PL:ILLUMINA -x $RNA_REF_INDEX --dta --rna-strandness RF -u $RNA_DATA_DIR2/SRR12424256.fastq.gz -S ./SRR12424256.sam

hisat2 -p 16 --rg-id=SRR12424257 --rg SM:INFECTED --rg LB:SRR12424257-INFECTED --summary-file SRR12424257.txt --rg PL:ILLUMINA -x $RNA_REF_INDEX --dta --rna-strandness RF -u $RNA_DATA_DIR2/SRR12424257.fastq.gz -S ./SRR12424257.sam
```

convert .sam files to .bam files and sort by aligned position 
```
	samtools sort -@ 8 -o SRR12424243.bam SRR12424243.sam
	samtools sort -@ 8 -o SRR12424244.bam SRR12424244.sam 
	samtools sort -@ 8 -o SRR12424245.bam SRR12424245.sam 
	samtools sort -@ 8 -o SRR12424255.bam SRR12424255.sam
	samtools sort -@ 8 -o SRR12424256.bam SRR12424256.sam
	samtools sort -@ 8 -o SRR12424257.bam SRR12424257.sam
```

Perform FastQC and MultiQC to .bam files
```
fastqc *.bam
multiqc .
```

Make a single BAM file combining all NORMAL data and another for all INFECTED data.
```
cd $RNA_HOME/alignments2/hisat2
java -Xmx2g -jar $RNA_HOME/tools/picard.jar MergeSamFiles OUTPUT=NORMAL.bam
INPUT=SRR12424243.bam INPUT=SRR12424244.bam INPUT=SRR12424245.bam
java -Xmx2g -jar $RNA_HOME/tools/picard.jar MergeSamFiles OUTPUT=INFECTED.bam
INPUT=SRR12424255.bam INPUT=SRR12424256.bam INPUT=SRR12424257.bam
```

For the purpose of viewing our alignments we need to index our BAM files. We will use samtools
index for this purpose. For convenience later, index all bam files.
```
echo $RNA_ALIGN_DIR2
cd $RNA_ALIGN_DIR2
find *.bam -exec echo samtools index {} \; | sh
```
## Expression
### StringTie
I used Stringtie to generate expression estimates from the SAM/BAM files generated by HISAT2.

StringTie takes as input spliced alignments in coordinate-sorted SAM/BAM/CRAM format and produces a GTF output which consists of assembled transcript structures and their estimated expression levels (FPKM/TPM and base coverage values). Firstly, we created directory as usual:
```
cd $RNA_HOME/
mkdir -p expression/stringtie/ref_only/
cd expression/stringtie/ref_only/
```

We’re going to use the ENSEMBL indexed hg38 human genome and human genome gtf files. With the following code, we’re going to produce the assembled transcripts GTF and gene abundance estimates for each sample.
```
stringtie -p 16 -G $RNA_REF_GTF -e -B -o SRR12424243/transcripts.gtf -A SRR12424243/gene_abundances.tsv $RNA_ALIGN_DIR2/SRR12424243.bam
stringtie -p 16 -G $RNA_REF_GTF -e -B -o SRR12424244/transcripts.gtf -A SRR12424244/gene_abundances.tsv $RNA_ALIGN_DIR2/SRR12424244.bam
stringtie -p 16 -G $RNA_REF_GTF -e -B -o SRR12424245/transcripts.gtf -A SRR12424245/gene_abundances.tsv $RNA_ALIGN_DIR2/SRR12424245.bam
stringtie -p 16 -G $RNA_REF_GTF -e -B -o SRR12424255/transcripts.gtf -A SRR12424255/gene_abundances.tsv $RNA_ALIGN_DIR2/SRR12424255.bam
stringtie -p 16 -G $RNA_REF_GTF -e -B -o SRR12424256/transcripts.gtf -A SRR12424256/gene_abundances.tsv $RNA_ALIGN_DIR2/SRR12424256.bam
stringtie -p 16 -G $RNA_REF_GTF -e -B -o SRR12424257/transcripts.gtf -A SRR12424257/gene_abundances.tsv $RNA_ALIGN_DIR2/SRR12424257.bam
```

What does the raw output from Stringtie look like? We can look raw output of all samples with below code:
```
less -S SRR12424243/transcripts.gtf
```

But it’s output going to be little complicated so instead we can run below code to limit the view to transcript records and their expression values (FPKM and TPM values)
```
awk '{if ($3=="transcript") print}' UHR_Rep1/transcripts.gtf | cut -f 1,4,9 | less
```

Create a tidy expression matrix files for the StringTie results both gene and transcript level and take into account FPKM, TPM
```
wget https://raw.githubusercontent.com/griffithlab/rnaseq_tutorial/master/scripts/stringtie_expression_matrix.pl
chmod +x stringtie_expression_matrix.pl

./stringtie_expression_matrix.pl --expression_metric=TPM --result_dirs='SRR12424243,SRR12424244,SRR12424245,SRR12424255,SRR12424256,SRR12424257' --transcript_matrix_file=transcript_tpm_all_samples.tsv --gene_matrix_file=gene_tpm_all_samples.tsv

./stringtie_expression_matrix.pl --expression_metric=FPKM --result_dirs='SRR12424243,SRR12424244,SRR12424245,SRR12424255,SRR12424256,SRR12424257' --transcript_matrix_file=transcript_fpkm_all_samples.tsv --gene_matrix_file=gene_fpkm_all_samples.tsv

./stringtie_expression_matrix.pl --expression_metric=Coverage --result_dirs='SRR12424243,SRR12424244,SRR12424245,SRR12424255,SRR12424256,SRR12424257' --transcript_matrix_file=transcript_coverage_all_samples.tsv --gene_matrix_file=gene_coverage_all_samples.tsv
```

Let’s look at inside the .tsv files with below codes:
```
head transcript_tpm_all_samples.tsv gene_tpm_all_samples.tsv
head transcript_fpkm_all_samples.tsv gene_fpkm_all_samples.tsv
```

### Htseq-Count
HTSeq-count counts the number of the reads from each bam file that map to the genomic features in the provided annotation file. For each feature (a gene for example) we will obtain a numerical value associated with the expression of that feature in our sample (i.e. the number of reads that were aligned to that gene).

Given a file with aligned sequencing reads and a list of genomic features, a common task is to count
how many reads map to each feature.

Run htseq-count on alignments instead to produce raw counts instead of FPKM/TPM values for
differential expression analysis and calculate gene-level counts with below code for each sample:
```
htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id /media/data03/cihan/workspace/rnaseq/alignments2/hisat2/SRR12424243.bam $RNA_REF_GTF > SRR12424243.tsv
htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id /media/data03/cihan/workspace/rnaseq/alignments2/hisat2/SRR12424244.bam $RNA_REF_GTF > SRR12424244.tsv
htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id /media/data03/cihan/workspace/rnaseq/alignments2/hisat2/SRR12424245.bam $RNA_REF_GTF > SRR12424245.tsv

htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id /media/data03/cihan/workspace/rnaseq/alignments2/hisat2/SRR12424255.bam $RNA_REF_GTF > SRR12424255.tsv
htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id /media/data03/cihan/workspace/rnaseq/alignments2/hisat2/SRR12424256.bam $RNA_REF_GTF > SRR12424256.tsv
htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id /media/data03/cihan/workspace/rnaseq/alignments2/hisat2/SRR12424257.bam $RNA_REF_GTF > SRR12424257.tsv
 ```
  
Merge results files into a single matrix for use in edgeR:
```
cd $RNA_HOME/expression2/htseq_counts/
join SRR12424243.tsv SRR12424244.tsv | join - SRR12424245.tsv | join - SRR12424255.tsv | join - SRR12424256.tsv | join -SRR12424257.tsv > gene_read_counts_table_all.tsv
echo "GeneID SRR12424243 SRR12424244 SRR12424245 SRR12424255 SRR12424256 SRR12424257" > header.txt
cat header.txt gene_read_counts_table_all.tsv | grep -v "__" | perl -ne 'chomp $_; $_ =~ s/\s+/\t/g; print "$_\n"' > gene_read_counts_table_all_final.tsv 
rm -f gene_read_counts_table_all.tsv header.txt
head gene_read_counts_table_all_final.tsv
```

View the .tsv file to look at the calculated gene-level counts of each in each sample:
```
head gene_read_counts_table_all_final.tsv
```

## Differential Expression
### Ballgown DE Analysis
Ballgown to compare the NORMAL and INFECTED conditions.
We begun with creating necassary working directory:
```
mkdir -p $RNA_HOME/de/ballgown/ref_only/
cd $RNA_HOME/de/ballgown/ref_only/
```

We are going to perform NORMAL vs. INFECTED comparison, using all replicates, for known
(reference only mode) transcripts:

First we created a file that lists our 6 expression files.
```
printf
"\"ids\",\"type\",\"path\"\n\"SRR12424243\",\"NORMAL\",\"$RNA_HOME/expression/stringtie/re
f_only/SRR12424243\"\n\"SRR12424244\",\"NORMAL\",\"$RNA_HOME/expression/stringtie/r
ef_only/SRR12424244\"\n\"SRR12424245\",\"NORMAL\",\"$RNA_HOME/expression/stringtie/
ref_only/SRR12424245\"\n\"SRR12424255\",\"INFECTED\",\"$RNA_HOME/expression/stringt
ie/ref_only/SRR12424255\"\n\"SRR12424256\",\"INFECTED\",\"$RNA_HOME/expression/stri
ngtie/ref_only/SRR12424256\"\n\"SRR12424257\",\"INFECTED\",\"$RNA_HOME/expression/
stringtie/ref_only/SRR12424257\"\n" > NORMAL_vs_INFECTED.csv
```

Let’s view that file:
```
cat NORMAL_vs_INFECTED.csv
```
Then start an R session where we will examine the NORMAL_vs_INFECTED.csv. You can reach the code by [Tutorial_Part1_ballgown.R]()

