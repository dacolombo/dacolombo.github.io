---
title: "Genomics Project"
date: 2021-04-30T10:25:32+02:00
draft: false
cover: "/genomics/Reference_annotation.png"
description: "In this project me and two of my collegues performed a de novo assembly of a chloroplast genome. We run the assembly on a server through the CLI, and then compared the results we obtained using short and long reads." 
keywords: ["Bioinformatics","Genomics","Assembly","Mapping","Chloroplast"]
tags: ["Bioinformatics","Genomics","Assembly","Mapping","Chloroplast"]
----


The main topics of this project are:
- [Data Mining](https://www.danielecolombo.me/projects/genomics/#data-mining)
- [Reads trimming](https://www.danielecolombo.me/projects/genomics/#reads-trimming)
- [Quality control](https://www.danielecolombo.me/projects/genomics/#quality-control)
- [Alignment](https://www.danielecolombo.me/projects/genomics/#alignment)
	- [bowtie2](https://www.danielecolombo.me/projects/genomics/#alignment-of-illumina-reads-bowtie2)
	- [minimap2](https://www.danielecolombo.me/projects/genomics/#alignment-of-nanopore-reads-minimap2)
- [Assembly](https://www.danielecolombo.me/projects/genomics/#the-assembly)
	- [ABySS](https://www.danielecolombo.me/projects/genomics/#abyss)
	- [Canu](https://www.danielecolombo.me/projects/genomics/#canu)
- [Annotation](https://www.danielecolombo.me/projects/genomics/#annotation)

# Evaluation of the differences in chloroplast genome assembly using Illumina and Oxford Nanopore reads
This is the project that me and two of my collegues had done for the evaluation of the Genomics course. We could access for the first time to a server that the University provided to us in order to do all the computational tasks.\
The project consisted in performing a chloroplast genome assembly and comparing the differences between results obtained with short and long reads.


# Data Mining
The first thing that we had to do was of course finding the data to work with. In order to do that, we searched in the taxonomy database of [NCBI](https://www.ncbi.nlm.nih.gov/).\
We also used the [pfaf](https://pfaf.org/user/Default.aspx) website (Plants For A Future) to have some inspiration in our search: we looked for plants that had edible and medicinal uses.\
After some research, we found a good candidate: [Prunus mandshurica](https://en.wikipedia.org/wiki/Prunus_mandshurica). In fact, the Prunus genus included many other species that had their chloroplast genome mapped to use as a reference down the line. The two SRA experiments that we chose for our analysis were:
- [ERR4762302](https://www.ncbi.nlm.nih.gov/sra/?term=ERR4762302): obtained with Illumina HiSeq 4000 sequencing.
- [ERR4656976](https://www.ncbi.nlm.nih.gov/sra/?term=ERR4656976): obtained with Oxford Nanopore MinION sequencing.

We then downloaded the fastq files of the runs with `fastq-dump`:
{{< code language="bash" >}}
fastq-dump --defline-qual "+" --split-files --gzip --clip ERR4762302
fastq-dump --defline-qual "+" --gzip --clip ERR4656976
{{< /code >}}

Note: after downloading, a folder `ncbi/` will appear in the working directory. This contains the .SRA files that were used by fastq-dump and can be removed in order to save space (space is very precious on a server where a lot of people will be working).


# Reads Trimming
When working with Illumina and Oxford Nanopore reads, the first thing to do is to trim them in order to remove the adapters (if present) used during the sequencing. 


#### Illumina Reads Trimming
To do this for Illumina reads, a FASTA file containing all of the possible adapters is needed. The command to perform the trimming is [`fastq-mcf`](https://github.com/ExpressionAnalysis/ea-utils/blob/wiki/FastqMcf.md), and was run with the following options:
{{< code language="bash" >}}
fastq-mcf -o trimmed.ERR4762302_1.fastq.gz -o trimmed.ERR4762302_2.fastq.gz\
«dir to adapters»/IlluminaAdapters.fasta\
ERR4762302_1.fastq.gz ERR4762302_2.fastq.gz
{{< /code >}}
In our case, the output of this command indicated that no adapters were found: this is due to the fact that the reads were already processed before being uploaded to the SRA database. This was also confirmed by the fact that the reads had a length distribution of mean 146 and variance around 15: this indicates that the reads are already trimmed, as Illumina reads are usually of 150 bp strictly.


#### Nanopore Reads Trimming
For Nanopore reads, the trimming can be performed with the command [`porechop`](https://github.com/rrwick/Porechop). This tool searches for the adapters at the ends and in the middle of a sample of the reads: the adapters with the higher mapping percentage are then searched in every read and removed. We run the following command:
{{< code language="bash" >}}
porechop -i ERR4656976.fastq.gz -o ERR4656976_trimmed.fastq.gz --threads 10
{{< /code >}}

The two original fastq files can then be safely removed as they will no longer be used.


# Quality Control
After downloading and trimming the reads, the first thing was to do some quality checks on the reads that we chose. To do this, `fastqc` was used. The output of this command is a .html file, that can be visualized in a browser. In particular, to run this command the only argument to pass is the fastq file.\
For a guide on how to read the output of fastqc, refer to the [fastqc documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/)


# Alignment
The next step of the project was to select the reads that are part of the chloroplast of P. mandshurica. To do that, a reference chloroplast was used. Chloroplasts genomes are usually highly conserved between species of the same genus. That's why we chose as reference some assembled chloroplast genomes of plants belonging to the Prunus genus, in particular:
- P. persica
- P. japonica
- P. mume
- P. avium
- P. dulcis

To use them in the server we had access to, we had to download the fasta files of these genomes locally, and then upload them through the command line with `scp`.\
We then performed the alignment of both Illumina and Nanopore reads (the how will be explained in the following section) for all the refereces, in order to select the reference for which the reads showed the highest coverage.\
The reads mapped to the references more or less with the same percentage (between 10% and 12%), also resulting in similar coverages. We chose P. persica as our reference because it showed slightly higher numbers.


#### Alignment of Illumina reads: bowtie2
One of the main tools to align Illumina reads to a reference is `bowtie2`. For an overview on how to use it, refer to its [github page](https://github.com/BenLangmead/bowtie2).\
The first thing to do in order to run bowtie2 is creating an index of the reference, which will be stored as a .fai file. This can be done with the `bowtie2-build` command:
{{< code language="bash" >}}
bowtie2-build Prunus_persica_chloroplast.fasta persica_CHL
{{< /code >}}

After creating the reference index, the alignment can be performed. As the output of bowtie2 is a .SAM file, which occupies a lot of space, it's suggested to directly convert it into the BAM format. This can be done with `samtools` (see [samtools manual](https://www.htslib.org/doc/samtools.html)). In particular, `samtools view` was used (by piping the two commands):
{{< code language="bash" >}}
bowtie2 -x «dir to reference»/persica_CHL -p 10\
	-1 cleaned.ERR4762302_1.fq.gz -2 cleaned.ERR4762302_2.fq.gz\
	| samtools view -@ 10 -b -o ERR4762302.bam -
{{< /code >}}

As a last step, the resulting BAM file can be better used (and stored in a more efficient way) by ordering the reads in it by position. This was done with `samtools sort`:
{{< code language="bash" >}}
samtools sort -o ERR4762302.sort.bam ERR4762302.bam
{{< /code >}}


#### Alignment of Nanopore reads: minimap2
The tools we decided to use for the alignment of Nanopore reads was `minimap2`. We referred to its [github readme](https://github.com/lh3/minimap2) and its help page to run it. Also in this case, we piped it with the `samtools view` command in order to convert the SAM output to the BAM format and we then sorted the file by position:
{{< code language="bash" >}}
minimap2 -t 10 -ax map-ont «dir to reference»/Prunus_persica_chloroplast.fasta\
ERR4656976_trimmed.fastq.gz | samtools view -@ 10 -b -o ERR4656976.bam -;

samtools sort -@ 10 -o ERR4656976.sort.bam ERR4656976.bam;
{{< /code >}}


# Reads Extraction
The mapped reads contained in the BAM file can be extracted with `samtools view`:
{{< code language="bash" >}}
#Extract Illumina reads
samtools view -@ 10 -b -F 4 -o ERR4762302_persica_mapped.sort.bam ERR4762302.sort.bam

#Extract Nanopore reads
samtools view -@ 10 -b -F 4 -o ERR4656976_persica_mapped.sort.bam ERR4656976.sort.bam
{{< /code >}}
This command will remove all the reads that are unmapped (mapped with the `4` flag) in the output BAM file.

After extracting the reads, the mapped reads can be counted with:
{{< code language="bash" >}}
samtools view -@ 10 -c «BAM file»
{{< /code >}}

To perform alignment evaluation and proceed with the assembly in further steps, we extracted the reads in the BAM file into a fastq file.\
This procedure was carried out with different tools for Illumina and Nanopore reads:
- For Illumina reads, `samtools fastq` was used:
{{< code language="bash" >}}
samtools fastq -@ 10 -1 ERR4762302_persica_mapped_1.fastq.gz\
	-2 ERR4762302_persica_mapped_2.fastq.gz -c 6\
	-N ERR4762302_persica_mapped.sort.bam
{{< /code >}}
- For Nanopore reads, `bedtools bamtofastq` was used (see [bamtofastq manual](https://bedtools.readthedocs.io/en/latest/content/tools/bamtofastq.html)):
{{< code language="bash">}}
bedtools bamtofastq -i ERR4656976_persica_mapped.sort.bam\
-fq ERR4656976_persica_mapped.fastq;

gzip ERR4656976_persica_mapped.fastq;
{{< /code >}}


# Alignment Evaluation: Coverage
The main parameter to evaluate the alignment of reads on a reference genome is their coverage.\
To compute it, we build the reference index with `samtools faidx` and then we used`bedtools genomecov` on it (see [genomecov manual](https://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html)):
{{< code language="bash" >}}
#Build the reference index
samtools faidx «dir to reference»/Prunus_persica_chloroplast.fasta

#Compute coverage table for Illumina reads
bedtools genomecov -pc -d -ibam ERR4762302_persica_mapped.sort.bam\
-g «dir to reference»/Prunus_persica_chl.fai > ERR4762302_persica_mapped.covbed.txt

#Compute coverage table for Nanopore reads
bedtools genomecov -d -ibam ERR4656976_persica_mapped.sort.bam\
-g «dir to reference»/Prunus_persica_chl.fai > ERR4656976_persica_mapped.covbed.txt
{{< /code >}}

The `.covbed.txt` file generated contains three columns:
1. The name of the reference sequence
2. The position of the base in the reference sequence
3. The coverage for that base

By importing this file in R, it's possible to use the data in the third column to plot the per-base-coverage and to compute its mean:
{{< code language="R" >}}
# Coverage of mapped Illumina reads on P. persica
Covbed_Illumina <- read.delim("ERR4762302_persica_mapped.covbed.txt", header=FALSE)
plot(Covbed_Illumina$V3, type="l", ylab="Coverage", main="Coverage of Illumina reads")
mean(Covbed_Illumina$V3)

# Coverage of mapped Nanopore reads on P. persica
Covbed_Nanopore <- read.delim("ERR4656976_persica_mapped.covbed.txt", header=FALSE)
plot(Covbed_Nanopore$V3, type="l", ylab="Coverage", main="Coverage of Nanopore reads")
mean(Covbed_Nanopore$V3)
{{< /code >}}

Here are the two resulting plots:

{{< image src="/genomics/Coverage_illumina.jpeg" style="max-width:45%; display:block; float:left; white-space:nowrap; margin-left:15px; margin-right:auto; padding-bottom:15px;" >}}
{{< image position="right" src="/genomics/Coverage_nanopore.jpeg" style="max-width:45%; display:block; white-space:nowrap; overflow:hidden; margin-left:auto; margin-right:15px; padding-bottom:15px">}}

We assumed that the coverages were more than sufficient to satisfy the conditions for a good assembly. We also noted that for the Nanopore reads, there were some high-coverage regions, that we hypothesized were repeated regions.\
We performed a BLAST of the reference to itself, in order to confirm our hypothesis.

{{< image src="/genomics/Reference_dotplot.png" position="center">}}

As can be seen in the resulting dotplot, there are two regions which are inversely repeated. These are exactly positioned in the regions for which the Nanopore reads have an increased coverage. With this we could observe that the reference genome presented Inversed Repeats (IR), which are quite common in chloroplasts regions.\
The fact that the coverage has high differences in specific regions may induce some problems in the assembly, and needs to be taken into account.


# The assembly
The main approaches for sequence assembly are two:
1. Overlap Layout Consensus (OLC): this technique is mainly used for low-coverage long reads.
2. De Bruijn Graph (DBG): this technique is mainly used for high coverage short reads.

To see in details the characteristics of these two techniques, you can read this [article](https://academic.oup.com/bfg/article/11/1/25/191455).

We decided to use these two tools for the assembly:
- [`ABySS`](https://github.com/bcgsc/abyss) to assembly Illumina reads. This tool is based on the DBG algorithm.
- [`canu`](https://canu.readthedocs.io/en/latest/) to assembly Nanopore reads. This tool is based on the OLC algorithm.

#### ABySS
Before running the actual assembly, some steps are required to optimize the parameter of the assembly itself.\
The first thing to do in order to have the best assembly possible is to find the best **kmer size**. This was done following the [ABySS readme](https://github.com/bcgsc/abyss#optimizing-the-parameter-k):
{{< code language="bash">}}
for k in $(seq 17 8 99);
do

mkdir k$k;

abyss-pe -C k$k name=pmandshurica k=$k\
in="../ERR4762302_persica_mapped_1.fastq.gz ../ERR4762302_persica_mapped_2.fastq.gz";

done
{{< /code >}}

After this is done, the following commands were run in order to evaluate the goodness of the assembly for each k value:
{{< code language="bash" >}}
abyss-fac k*/pmandshurica-unitigs.fa
{{< /code >}}
It's possible that the command abyss-fac is not found. In that case, you need to specify the ABySS directory in which the program is installed, that should be:
{{< code language="bash" >}}
/usr/lib/abyss/abyss-fac k*/pmandshurica-unitigs.fa
{{< /code >}}

We obtained the following table as a result:
{{< code language="bash" >}}
n       n:500   L50     min     N80     N50    N20      E-size  max     sum     name
1276    75      23      500     627     1109   1895     1316    4426    72063   k17/pmandshurica-unitigs.fa
132     37      8       540     2764    5892   9827     6655    17362   130085  k25/pmandshurica-unitigs.fa
36      14      4       1600    8245    13415  18115    14100   25583   131013  k33/pmandshurica-unitigs.fa
16      12      4       1397    10908   14789  18142    15250   24549   131714  k41/pmandshurica-unitigs.fa
14      9       2       1236    15651   24528  51139    30162   51139   131577  k49/pmandshurica-unitigs.fa
10      8       3       1870    15659   19689  40082    24657   40082   131557  k57/pmandshurica-unitigs.fa
7       7       2       1870    19177   35228  40082    28336   40082   131553  k65/pmandshurica-unitigs.fa
11      10      3       1610    9716    18121  40082    23307   40082   131690  k73/pmandshurica-unitigs.fa
11      10      3       1610    9770    18121  40082    23305   40082   131760  k81/pmandshurica-unitigs.fa
13      11      3       1492    9778    18296  40082    23165   40082   131865  k89/pmandshurica-unitigs.fa
14      12      3       515     9786    18304  40083    23047   40083   132724  k97/pmandshurica-unitigs.fa
{{< /code >}}
We observed that the best kmer size to use was `65`, as it led to the lowest number of unitigs.


After that, the following thing to do was to find the best sample size of reads to use in order to obtain the best assembly. Given the average coverage of ~1200X of all the reads, we decided to perform some sampling in order to obtain the following estimated coverages:
- 50X: with a sample of 50.000 reads
- 100X: with a sample of 110.000 reads
- 200X: with a sample of 220.000 reads
- 350X: with a sample of 385.000 reads
- 500X: with a sample of 550.000 reads
- 750X: with a sample of 820.000 reads
- 1000X: with a sample of 1.100.000 reads

The amount of reads needed to obtain a wanted coverage was computed through the Lander/Waterman equation (the formula can be viewed in [this document](https://www.illumina.com/documents/products/technotes/technote_coverage_calculation.pdf)).

We performed the sampling with `seqtk sample` (see [seqtk github page](https://github.com/lh3/seqtk)):
{{< code language="bash" >}}
samples=(25000 55000 110000 192500 275000 410000 550000);

for size in "${samples[@]}";
do

seqtk sample -s1000 ERR4762302_persica_mapped_1.fastq.gz ${size} > samples/ERR4762302_$((size/1000))K_1.fastq;
gzip samples/ERR4762302_$((size/1000))K_1.fastq

seqtk sample -s1000 ERR4762302_persica_mapped_2.fastq.gz ${size} > samples/ERR4762302_$((size/1000))K_2.fastq;
gzip samples/ERR4762302_$((size/1000))K_2.fastq;

done
{{< /code >}}

To evaluate which is the best sample size, ABySS was run on all the samples:
{{< code language="bash" >}}
samples=(25000 55000 110000 192500 275000 410000 550000);

for size in "${samples[@]}";
do

size=$((size/1000))

mkdir ${size}K;
cd ${size}K;

abyss-pe name="MandshuricaCHL_AbPE_${size}K" k=65\
in="../samples/ERR4762302_${size}K_1.fastq.gz ../samples/ERR4762302_${size}K_2.fastq.gz";

cd ..;

done;
{{< /code >}}

To check the best assembly, the command `abyss-fac` can be run on the `scaffolds.fa` files as previously shown.\
We observed that the best assembly was performed with the sample of 820.000 reads (750X coverage), as it was assembled in a single scaffold of `131824`bp.


#### Canu
As canu requires lower coverages to perform an assembly, samples of lower sizes were selected by running `seqtk sample`:
{{< code language="bash">}}
samples=(5000 10000 20000 30000);

for size in "${samples[@]}";
do

seqtk sample -s1000 ERR4656976_persica_mapped.fastq.gz ${size} > samples/ERR4656976_$((size/1000))K.fastq;
gzip samples/ERR4656976_$((size/1000))K.fastq

done
{{< /code >}}

Then, we run canu on all the samples with the following comand:
{{< code language="bash">}}
samples=(5000 10000 20000 30000);

for size in "${samples[@]}";
do

size=$((size/1000))

canu -d canu${size}K -p MandshuricaCHL${size} genomeSize=158k\
-nanopore ../samples/ERR4656976_${size}.fastq.gz;

done
{{< /code >}}

We evaluated the best canu assembly with `FastaSeqStats`, a tool that is part of [`GenoToolBox`](https://github.com/aubombarely/GenoToolBox).\
We observed that the best assembly was performed with the sample of 20.000 reads: this resulted in a single contig of `197983`bp.

## Differences in the assembly
The two assemblies resulted in quite different sequences in term of length. We thought that this could be due to the presence of repeated sequences or the presence of overlapped sequences in the ends (this is quite frequent with chloroplast genomes, that are circular).\
The `suggestCircular` flag in the contigs.fasta file resulting from canu assembly was equal to `no`, so we thought to BLAST the sequence to itself in order to do further investigations. This was the resulting dotplot:

{{< image src="/genomics/Canu_dotplot.png" position="center">}}

In this image, the presence of Inversed Reapeats (IR) can be observed. These sequences have probably led to problems in the assembly of Nanopore reads.


## Annotation
The last step was to annotate the assembly, in order to observe the presence of chloroplast genes and evaluate the goodness of the assembly.\
To do this, the tool [`Chlorobox GeSeq`](https://chlorobox.mpimp-golm.mpg.de/geseq.html) was used.

We obtained the following annotations for the two assemblies:

{{< figure position="center" src="/genomics/Abyss_annotation.jpeg"  style="width:75%;" caption="Annotation of ABySS assembly" captionStyle="width:75%;" >}}
{{< figure src="/genomics/Canu_annotation.jpeg" position="center" style="width:75%;" caption="Annotation of Canu assembly" captionStyle="width:75%;" >}}

The assembly performed with ABySS shows a similar annotation to the reference chloroplast genome (that you can observe in the cover image of this post), even though the total length of the assembly is lower than the one of the reference.\
Meanwhile, the assembly performed with canu shows regions repeated two or more times in opposite strands, confirming the presence of multiple IRs in the assembly.


## Conclusions
Even though in general longer reads improve the quality of de novo assembly, in our case it led to problems due to the presence of regions on the reference genome that had an higher local coverage: we identified these regions as Inverted Reapeats (IR). This fact caused Nanopore reads to be assembled in a sequence which is longer and presents multiple IRs. 