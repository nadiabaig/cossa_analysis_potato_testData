## cossa_analysis_potato
Modified pipeline to get the kmers and SNPs for a potato test dataset
```diff
@@ Main steps of the pipeline inludes: @@ 
```

* Kmer tables building using Glistcompare
* Set algebra based comparisons to get genotype specific, shared Kmers
* Getting basic statistics of the kmers list
* Kmer list conversion into fastq file
* Mapping of kmers on thwe reference genome
* Proprocessing and filtering of aligned reads
* Counting number of kmers aligned to different chromosomes
* Generating plots for the aligned kmers to chromosomes
* Variant calling and filtering
* Annotation of filtered variants

 -  (&#x1F534;)  Considering 200X coverage for the genotypes, the whole pipeline might take 4-5 working days on a linux server with 256GB RAM !

   <!-- GETTING STARTED -->
## :desktop_computer: Getting Started

### :keyboard: Installation
1. bedtools
```sh
wget https://github.com/arq5x/bedtools2/releases/download/v2.31.0/bedtools.static
mv bedtools.static bedtools
chmod a+x bedtools 
```
2. TRIMMOMATIC
```sh
apt-get download trimmomatic
mv trimmomatic_0.39+dfsg-2_all.deb ../
mv trimmomatic_0.39+dfsg-2_all.deb ../
dpkg -x trimmomatic_0.39+dfsg-2_all.deb /media/nadia.baig/Kmer
```

3. SRA toolkit
```sh
wget --output-document sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
tar -vxzf sratoolkit.tar.gz
cd sratoolkit.3.0.6-ubuntu64/
export PATH=$PATH:/mnt/d/Softwares/sratoolkit.3.0.6-ubuntu64/bin
which fastq-dump
```
4. bcftools
```sh
sudo apt-get -y install bcftools
```

5. Picard
```sh
 git clone https://github.com/broadinstitute/picard.git
 cd picard/
 ./gradle shadowJar # installing all dependencies
 java -jar build/libs/picard.jar 
```

6. Freebayes
   
```sh
sudo apt-get install freebayes
```

7. IGV (Interative Genome Viewer)
```sh
 sudo apt -y install igv  
```

<!-- USAGE EXAMPLES -->
## Usage
### Setting paths of the required tools and input data
```sh
!/bin/bash
export PATH=$PATH:/mnt/d/Softwares/Glistmaker/bin
path=/mnt/d/Cossa_pipeline_Potatotools
fastq1=/mnt/d/PotatoTools/Laura/Laura_S10_L003_R1_001.fastq.gz
fastq2=/mnt/d/PotatoTools/Laura/Laura_S10_L003_R2_001.fastq.gz
gunzip $fastq1
gunzip $fastq2

wait
fastq11=/mnt/d/PotatoTools/Otolia/Otolia_S3_L001_R1_001.fastq
fastq22=/mnt/d/PotatoTools/Otolia/Otolia_S3_L001_R2_001.fastq

wait
glistmaker $fastq11 -w 31 -o $path/Otolia.all.read1

glistmaker $fastq22 -w 31 -o $path/Otolia.all.read2
echo {"step1 completed"}
#repeat the same with other datasets too
```
### Kmers table generation
```sh
!/bin/bash
glistcompare Otolia.all.read1_31.list Otolia.all.read2_31.list -u -o Otolia.all
wait
glistcompare Otolia.all_31_union.list Otolia.all_31_union.list -i -c 20 -o Otolia.all.cutoff2
wait
glistquery Otolia.all.cutoff2_31_intrsec.list > Otolia.file.kmer
```
### Set alegebra based comparisons
```sh
!/bin/bash
glistcompare  Otolia.all.cutoff2_31_intrsec.list Kuba.all.cutoff2_31_intrsec.list -d -o 1.Otolia-bulk-specific-kmers.cutoff10to20-not-in-Kuba &

glistcompare Otolia.all.cutoff2_31_intrsec.list  1.Otolia-bulk-specific-kmers.cutoff10to20-not-in-Kuba_31_0_diff1.list -d -o 2.Otolia-bulk-specific-kmers.cutoff10to20-in-common-with-Kuba &

glistcompare  Otolia.all.cutoff2_31_intrsec.list Laura.all.cutoff2_31_intrsec.list -d -o 3.Otolia-bulk-specific-kmers.cutoff10to20-not-in-Laura &

glistcompare Otolia.all.cutoff2_31_intrsec.list  3.Otolia-bulk-specific-kmers.cutoff10to20-not-in-Laura_31_0_diff1.list -d -o 4.Otolia-bulk-specific-kmers.cutoff10to20-in-common-with-Laura &

glistcompare  1.Otolia-bulk-specific-kmers.cutoff10to20-not-in-Kuba_31_0_diff1.list 3.Otolia-bulk-specific-kmers.cutoff10to20-not-in-Laura_31_0_diff1.list -d -o 5.Otolia-bulk-specific-kmers.cutoff10to20-from-Laura &

glistcompare 3.Otolia-bulk-specific-kmers.cutoff10to20-not-in-Laura_31_0_diff1.list 1.Otolia-bulk-specific-kmers.cutoff10to20-not-in-Kuba_31_0_diff1.list -d -o 6.Otolia-bulk-specific-kmers.cutoff10to20-from-Kuba &

glistcompare 3.Otolia-bulk-specific-kmers.cutoff10to20-not-in-Laura_31_0_diff1.list 2.Otolia-bulk-specific-kmers.cutoff10to20-in-common-with-Kuba_31_0_diff1.list -d -o 7.Otolia-bulk-specific-kmers.cutoff10to20-from-none &

glistcompare 2.Otolia-bulk-specific-kmers.cutoff10to20-in-common-with-Kuba_31_0_diff1.list 4.Otolia-bulk-specific-kmers.cutoff10to20-in-common-with-Laura_31_0_diff1.list -i -o 8.Otolia-bulk-specific-kmers.cutoff10to20-from-both

```
### Basic stats for unique kmers in each comparison followed by kmer to FASTA conversion
```sh
!/bin/bash
for i in $(ls *diff1.list)
do
glistquery ${i} -stat > ${i}.stat
done
wait
#transforming binary file into text files
for i in $(ls *diff1.list)
do
glistquery ${i} > ${i}.kmer
done

wait

#converting kmers into fastq file
for j in $(ls *.kmer)
do
       python3 kmer_to_fastq.py ${j} > ${j}.fq  #get this script from  (https://github.com/cprodhom/CoSSA-workflows) 
done
wait
```
### Kmers alignment on a reference genome, followed by proprocessing and read filtering

```sh
!/bin/bash
reference="/mnt/d/PotatoTools/Agria-Ref/Agria_21082020/Potato/dAg_v1.0/Agria_assembly_final_2020_21_08.fasta"
bwa index $reference
wait

for k in $(ls *list.kmer.fq)
do
       bwa aln $reference ${k} > ${k}.sai
done

wait
input_dir="/mnt/d/Cossa_pipeline_Potatotools"
l1="1.Otolia-bulk-specific-kmers.cutoff10to20-not-in-Kuba_31_0_diff1.list.kmer.fq.sai"
l5="5.Otolia-bulk-specific-kmers.cutoff10to20-from-Laura_31_0_diff1.list.kmer.fq.sai"
l2="2.Otolia-bulk-specific-kmers.cutoff10to20-in-common-with-Kuba_31_0_diff1.list.kmer.fq.sai"
l6="6.Otolia-bulk-specific-kmers.cutoff10to20-from-Kuba_31_0_diff1.list.kmer.fq.sai"
l3="3.Otolia-bulk-specific-kmers.cutoff10to20-not-in-Laura_31_0_diff1.list.kmer.fq.sai"
l7="7.Otolia-bulk-specific-kmers.cutoff10to20-from-none_31_0_diff1.list.kmer.fq.sai"
l4="4.Otolia-bulk-specific-kmers.cutoff10to20-in-common-with-Laura_31_0_diff1.list.kmer.fq.sai"

#not using for loop as it takes time other alternative is to run multiple commands in parallel with & operator

bwa samse ${reference} $l1 1.Otolia-bulk-specific-kmers.cutoff10to20-not-in-Kuba_31_0_diff1.list.kmer.fq > 1.Otolia-bulk-specific-kmers.cutoff10to20-not-in-Kuba_31_0_diff1.list.kmer.fq.sai.tmp &

bwa samse ${reference} $l2 2.Otolia-bulk-specific-kmers.cutoff10to20-in-common-with-Kuba_31_0_diff1.list.kmer.fq > 2.Otolia-bulk-specific-kmers.cutoff10to20-in-common-with-Kuba_31_0_diff1.list.kmer.fq.sai.tmp &

bwa samse  ${reference} $l3 3.Otolia-bulk-specific-kmers.cutoff10to20-not-in-Laura_31_0_diff1.list.kmer.fq > 3.Otolia-bulk-specific-kmers.cutoff10to20-not-in-Laura_31_0_diff1.list.kmer.fq.sai.tmp &

bwa samse  ${reference} $l4 4.Otolia-bulk-specific-kmers.cutoff10to20-in-common-with-Laura_31_0_diff1.list.kmer.fq  > 4.Otolia-bulk-specific-kmers.cutoff10to20-in-common-with-Laura_31_0_diff1.list.kmer.fq.sai.tmp &

bwa samse  ${reference} $l5 5.Otolia-bulk-specific-kmers.cutoff10to20-from-Laura_31_0_diff1.list.kmer.fq  > 5.Otolia-bulk-specific-kmers.cutoff10to20-from-Laura_31_0_diff1.list.kmer.fq.sai.tmp &

bwa samse  ${reference} $l6 6.Otolia-bulk-specific-kmers.cutoff10to20-from-Kuba_31_0_diff1.list.kmer.fq > 6.Otolia-bulk-specific-kmers.cutoff10to20-from-Kuba_31_0_diff1.list.kmer.fq.sai.tmp &

bwa samse  ${reference} $l7 7.Otolia-bulk-specific-kmers.cutoff10to20-from-none_31_0_diff1.list.kmer.fq > 7.Otolia-bulk-specific-kmers.cutoff10to20-from-none_31_0_diff1.list.kmer.fq.sai.tmp &

wait

samtools view -Sbh 1.Otolia-bulk-specific-kmers.cutoff10to20-not-in-Kuba_31_0_diff1.list.kmer.fq.sai.tmp -o  1.Otolia-bulk-specific-kmers.cutoff10to20-not-in-Kuba_31_0_diff1.list.kmer.fq.sai.tmp.bam &

samtools view -Sbh 2.Otolia-bulk-specific-kmers.cutoff10to20-in-common-with-Kuba_31_0_diff1.list.kmer.fq.sai.tmp -o 2.Otolia-bulk-specific-kmers.cutoff10to20-in-common-with-Kuba_31_0_diff1.list.kmer.fq.sai.tmp.bam &

samtools view -Sbh 3.Otolia-bulk-specific-kmers.cutoff10to20-not-in-Laura_31_0_diff1.list.kmer.fq.sai.tmp -o 3.Otolia-bulk-specific-kmers.cutoff10to20-not-in-Laura_31_0_diff1.list.kmer.fq.sai.tmp.bam &

samtools view -Sbh 4.Otolia-bulk-specific-kmers.cutoff10to20-in-common-with-Laura_31_0_diff1.list.kmer.fq.sai.tmp -o 4.Otolia-bulk-specific-kmers.cutoff10to20-in-common-with-Laura_31_0_diff1.list.kmer.fq.sai.tmp.bam &

samtools view -Sbh 5.Otolia-bulk-specific-kmers.cutoff10to20-from-Laura_31_0_diff1.list.kmer.fq.sai.tmp -o 5.Otolia-bulk-specific-kmers.cutoff10to20-from-Laura_31_0_diff1.list.kmer.fq.sai.tmp.bam &

samtools view -Sbh 6.Otolia-bulk-specific-kmers.cutoff10to20-from-Kuba_31_0_diff1.list.kmer.fq.sai.tmp -o 6.Otolia-bulk-specific-kmers.cutoff10to20-from-Kuba_31_0_diff1.list.kmer.fq.sai.tmp.bam &

samtools view -Sbh 7.Otolia-bulk-specific-kmers.cutoff10to20-from-none_31_0_diff1.list.kmer.fq.sai.tmp -o 7.Otolia-bulk-specific-kmers.cutoff10to20-from-none_31_0_diff1.list.kmer.fq.sai.tmp.bam &

wait

samtools sort -T 1.Otolia-bulk-specific-kmers.cutoff10to20-not-in-Kuba_31_0_diff1.list.kmer.fq.sai.sort.tmp 1.Otolia-bulk-specific-kmers.cutoff10to20-not-in-Kuba_31_0_diff1.list.kmer.fq.sai.tmp.bam -o 1.Otolia-bulk-specific-kmers.cutoff10to20-not-in-Kuba_31_0_diff1.list.kmer.fq.sai.tmp.srt.bam &

samtools sort -T 2.Otolia-bulk-specific-kmers.cutoff10to20-in-common-with-Kuba_31_0_diff1.list.kmer.fq.sai.sort.tmp 2.Otolia-bulk-specific-kmers.cutoff10to20-in-common-with-Kuba_31_0_diff1.list.kmer.fq.sai.tmp.bam -o 2.Otolia-bulk-specific-kmers.cutoff10to20-in-common-with-Kuba_31_0_diff1.list.kmer.fq.sai.tmp.srt.bam &

samtools sort -T 3.Otolia-bulk-specific-kmers.cutoff10to20-not-in-Laura_31_0_diff1.list.kmer.fq.sai.sort.tmp 3.Otolia-bulk-specific-kmers.cutoff10to20-not-in-Laura_31_0_diff1.list.kmer.fq.sai.tmp.bam -o 3.Otolia-bulk-specific-kmers.cutoff10to20-not-in-Laura_31_0_diff1.list.kmer.fq.sai.tmp.srt.bam &

samtools sort -T 4.Otolia-bulk-specific-kmers.cutoff10to20-in-common-with-Laura_31_0_diff1.list.kmer.fq.sai.sort.tmp 4.Otolia-bulk-specific-kmers.cutoff10to20-in-common-with-Laura_31_0_diff1.list.kmer.fq.sai.tmp.bam -o 4.Otolia-bulk-specific-kmers.cutoff10to20-in-common-with-Laura_31_0_diff1.list.kmer.fq.sai.tmp.srt.bam

samtools sort -T 5.Otolia-bulk-specific-kmers.cutoff10to20-from-Laura_31_0_diff1.list.kmer.fq.sai.sort.tmp 5.Otolia-bulk-specific-kmers.cutoff10to20-from-Laura_31_0_diff1.list.kmer.fq.sai.tmp.bam -o 5.Otolia-bulk-specific-kmers.cutoff10to20-from-Laura_31_0_diff1.list.kmer.fq.sai.tmp.srt.bam &

samtools sort -T 6.Otolia-bulk-specific-kmers.cutoff10to20-from-Kuba_31_0_diff1.list.kmer.fq.sai.sort.tmp 6.Otolia-bulk-specific-kmers.cutoff10to20-from-Kuba_31_0_diff1.list.kmer.fq.sai.tmp.bam -o 6.Otolia-bulk-specific-kmers.cutoff10to20-from-Kuba_31_0_diff1.list.kmer.fq.sai.tmp.srt.bam &

samtools sort -T 7.Otolia-bulk-specific-kmers.cutoff10to20-from-none_31_0_diff1.list.kmer.fq.sai.sort.tmp 7.Otolia-bulk-specific-kmers.cutoff10to20-from-none_31_0_diff1.list.kmer.fq.sai.tmp.bam -o 7.Otolia-bulk-specific-kmers.cutoff10to20-from-none_31_0_diff1.list.kmer.fq.sai.tmp.srt.bam

```
### Counting kmers mapped to the individual chromosomes of the genome in 1MB windows
```sh
!/bin/bash
bedtool=/mnt/d/Softwares/bedtools
for m in $(ls *.tmp.srt.bam)
do
       $bedtool makewindows -b Chr01.bed -w 1000000 | $bedtool intersect -a - -b $m -c -bed > $m.Chr01.bins   #repeat the same if you want to extract the kmers for other chromosomes too
done

```

### Plot the kmers obtained for all the set algebra operation
```sh
#creating plots
import pandas as pd
import matplotlib.pyplot as plt

# Load multiple dataframes
dataframes = [pd.read_csv("merged_kmers.txt", sep="\t"),
              pd.read_csv("merged_2_otolia.txt", sep="\t")]
              #pd.read_csv("merged_2_otolia_3.txt", sep="\t")



# Create a list of colors for differentiating dataframes
colors = ['green', 'blue']

# Get unique chromosome names from the first dataframe (assuming they are the same for all dataframes)
chromosomes = dataframes[0]['chromosome'].unique()

# Calculate the number of rows and columns for subplots
num_rows = (len(chromosomes) + 3) // 3  # Ensure at least 3 per row
num_cols = min(len(chromosomes), 3)

# Create subplots
fig, axes = plt.subplots(num_rows, num_cols, figsize=(18, 4 * num_rows), sharex=True)

# Flatten the axes array for easy iteration
axes = axes.flatten()

# Iterate through each chromosome and plot the number of unique Kmers for each dataframe
for i, chromosome in enumerate(chromosomes):
    ax = axes[i]
    for j, dataframe in enumerate(dataframes):
        chromosome_data = dataframe[dataframe['chromosome'] == chromosome]
        ax.plot(chromosome_data['start'], chromosome_data['unique_kmers'], linestyle='-', color=colors[j], label=f'Dataframe {j+1}')
    ax.set_ylabel(f'Number of unique Kmers ({chromosome})')
    ax.set_xlabel('Genome Position')
    ax.legend()

# Remove any empty subplots if the number of chromosomes is not a multiple of 4
for i in range(len(chromosomes), num_rows * num_cols):
    fig.delaxes(axes[i])

# Adjust layout for better spacing
plt.tight_layout()

# Show the plots
plt.show()
```

### Variant calling, preprocessing and quality filtering
```sh
!/bin/bash

for i in $(ls *.srt.bam)
do
freebayes -f $reference --ploidy 4 -C 5 -m 30 -q 20 -S 0 -R 0 -F 0.2 -J -H -P 0.95 --use-best-n-alleles 4 ${i} > ${i}.vcf
done

wait
#normalizing vcf files
for i in $(ls *.vcf)
do
        bcftools norm -f $reference -o ${i}.norm ${i}
done

wait

#filtering normalized vcfs followed by annotation
for m in $(ls *.vcf.norm)
do
        bgzip ${m}
done

wait
for n in $(ls *.norm.gz)
do
        tabix ${n}
done

wait

#vcflib based filrering

wait
export PATH=/mnt/d/Softwares:$PATH


for l in $(ls *.norm.gz);
do

 vcfbiallelic | vcffilter -f "QUAL > 19 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1 & EPP > 2 & SRP > 2 & SAP > 2 & MQM / MQMR > 0.39 & MQM / MQMR < 2.51" -g "GQ > 20"  ${l} >  ${l}_filtered.vcf

done

```
### Variant Annotation
```sh
!/bin/bash
aaa
```
```diff
+ Note: Save the afforementioned script in a single bash file and use wait command after each step to run it as a pipeline. Properly set your input/output, installed softwares paths to avoid errors. If the computational power (RAM) of your server isn"t high, then run the scripts in chunks (i.e one step at a time) to avoid out of memory errors.

# References
1. https://github.com/cprodhom/CoSSA-workflows





