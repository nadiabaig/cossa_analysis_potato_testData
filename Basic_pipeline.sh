tool="/home/baign/anaconda3/bin/glistmaker"

fastq1="/1data/Nadia/10X_Genomics_data/10XGenomics_data/batch2/usftp1.novogene.com/H202SC19122466/raw_data/HMK3HDSXX/Otolia/Otolia_S3_L001_R1_001.fastq"
fastq2="/1data/Nadia/10X_Genomics_data/10XGenomics_data/batch2/usftp1.novogene.com/H202SC19122466/raw_data/HMK3HDSXX/Otolia/Otolia_S3_L001_R2_001.fastq"

#$tool $fastq1 -w 31 -o otolia.all.read1
#$tool $fastq2 -w 31 -o otolia.all.read2

#/home/baign/anaconda3/bin/glistcompare otolia.all.read1_31.list otolia.all.read2_31.list -u -o otolia.all
#/home/baign/anaconda3/bin/glistcompare otolia.all_31_union.list otolia.all_31_union.list -i -c 20 -o otolia.all.cutoff2
#wait
#/home/baign/anaconda3/bin/glistquery otolia.all.cutoff2_31_intrsec.list > kmer.file.kmer
#now keep the kmers that are shared between both datasets
#kmers present only in Otolia
#kmers present in the dataset downloaded

#python /home/baign/anaconda3/bin/kmer_to_fastq.py kmer.file.kmer > kmer.file.fq

#comparing two files
#1) Otolia specific kmers not present in SenSephir
#python optimized.py
#wait
#conversion of kmers to fastq and further analysis 
#python /home/baign/anaconda3/bin/kmer_to_fastq.py Otolia_only.txt > Otolia_specific.fq
#python /home/baign/anaconda3/bin/kmer_to_fastq.py Sen3Sephir_only.txt > Sen3Sephir_only.fq
#python /home/baign/anaconda3/bin/kmer_to_fastq.py common_otolia_sen3.txt > common_otolia_sen3.fq

#wait
reference="/1data/Nadia/Project1_SNP_array/Potato_samples_processing/refdata-Agria_assembly_final_2020_21_08/Agria_assembly_final_2020_21_08.fasta
"

#bwa index $reference
#bwa aln $reference Otolia_specific.fq > Otolia_specific.fq.sai
#bwa aln $reference Sen3Sephir_only.fq > Sen3Sephir_only.fq.sai
#bwa aln $reference common_otolia_sen3.fq > common_otolia_sen3.fq.sai

#wait
#bwa samse $reference Sen3Sephir_only.fq.sai Sen3Sephir_only.fq > Sen3Sephir_only.fq.samse.tmp
#wait
#samtools view -Sbh Sen3Sephir_only.fq.samse.tmp -o Sen3Sephir_only.fq.samse.bam
#wait
#samtools sort -T Sen3Sephir_only.fq.samse.sort.tmp Sen3Sephir_only.fq.samse.bam -o Sen3Sephir_only.fq.samse.srt.bam
#wait
#rm Sen3Sephir_only.fq.samse.tmp Sen3Sephir_only.fq.samse.bam


#bwa samse $reference common_otolia_sen3.fq.sai common_otolia_sen3.fq > common_otolia_sen3.fq.samse.tmp
#wait
#samtools view -Sbh common_otolia_sen3.fq.samse.tmp -o common_otolia_sen3.fq.samse.bam
#wait
#samtools sort -T common_otolia_sen3.fq.samse.sort.tmp common_otolia_sen3.fq.samse.bam -o common_otolia_sen3.fq.samse.srt.bam
#wait
#rm common_otolia_sen3.fq.samse.tmp common_otolia_sen3.fq.samse.bam
#wait
#for i in $(ls *.fq.samse.srt.bam)
#do
#	samtools index ${i} 
#done
#wait

#extracting Chr10 specific reads only aS sEN3 Gene is on Chr10
#samtools view -b Otolia_specific.fq.samse.srt.bam Chr10 >  Otolia_specific_Chr10.bam
#samtools view -b Sen3Sephir_only.fq.samse.srt.bam Chr10 > Sen3Sephir_only_Chr10.bam
#samtools view -b common_otolia_sen3.fq.samse.srt.bam Chr10 > common_otolia_Chr10.bam
#wait

#counting how many kmers are mapped to the reference, how many kmers mapped to each 1MB chromosome bin

#bedtools makewindows -b chr10.bed -w 1000000 | bedtools intersect -a - -b Otolia_specific_Chr10.bam -c -bed > Otolia_specific.chr10.txt
#bedtools makewindows -b chr10.bed -w 1000000 | bedtools intersect -a - -b Sen3Sephir_only_Chr10.bam -c -bed > Sen3Sephir_specific.chr10.txt
#bedtools makewindows -b chr10.bed -w 1000000 | bedtools intersect -a - -b common_otolia_Chr10.bam -c -bed > Common_otolia_sen3.chr10.txt
#wait
#freebayes -f $reference --ploidy 4 -C 5 -m 30 -q 20 -S 0 -R 0 -F 0.2 -J -H -P 0.95 --use-best-n-alleles 4 Otolia_specific_Chr10.bam > Otolia_specific_Chr10.vcf
#freebayes -f $reference --ploidy 4 -C 5 -m 30 -q 20 -S 0 -R 0 -F 0.2 -J -H -P 0.95 --use-best-n-alleles 4 Sen3Sephir_only_Chr10.bam > Sen3Sephir_only_Chr10.vcf
#freebayes -f $reference --ploidy 4 -C 5 -m 30 -q 20 -S 0 -R 0 -F 0.2 -J -H -P 0.95 --use-best-n-alleles 4 common_otolia_Chr10.bam > common_otolia_Chr10.vcf
#wait
#normalizing vcfs/getting biallelic
#for i in $(ls *_Chr10.vcf)
#do
#	bcftools norm -f $reference -o ${i}.norm  ${i}
#done

#wait
#annotation

#java -jar /1data/Nadia/10X_Genomics_data/10XGenomics_data/batch2/usftp1.novogene.com/H202SC19122466/raw_data/HMK3HDSXX/Otolia/kmer_analysis/snpEff/snpEff.jar download Solanum_tuberosum

#wait

for m in $(ls *.norm)
do

java -jar /1data/Nadia/10X_Genomics_data/10XGenomics_data/batch2/usftp1.novogene.com/H202SC19122466/raw_data/HMK3HDSXX/Otolia/kmer_analysis/snpEff/snpEff.jar Solanum_tuberosum ${m} > ${m}_annotated.vcf
#wait

#java -jar /1data/Nadia/10X_Genomics_data/10XGenomics_data/batch2/usftp1.novogene.com/H202SC19122466/raw_data/HMK3HDSXX/Otolia/kmer_analysis/snpEff/snpEff.jar stats -v ${m}.log > ${m}.html
done

#keeping only snp,ins and del
#bcftools filter -i "TYPE='snp'" input.vcf > snp
