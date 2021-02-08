#!/bin/bash
#SBATCH --cpus-per-task=8

# module load proteinortho
# module load cufflinks

### Code Chunk 1 - starts

#mapping samples to ANT genome

#index the ant genome, in parent dir
#hisat2-build -f camp_genome.fna camp_index
#DONE

#fresh ant index, with exons
#simply in console:
#module load cufflinks
#gffread "/home/billu/cflo_genome/GCF_003227725.1_Cflo_v7.5_genomic.gff" -T -o camp.gtf
#module load hisat2
#hisat2_extract_splice_sites.py "/home/billu/cflo_genome/camp.gtf" > camp_splicesites.txt
#hisat2_extract_exons.py "/home/billu/cflo_genome/camp.gtf" > camp_exons.txt
#DONE

#echo "BUILDING THE ANT INDEX WITH SPLICE AND EXONS SITE"
#index the exon site version of the ANT genome
#done, resubmitting mapping
#hisat2-build -f --ss camp_splicesites.txt --exon camp_exons.txt ./camp_genome.fna camp_exons_index
#DONE

### chunk ends 


### Code Chunk 2 - starts
#map timecourse rnaseq samples to (indexed) genome
#Foragers output --> [id]_[timecourse#]_{run#].sam
#hisat2 -p 8 --new-summary -q --dta-cufflinks -x /home/billu/cflo_genome/camp_exons_index -U /home/billu/TC5_seq_reads/Foragers/Trimmed-2F_Galaxy96.fastq.gz -S /home/billu/TC5_seq_reads/2F_time5_1.sam

### chunk ends

### Code Chunk 3 - starts
module load samtools
#using .bam from hisat2 mapping>samtools sort and convert and Genbanks's supplied .gff3 annotations. 
#NO CUFFLINKS step used (nor cuffquant).
#-b ant genome
#sort and convert Hisat2 sams to sorted bams
samtools sort -o ./2F_time5_1.bam ./2F_time5_1.sam
samtools sort -o ./4F_time5_1.bam ./4F_time5_1.sam

cuffdiff -o Foragers -b "/home/billu/cflo_genome/camp_genome.fna" -p 16 -u "/home/billu/cflo_genome/GCF_003227725.1_Cflo_v7.5_genomic.gff" \2F_time5_1.bam \4F_time5_1.bam

### chunk ends