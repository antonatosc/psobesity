#!/bin/bash

#STAR \
#--runThreadN 6 \
#--runMode genomeGenerate \
#--genomeDir ~/Psobesity/Indexes \
#--sjdbGTFfile ~/Psobesity/Indexes/Homo_sapiens.GRCh38.107.gtf \
#--genomeFastaFiles ~/Psobesity/Indexes/Homo_sapiens.GRCh38.dna.primary_assembly.fa

#Load Genome
STAR \
--genomeLoad LoadAndExit \
--genomeDir ~/Indexes/GenomeDir

#get file names
cd ~/Psobesity/datasets
datasets="*.txt"
for f in $datasets
do
#get sample names
while read -r line
do
STAR \
--genomeDir ~/Psobesity/GenomeDir \
--genomeLoad LoadAndKeep \
--readFilesIn ~/Psobesity/$f/$line/*.fq.gz \
--runThreadN 6 \
--readFilesCommand zcat \
--outSAMtype BAM SortedByCoordinate \
--limitBAMsortRAM 10000000000 \
--outFilterMultimapNmax 1 \
--outFilterScoreMinOverLread 0 \
--outFilterMatchNminOverLread 0 \
--outBAMsortingBinsN 80 \
--outFileNamePrefix ~/Psobesity/$f/$line.final
echo  "Completed alignment of $line in project $f."
done < $f
echo "Completed alignment of $project $f."
done

STAR \
--genomeDir ~/Psobesity/GenomeDir \
--genomeLoad Remove
