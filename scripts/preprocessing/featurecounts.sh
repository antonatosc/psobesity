#!/bin/bash
cd ~/Psobesity/SRP055813.txt
featureCounts -T 6 \
-a ~/Psobesity/Indexes/Homo_sapiens.GRCh38.107.gtf \
-o ~/Psobesity/SRP055813.txt/final_count \
~/Psobesity/SRP055813.txt/*.finalAligned.out.bam

cd ~/Psobesity/SRP065812.txt
featureCounts -T 6 -p -s 1 \
--countReadPairs \
-a ~/Psobesity/Indexes/Homo_sapiens.GRCh38.107.gtf \
-o ~/Psobesity/SRP065812.txt/final_count \
~/Psobesity/SRP065812.txt/*.finalAligned.out.bam

cd ~/Psobesity/SRP132990.txt
featureCounts -T 6 -p \
--countReadPairs \
-a ~/Psobesity/Indexes/Homo_sapiens.GRCh38.107.gtf \
-o ~/Psobesity/SRP132990.txt/final_count \
~/Psobesity/SRP132990.txt/*.finalAligned.out.bam

cd ~/Psobesity/SRP145260.txt
featureCounts -T 6 -s 1 \
-a ~/Psobesity/Indexes/Homo_sapiens.GRCh38.107.gtf \
-o ~/Psobesity/SRP145260.txt/final_count \
~/Psobesity/SRP145260.txt/*.finalAligned.out.bam

cd ~/Psobesity/SRP165679.txt
featureCounts -T 6 -p -s 2 \
--countReadPairs \
-a ~/Psobesity/Indexes/Homo_sapiens.GRCh38.107.gtf \
-o ~/Psobesity/SRP165679.txt/final_count \
~/Psobesity/SRP165679.txt/*.finalAligned.out.bam

cd ~/Psobesity/SRP268322.txt
featureCounts -T 6 -s 1 \
-a ~/Psobesity/Indexes/Homo_sapiens.GRCh38.107.gtf \
-o ~/Psobesity/SRP268322.txt/final_count \
~/Psobesity/SRP268322.txt/*.finalAligned.out.bam

cd ~/Psobesity/SRP278883.txt
featureCounts -T 6 -p -s 2 \
--countReadPairs \
-a ~/Psobesity/Indexes/Homo_sapiens.GRCh38.107.gtf \
-o ~/Psobesity/SRP278883.txt/final_count \
~/Psobesity/SRP278883.txt/*.finalAligned.out.bam

cd ~/Psobesity/SRP295864.txt
featureCounts -T 6 \
-a ~/Psobesity/Indexes/Homo_sapiens.GRCh38.107.gtf \
-o ~/Psobesity/SRP295864.txt/final_count \
~/Psobesity/SRP295864.txt/*.finalAligned.out.bam

cd ~/Psobesity/SRP304398.txt
featureCounts -T 6 -p -s 1 \
--countReadPairs \
-a ~/Psobesity/Indexes/Homo_sapiens.GRCh38.107.gtf \
-o ~/Psobesity/SRP304398.txt/final_count \
~/Psobesity/SRP304398.txt/*.finalAligned.out.bam
