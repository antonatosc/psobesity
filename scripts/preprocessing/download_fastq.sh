#!/bin/bash
cd ~/Psobesity/datasets
datasets="*.txt" #Datasets and included samples are provided at supplementary table 1

for f in $datasets
do
while read -r line
do
prefetch $line \
--output-directory ~/Psobesity/$f \
--verify yes \
--max-size u \
-p
echo "$line of the $f dataset has been downloaded at ~/Psobesity/$f/$line"
done < $f
echo "Completed downloading $f dataset."
done

for f in $datasets
do
while read -r line
do
cd ~/Psobesity/$f/$line
fasterq-dump $line.sra \
--threads 6
rm $line.sra
pigz *.fastq \
--processes 6 
cd ~/Psobesity/datasets
echo "Converted $line of $f dataset into fastq.gz file."
done < $f
echo "Completed converting $f dataset."
done