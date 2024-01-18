#!/bin/bash
####single
cd ~/Psobesity/datasets/single
datasets="*.txt"
for f in $datasets
do
while read -r line
do
trim_galore --fastqc ~/Psobesity/$f/$line/*.gz
done < $f
done
#####paired
cd ~/Psobesity/datasets/paired
datasets="*.txt"
for f in $datasets
do
while read -r line
do
trim_galore --paired --fastqc ~/Psobesity/$f/$line/*.gz
done < $f
done
