#!/bin/bash
mkdir -p genome/assemblies genome/annotations genome/proteins
for full_id in $(cat asm-ids-full.txt);do
  gc_id=${full_id:0:15}
  asm_id=${full_id:16}
  gc_id_1=${gc_id:0:3}
  gc_id_2=${gc_id:4:3}
  gc_id_3=${gc_id:7:3}
  gc_id_4=${gc_id:10:3}
  #echo $gc_id_1 $gc_id_2 $gc_id_3 $gc_id_4
  url=https://ftp.ncbi.nlm.nih.gov/genomes/all/$gc_id_1/$gc_id_2/$gc_id_3/$gc_id_4/$full_id/${full_id}_genomic.fna.gz
  if [ -s genome/assemblies/${full_id}_genomic.fna ];then
     echo "$url already downloaded ."  
  else
     echo "downloading $url ..."
     wget -O genome/assemblies/${full_id}_genomic.fna.gz $url
     gunzip genome/assemblies/${full_id}_genomic.fna.gz
     samtools faidx genome/assemblies/${full_id}_genomic.fna
  fi
  url=https://ftp.ncbi.nlm.nih.gov/genomes/all/$gc_id_1/$gc_id_2/$gc_id_3/$gc_id_4/$full_id/${full_id}_genomic.gff.gz
  if [ -s genome/annotations/${full_id}_genomic.gff ];then
    echo "$url already downloaded ."
  else
    echo "downloading $url ..."
    wget -O genome/annotations/${full_id}_genomic.gff.gz $url
    gunzip genome/annotations/${full_id}_genomic.gff.gz
  fi
  url=https://ftp.ncbi.nlm.nih.gov/genomes/all/$gc_id_1/$gc_id_2/$gc_id_3/$gc_id_4/$full_id/${full_id}_protein.faa.gz
  if [ -s genome/proteins/${full_id}_protein.faa ];then
    echo "$url already downloaded ."
  else
    echo "downloading $url ..."
    wget -O genome/proteins/${full_id}_protein.faa.gz  $url
    gunzip genome/proteins/${full_id}_protein.faa.gz
  fi
done
