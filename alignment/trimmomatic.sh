#Script used for treaming the raw fasta files. Raw reads were trimmed with trimmomatic version 0.33.

java -jar trimmomatic-0.33.jar PE $path$filea1 $path$filea2 $path$filea/output_forward_paired.fq.gz $path$filea/output_forward_unpaired.fq.gz $path$filea/output_reverse_paired.fq.gz $path$filea/output_reverse_unpaired.fq.gz ILLUMINACLIP:../../soft/Trimmomatic-0.33/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:65
