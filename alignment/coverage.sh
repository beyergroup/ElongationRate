#!/bin/bash -l
#SBATCH --account=cdebes
#SBATCH --error=/data/public/cdebes/log_mm.log
#SBATCH --job-name=mm_mpi
#SBATCH --ntasks=2
#SBATCH --partition=single

module unload bedtools-2.25.0
module load bedtools-2.25.0-master

#path=/cellnet/SyBACol/data_time_series/time_series/hsapiens/
#genome=/data/public/cdebes/workspace/genomes/hsapiens/

#FILES=($(find $path -name SN*.bam))

#for f in ${FILES[@]:0:6};
#do
#    o=(`basename $f | cut -d. -f1`)
#    samtools view -bSq 50 $path$o/$o\.sorted.bam > $path$o/$o\.prim.bam
#    bedtools genomecov  -ibam $path$o/$o\.prim.bam -g $genome/m.genome -strand - -bga -du -fs > $path$o/$o\_plus.bed
#    bedtools genomecov  -ibam $path$o/$o\.prim.bam -g $genome/m.genome -strand + -bga -du -fs > $path$o/$o\_minus.bed
#done

path=/cellnet/SyBACol/data_time_series/time_series/mmusculus/
path=/cellnet/SyBACol/data_time_series/time_series/mmusculus_kidney_CR/
genome=/data/public/cdebes/workspace/genomes/mmusculus/

FILES=($(find $path -name *.bam))

for k in ${FILES[*]}
do
    echo $k
    #samtools sort $k -o $k\.sorted
    samtools index $k\.sorted
done

for f in ${FILES[@]:0:8};
do
     #o=(`basename $f | cut -d. -f1`)
     #samtools view -bSq 50 $path/$o\.bam > $path/$o\.prim.bam
     bedtools genomecov  -split -ibam $f\.sorted -g $genome/m.genome -strand - -bga -du -fs > $k\_plus_splic.bed 
done


for f in ${FILES[@]:0:8};
do
     #o=(`basename $f | cut -d. -f1`)
     #samtools view -bSq 50 $path/$o\.bam > $path/$o\.prim.bam
     bedtools genomecov  -split -ibam $f\.sorted -g $genome/m.genome -strand + -bga -du -fs > $k\_minus_splic.bed 
done


for f in ${FILES[@]:0:8};
do
     o=(`basename $f | cut -d. -f1`)
     #samtools view -bSq 50 $path/$o\.bam > $path/$o\.prim.bam
     bedtools genomecov  -split -ibam $path/$o/$o\.sorted -g $genome/m.genome -strand - -bga -du -fs > $path/$o\_plus_splic.bed 
done


for f in ${FILES[@]:0:8};
do
     o=(`basename $f | cut -d. -f1`)
     #samtools view -bSq 50 $path/$o\.bam > $path/$o\.prim.bam
     bedtools genomecov  -split -ibam $path/$o/$o\.sorted -g $genome/m.genome -strand - -bga -du -fs > $path/$o\_plus_splic.bed 
done&

for f in ${FILES[@]:0:8};
do
    o=(`basename $f | cut -d. -f1`)
    bedtools genomecov  -split -ibam $path/$o/$o\.sorted.bam -g $genome/m.genome -strand + -bga -du -fs > $path/$o\_minus_splic.bed
done&

path=/cellnet/SyBACol/data_time_series/time_series/drosophila_d10/FastQ/
genome=/data/public/cdebes/workspace/genomes/dmelanogaster/

FILES=($(find $path -name *.sorted.bam))

for f in ${FILES[@]:0:6};
do
     o=(`basename $f | cut -d. -f1`)
     #samtools view -bSq 50 $path/$o\.bam > $path/$o\.prim.bam
     bedtools genomecov  -split -ibam $path/$o/$o\.sorted.bam -g $genome/m.genome -strand - -bga > $path/$o\_plus_splic.bed     
done&

for f in ${FILES[@]:0:6};
do
    o=(`basename $f | cut -d. -f1`)
    bedtools genomecov  -split -ibam $path/$o/$o\.sorted.bam -g $genome/m.genome -strand + -bga  > $path/$o\_minus_splic.bed
done&

for f in ${FILES[@]:7:16};
do
     o=(`basename $f | cut -d. -f1`)
     #samtools view -bSq 50 $path/$o\.bam > $path/$o\.prim.bam
     bedtools genomecov  -split -ibam $path/$o/$o\.sorted.bam -g $genome/m.genome -strand - -bga > $path/$o\_plus_splic.bed     
done&

for f in ${FILES[@]:7:16};
do
    o=(`basename $f | cut -d. -f1`)
    bedtools genomecov  -split -ibam $path/$o/$o\.sorted.bam -g $genome/m.genome -strand + -bga  > $path/$o\_minus_splic.bed
done&


for f in ${FILES[@]:17:23};
do
     o=(`basename $f | cut -d. -f1`)
     #samtools view -bSq 50 $path/$o\.bam > $path/$o\.prim.bam
     bedtools genomecov  -split -ibam $path/$o/$o\.sorted.bam -g $genome/m.genome -strand - -bga > $path/$o\_plus_splic.bed     
done&

for f in ${FILES[@]:17:23};
do
    o=(`basename $f | cut -d. -f1`)
    bedtools genomecov  -split -ibam $path/$o/$o\.sorted.bam -g $genome/m.genome -strand + -bga  > $path/$o\_minus_splic.bed
done&

for f in ${FILES[@]:24:30};
do
     o=(`basename $f | cut -d. -f1`)
     #samtools view -bSq 50 $path/$o\.bam > $path/$o\.prim.bam
     bedtools genomecov  -split -ibam $path/$o/$o\.sorted.bam -g $genome/m.genome -strand - -bga > $path/$o\_plus_splic.bed     
done&

for f in ${FILES[@]:24:30};
do
    o=(`basename $f | cut -d. -f1`)
    bedtools genomecov  -split -ibam $path/$o/$o\.sorted.bam -g $genome/m.genome -strand + -bga  > $path/$o\_minus_splic.bed
done&

for f in ${FILES[@]:31:52};
do
     o=(`basename $f | cut -d. -f1`)
     #samtools view -bSq 50 $path/$o\.bam > $path/$o\.prim.bam
     bedtools genomecov  -split -ibam $path/$o/$o\.sorted.bam -g $genome/m.genome -strand - -bga > $path/$o\_plus_splic.bed     
done&

for f in ${FILES[@]:31:52};
do
    o=(`basename $f | cut -d. -f1`)
    bedtools genomecov  -split -ibam $path/$o/$o\.sorted.bam -g $genome/m.genome -strand + -bga  > $path/$o\_minus_splic.bed
done&
