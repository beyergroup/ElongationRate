#!/bin/bash -l

#SBATCH --account=cdebes
#SBATCH --error=/data/public/cdebes/star_circ.log
#SBATCH --job-name=star_circ
#SBATCH --partition=all
#SBATCH --ntasks=32

module load STAR-2.5.2b
module unload R-3.2.3
module load R-3.3.1

path=$1
filea=$2

filea1=$filea\_1.fq.gz
filea2=$filea\_2.fq.gz

#filea1=$filea\_1.fastq.gz
#filea2=$filea\_2.fastq.gz

#filea1=$filea\_1_sequence.fq.gz
#filea2=$filea\_2_sequence.fq.gz

#echo $path$filea1
#echo $path$filea2

#filea1=$filea\_R1_001_f.fastq.gz
#filea2=$filea\_R2_001_f.fastq.gz

#filea1=$filea\.fastq.gz
#filea1=$filea\.fq.gz
#filea1=$filea\.fastq

#filea1=$filea\_1.fastq
#filea2=$filea\_2.fastq

#mkdir $path/$filea/
#java -jar /data/public/cdebes/soft/Trimmomatic-0.33/trimmomatic-0.33.jar PE $path$filea1 $path$filea2 $path$filea/output_forward_paired.fq.gz $path$filea/output_forward_unpaired.fq.gz $path$filea/output_reverse_paired.fq.gz $path$filea/output_reverse_unpaired.fq.gz ILLUMINACLIP:../../soft/Trimmomatic-0.33/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:45

#java -jar /data/public/cdebes/soft/Trimmomatic-0.33/trimmomatic-0.33.jar PE $path$filea1 $path$filea2 $path$filea/output_forward_paired.fq.gz $path$filea/output_forward_unpaired.fq.gz $path$filea/output_reverse_paired.fq.gz $path$filea/output_reverse_unpaired.fq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:65

#java -jar /data/public/cdebes/soft/Trimmomatic-0.33/trimmomatic-0.33.jar SE $path$filea1 $path$filea/output_forward_paired.fq.gz ILLUMINACLIP:../../soft/Trimmomatic-0.33/adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:45

#STAR --runThreadN 54 --genomeDir /data/public/mgarmha1/Assemblies/Mouse/Starindex_75 --readFilesIn $path$filea1 $path$filea2 --outFileNamePrefix $path/$filea/ --readFilesCommand zcat --outSAMtype BAM Unsorted
#--outWigType bedGraph --outWigNorm RPM

#mkdir /data/public/cdebes/workspace/scripts/hsapiens_nascent/$filea\_circ/
#STAR --genomeDir /data/public/cdebes/workspace/genomes/hsapiens/ --runThreadN 32 --readFilesIn $path$filea/output_forward_paired.fq.gz $path$filea/output_reverse_paired.fq.gz --outFileNamePrefix /data/public/cdebes/workspace/scripts/hsapiens_nascent/$filea\_circ/ --chimSegmentMin 15 --readFilesCommand zcat --outSJfilterOverhangMin 15 15 15 15 --alignSJoverhangMin 15 --alignSJDBoverhangMin 15 --seedSearchStartLmax 30 --outFilterMultimapNmax 20 --outFilterScoreMin 1 --outFilterMatchNmin 1 --outFilterMismatchNmax 2 --chimScoreMin 15 --chimScoreSeparation 10 --chimJunctionOverhangMin 15
#mkdir /data/public/cdebes/workspace/scripts/hsapiens_nascent/$filea\_circ1/

#STAR --runThreadN 32 --genomeDir [genome] --outSAMtype None --readFilesIn Sample1_1.fastq.gz 
#       --readFilesCommand zcat \
#       --outFileNamePrefix [sample prefix] \
#       --outReadsUnmapped Fastx \
#       --outSJfilterOverhangMin 15 15 15 15 \
#       --alignSJoverhangMin 15 \
#       --alignSJDBoverhangMin 15 \
#       --seedSearchStartLmax 30 \
#       --outFilterMultimapNmax 20 \
#       --outFilterScoreMin 1 \
#       --outFilterMatchNmin 1 \
#       --outFilterMismatchNmax 2 \
#       --chimSegmentMin 15 \
#       --chimScoreMin 15 \
#       --chimScoreSeparation 10 \
#       --chimJunctionOverhangMin 15 \

#TO DO
mkdir /data/public/cdebes/workspace/scripts/hsapiens_prog/$filea\_circ/
STAR --genomeDir /data/public/cdebes/workspace/genomes/hsapiens/ --runThreadN 24 --readFilesIn $path$filea/output_forward_paired.fq.gz $path$filea/output_reverse_paired.fq.gz --outFileNamePrefix /data/public/cdebes/workspace/scripts/hsapiens_prog/$filea\_circ/ --outSAMtype BAM Unsorted SortedByCoordinate --chimSegmentMin 15 --readFilesCommand zcat --outSJfilterOverhangMin 15 15 15 15 --alignSJoverhangMin 15 --alignSJDBoverhangMin 15 --seedSearchStartLmax 30 --outFilterMultimapNmax 20 --outFilterScoreMin 1 --outFilterMatchNmin 1 --outFilterMismatchNmax 2 --chimScoreMin 15 --chimScoreSeparation 10 --chimJunctionOverhangMin 15

#mkdir /data/public/cdebes/workspace/scripts/hsapiens/$filea\_circ/
#STAR --genomeDir /data/public/cdebes/workspace/genomes/hsapiens/ --runThreadN 24 --readFilesIn $path$filea/output_forward_paired.fq.gz  --outFileNamePrefix /data/public/cdebes/workspace/scripts/fuchs_labeling/$filea --outSAMtype BAM Unsorted SortedByCoordinate --sjdbGTFfile /data/public/cdebes/workspace/genomes/hsapiens/Homo_sapiens.GRCh38.85.gtf --sjdbGTFchrPrefix chr --readFilesCommand zcat --outFilterMatchNmin 20 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMismatchNmax 35 --outWigType bedGraph --outWigNorm None

#mkdir /data/public/cdebes/workspace/scripts/dmelanogaster_mut/$filea\_circ/
#STAR --genomeDir /data/public/cdebes/workspace/genomes/dmelanogaster/ --runThreadN 32 --readFilesIn $path$filea/output_forward_paired.fq.gz $path$filea/output_reverse_paired.fq.gz  --outFileNamePrefix /data/public/cdebes/workspace/scripts/dmelanogaster_mut/$filea\_circ/ --sjdbGTFchrPrefix chr --readFilesCommand zcat --chimSegmentMin 15 --outSJfilterOverhangMin 15 15 15 15 --alignSJoverhangMin 15 --alignSJDBoverhangMin 15 --seedSearchStartLmax 30 --outFilterMultimapNmax 20 --outFilterScoreMin 1 --outFilterMatchNmin 1 --outFilterMismatchNmax 2 --chimScoreMin 15 --chimScoreSeparation 10 --chimJunctionOverhangMin 15

#mkdir /data/public/cdebes/workspace/scripts/dmelanogaster_merged/$filea\_circ/
#STAR --genomeDir /data/public/cdebes/workspace/genomes/dmelanogaster/ --runThreadN 24 --readFilesIn  $path$filea.fastq.gz --outFileNamePrefix /data/public/cdebes/workspace/scripts/dmelanogaster_merged/$filea\_circ/ --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --outWigType bedGraph --outFilterType BySJout --outWigNorm None --chimSegmentMin 15 --outSJfilterOverhangMin 15 15 15 15 --alignSJoverhangMin 15 --alignSJDBoverhangMin 15 --seedSearchStartLmax 30 --outFilterMultimapNmax 20 --outFilterScoreMin 1 --outFilterMatchNmin 1 --outFilterMismatchNmax 2 --chimScoreMin 15 --chimScoreSeparation 10 --chimJunctionOverhangMin 15

#mkdir /data/public/cdebes/workspace/scripts/rnorvegicus/$filea\_circ/
#STAR --genomeDir /data/public/cdebes/workspace/genomes/rnorvegicus/ --runThreadN 32 --readFilesIn $path$filea.fastq  --outFileNamePrefix /data/public/cdebes/workspace/scripts/rnorvegicus/$filea\_circ/ --outSAMtype BAM Unsorted SortedByCoordinate --outWigType bedGraph  --outFilterType BySJout --outWigNorm None --chimSegmentMin 15 --outSJfilterOverhangMin 15 15 15 15 --alignSJoverhangMin 15 --alignSJDBoverhangMin 15 --seedSearchStartLmax 30 --outFilterMultimapNmax 20 --outFilterScoreMin 1 --outFilterMatchNmin 1 --outFilterMismatchNmax 2 --chimScoreMin 15 --chimScoreSeparation 10 --chimJunctionOverhangMin 15

#STAR --genomeDir /data/public/cdebes/workspace/genomes/mmusculus_ensembl/ --runThreadN 50 --readFilesIn $path$filea1  $path$filea2 --outFileNamePrefix /data/public/cdebes/workspace/scripts/mmusculus_hypoxia/$filea --outSAMtype BAM Unsorted SortedByCoordinate --sjdbGTFfile /data/public/cdebes/workspace/genomes/mmusculus/Mus_musculus.GRCm38.87.gtf --sjdbGTFchrPrefix chr --readFilesCommand zcat --outWigType bedGraph  --outWigNorm RPM

#mkdir /data/public/cdebes/workspace/scripts/mmusculus_hypothalamus/$filea\_circ/
#STAR --genomeDir /data/public/cdebes/workspace/genomes/mmusculus/ --runThreadN 32 --readFilesIn $path$filea.fastq.gz --outFileNamePrefix /data/public/cdebes/workspace/scripts/mmusculus_hypothalamus/$filea\_circ/  --outSAMtype BAM Unsorted SortedByCoordinate --sjdbGTFfile /data/public/cdebes/workspace/genomes/mmusculus/Mus_musculus.GRCm38.87.gtf --sjdbGTFchrPrefix chr --readFilesCommand zcat --outWigType bedGraph  --outFilterType BySJout --outWigNorm None --chimSegmentMin 15 --outSJfilterOverhangMin 15 15 15 15 --alignSJoverhangMin 15 --alignSJDBoverhangMin 15 --seedSearchStartLmax 30 --outFilterMultimapNmax 20 --outFilterScoreMin 1 --outFilterMatchNmin 1 --outFilterMismatchNmax 2 --chimScoreMin 15 --chimScoreSeparation 10 --chimJunctionOverhangMin 15

#STAR --genomeDir /data/public/cdebes/workspace/genomes/mmusculus/ --runThreadN 32 --readFilesIn $path$filea/output_forward_paired.fq.gz $path$filea/output_reverse_paired.fq.gz  --outFileNamePrefix /data/public/cdebes/workspace/scripts/mmusculus_kidney_CR/$filea --readFilesCommand zcat --sjdbGTFfile /data/public/cdebes/workspace/genomes/mmusculus/Mus_musculus.GRCm38.87.gtf --sjdbGTFchrPrefix chr --outSAMtype BAM Unsorted SortedByCoordinate --readFilesCommand zcat --outWigType bedGraph   --outFilterType BySJout --outWigNorm None --outSJfilterOverhangMin 10 10 10 10 --outSJfilterCountUniqueMin 1 1 1 1 --outSJfilterIntronMaxVsReadN 400000 400000 400000 400000 --outSJfilterDistToOtherSJmin 0 0 0 0

#STAR --genomeDir /data/public/cdebes/workspace/genomes/mmusculus/ --runThreadN 32 --readFilesIn $path$filea/output_forward_paired.fq.gz $path$filea/output_reverse_paired.fq.gz  --outFileNamePrefix /data/public/cdebes/workspace/scripts/mmusculus/$filea --readFilesCommand zcat --sjdbGTFfile /data/public/cdebes/workspace/genomes/mmusculus/Mus_musculus.GRCm38.87_havana.gtf --sjdbGTFchrPrefix chr --outSAMtype BAM Unsorted SortedByCoordinate --outWigType bedGraph --outWigNorm None  --outFilterType BySJout --outSJfilterOverhangMin 10 10 10 10 --outSJfilterCountUniqueMin 1 1 1 1 --outSJfilterIntronMaxVsReadN 400000 400000 400000 400000 --outSJfilterDistToOtherSJmin 0 0 0 0

#mkdir /data/public/cdebes/workspace/scripts/mmusculus_aldr/$filea\_circ/
#STAR --genomeDir /data/public/cdebes/workspace/genomes/mmusculus/  --runThreadN 24 --readFilesIn $path$filea1 $path$filea2 --outFileNamePrefix /data/public/cdebes/workspace/scripts/mmusculus_aldr/$filea\_circ/ --readFilesCommand zcat --sjdbGTFchrPrefix chr --chimSegmentMin 15 --outSJfilterOverhangMin 15 15 15 15 --alignSJoverhangMin 15 --alignSJDBoverhangMin 15 --seedSearchStartLmax 30 --outFilterMultimapNmax 20 --outFilterScoreMin 1 --outFilterMatchNmin 1 --outFilterMismatchNmax 2 --chimScoreMin 15 --chimScoreSeparation 10 --chimJunctionOverhangMin 15

#mkdir /data/public/cdebes/workspace/scripts/mmusculus_kidney_CR/$filea\_circ/
#STAR --genomeDir /data/public/cdebes/workspace/genomes/mmusculus/ --runThreadN 24 --readFilesIn $path$filea/output_forward_paired.fq.gz $path$filea/output_reverse_paired.fq.gz --outFileNamePrefix /data/public/cdebes/workspace/scripts/mmusculus_kidney_CR/$filea\_circ/ --readFilesCommand zcat --sjdbGTFchrPrefix chr --chimSegmentMin 15 --outSJfilterOverhangMin 15 15 15 15 --alignSJoverhangMin 15 --alignSJDBoverhangMin 15 --seedSearchStartLmax 30 --outFilterMultimapNmax 20 --outFilterScoreMin 1 --outFilterMatchNmin 1 --outFilterMismatchNmax 2 --chimScoreMin 15 --chimScoreSeparation 10 --chimJunctionOverhangMin 15

#mkdir /data/public/cdebes/workspace/scripts/mmusculus/$filea\_circ/
#STAR --genomeDir /data/public/cdebes/workspace/genomes/mmusculus/  --runThreadN 24 --readFilesIn $path$filea/output_forward_paired.fq.gz $path$filea/output_reverse_paired.fq.gz --outFileNamePrefix /data/public/cdebes/workspace/scripts/mmusculus/$filea\_circ/ --readFilesCommand zcat --sjdbGTFchrPrefix chr --chimSegmentMin 15 --outSJfilterOverhangMin 15 15 15 15 --alignSJoverhangMin 15 --alignSJDBoverhangMin 15 --seedSearchStartLmax 30 --outFilterMultimapNmax 20 --outFilterScoreMin 1 --outFilterMatchNmin 1 --outFilterMismatchNmax 2 --chimScoreMin 15 --chimScoreSeparation 10 --chimJunctionOverhangMin 15

#STAR --genomeDir /data/public/cdebes/workspace/genomes/mmusculus/ --runThreadN 24 --readFilesIn $path$filea/output_forward_paired.fq.gz --outFileNamePrefix /data/public/cdebes/workspace/scripts/mmusculus_whiteadipose_polya/$filea --readFilesCommand zcat --sjdbGTFfile /data/public/cdebes/workspace/genomes/mmusculus/Mus_musculus.GRCm38.87.gtf --sjdbGTFchrPrefix chr --outWigType bedGraph --outWigNorm None --outSJfilterOverhangMin 10 10 10 10 --outSJfilterCountUniqueMin 1 1 1 1 --outSJfilterIntronMaxVsReadN 400000 400000 400000 400000 --outSJfilterDistToOtherSJmin 0 0 0 0 --outSAMstrandField intronMotif

#STAR --genomeDir /data/public/cdebes/workspace/genomes/mmusculus/ --runThreadN 24 --readFilesIn $path$filea1 $path$filea2 --outFileNamePrefix /data/public/cdebes/workspace/scripts/fromElifeLys/$filea --sjdbGTFfile /data/public/cdebes/workspace/genomes/mmusculus/Mus_musculus.GRCm38.83.gtf --sjdbGTFchrPrefix chr --outWigType bedGraph --outWigNorm None --outSJfilterOverhangMin 10 10 10 10 --outSJfilterCountUniqueMin 1 1 1 1 --outSJfilterIntronMaxVsReadN 400000 400000 400000 400000 --outSJfilterDistToOtherSJmin 0 0 0 0 --outSAMstrandField intronMotif

#mkdir /data/public/cdebes/workspace/scripts/celegans_mut/$filea\_circ/
#STAR --genomeDir /data/public/cdebes/workspace/genomes/celegans/ --runThreadN 32 --readFilesIn $path$filea/output_forward_paired.fq.gz $path$filea/output_reverse_paired.fq.gz  --outFileNamePrefix /data/public/cdebes/workspace/scripts/celegans_mut/$filea\_circ/ --sjdbGTFchrPrefix chr --readFilesCommand zcat --chimSegmentMin 15 --outSJfilterOverhangMin 15 15 15 15 --alignSJoverhangMin 15 --alignSJDBoverhangMin 15 --seedSearchStartLmax 30 --outFilterMultimapNmax 20 --outFilterScoreMin 1 --outFilterMatchNmin 1 --outFilterMismatchNmax 2 --chimScoreMin 15 --chimScoreSeparation 10 --chimJunctionOverhangMin 15

#mkdir /data/public/cdebes/workspace/scripts/celegans/$filea\_circ/
#STAR --genomeDir /data/public/cdebes/workspace/genomes/celegans/ --runThreadN 32 --readFilesIn $path$filea/output_forward_paired.fq.gz $path$filea/output_reverse_paired.fq.gz  --outFileNamePrefix /data/public/cdebes/workspace/scripts/celegans/$filea\_circ/ --sjdbGTFchrPrefix chr --readFilesCommand zcat --chimSegmentMin 15 --outSJfilterOverhangMin 15 15 15 15 --alignSJoverhangMin 15 --alignSJDBoverhangMin 15 --seedSearchStartLmax 30 --outFilterMultimapNmax 20 --outFilterScoreMin 1 --outFilterMatchNmin 1 --outFilterMismatchNmax 2 --chimScoreMin 15 --chimScoreSeparation 10 --chimJunctionOverhangMin 15

#STAR --genomeDir /data/public/cdebes/workspace/genomes/celegans/ --runThreadN 32 --readFilesIn $path$filea/output_forward_paired.fq.gz $path$filea/output_reverse_paired.fq.gz --outFileNamePrefix /data/public/cdebes/workspace/scripts/celegans/$filea --outSAMtype BAM Unsorted SortedByCoordinate --sjdbGTFfile /data/public/cdebes/workspace/genomes/celegans/Caenorhabditis_elegans.WBcel235.87.gtf --sjdbGTFchrPrefix chr --readFilesCommand zcat  --outWigType bedGraph --outWigNorm None --outSJfilterOverhangMin 10 10 10 10 --outSJfilterCountUniqueMin 1 1 1 1 --outSJfilterIntronMaxVsReadN 400000 400000 400000 400000 --outSJfilterDistToOtherSJmin 0 0 0 0 --outSAMstrandField intronMotif 

#STAR --genomeDir /data/public/cdebes/workspace/genomes/celegans/ --runThreadN 24 --readFilesIn $path$filea/output_forward_paired.fq.gz $path$filea/output_reverse_paired.fq.gz --outFileNamePrefix /data/public/cdebes/workspace/scripts/ssl1/$filea --outSAMtype BAM Unsorted SortedByCoordinate --sjdbGTFfile /data/public/cdebes/workspace/genomes/celegans/Caenorhabditis_elegans.WBcel235.87.gtf --sjdbGTFchrPrefix chr --readFilesCommand zcat  --outWigType bedGraph --outWigNorm None --outSAMstrandField intronMotif

#STAR --genomeDir /data/public/cdebes/workspace/genomes/celegans/ --runThreadN 24 --readFilesIn $path$filea/output_forward_paired.fq.gz $path$filea/output_reverse_paired.fq.gz --outFileNamePrefix /data/public/cdebes/workspace/scripts/celegans_mut/$filea --outSAMtype BAM Unsorted SortedByCoordinate --sjdbGTFfile /data/public/cdebes/workspace/genomes/celegans/Caenorhabditis_elegans.WBcel235.87.gtf --sjdbGTFchrPrefix chr --readFilesCommand zcat  --outWigType bedGraph --outWigNorm None --outSJfilterOverhangMin 10 10 10 10 --outSJfilterCountUniqueMin 1 1 1 1 --outSJfilterIntronMaxVsReadN 400000 400000 400000 400000 --outSJfilterDistToOtherSJmin 0 0 0 0 --outSAMstrandField intronMotif 

#STAR --genomeDir /data/public/cdebes/workspace/genomes/celegans/ --runThreadN 24 --readFilesIn $path$filea/output_forward_paired.fq.gz $path$filea/output_reverse_paired.fq.gz --outFileNamePrefix /data/public/cdebes/workspace/scripts/celegans_h2az/$filea --outSAMtype BAM Unsorted SortedByCoordinate --sjdbGTFfile /data/public/cdebes/workspace/genomes/celegans/Caenorhabditis_elegans.WBcel235.87.gtf --sjdbGTFchrPrefix chr --readFilesCommand zcat  --outWigType bedGraph --outWigNorm None --outSJfilterOverhangMin 10 10 10 10 --outSJfilterCountUniqueMin 1 1 1 1 --outSJfilterIntronMaxVsReadN 400000 400000 400000 400000 --outSJfilterDistToOtherSJmin 0 0 0 0 --outSAMstrandField intronMotif 


