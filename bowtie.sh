#!/bin/bash -l

#SBATCH --account=cdebes
#SBATCH --error=/data/public/cdebes/bowtie.log
#SBATCH --job-name=bowtie
#SBATCH --partition=all
#SBATCH --ntasks=32


module load bowtie2-2.3.4.1

path=$1
filea=$2

#filea1=$filea\_1_sequence.fq.gz
#filea2=$filea\_2_sequence.fq.gz

filea1=$filea\_1.fq.gz
filea2=$filea\_2.fq.gz

#filea1=$filea\.fastq.gz
#filea2=$filea\_2.fastq.gz
#mkdir $path$filea
#java -jar /data/public/cdebes/soft/Trimmomatic-0.33/trimmomatic-0.33.jar PE $path$filea1 $path$filea2 $path$filea/output_forward_paired.fq.gz $path$filea/output_forward_unpaired.fq.gz $path$filea/output_reverse_paired.fq.gz $path$filea/output_reverse_unpaired.fq.gz ILLUMINACLIP:../../soft/Trimmomatic-0.33/adapters/TruSeq3-PE.fa:2:30:10 LEADING:2 TRAILING:2 SLIDINGWINDOW:4:15 MINLEN:55

#bowtie2 --threads 32 -x /data/public/cdebes/workspace/genomes/celegans/Caenorhabditis_elegans/UCSC/ce10/Sequence/Bowtie2Index/genome --fr --no-discordant -1 $path$filea/output_forward_paired.fq.gz -2 $path$filea/output_reverse_paired.fq.gz -S /data/public/cdebes/workspace/scripts/ssl1/$filea.sam 2> /data/public/cdebes/workspace/scripts/ssl1/$filea.log
#/data/programs/bowtie2-2.3.4.1/bowtie2 --threads 32 -x  /data/public/cdebes/workspace/genomes/celegans/bowtie_ensembl_WBcel235  --fr --no-discordant -1 $path$filea/output_forward_paired.fq.gz -2 $path$filea/output_reverse_paired.fq.gz -S /data/public/cdebes/workspace/scripts/ssl1/$filea.sam 2> /data/public/cdebes/workspace/scripts/ssl1/$filea.WBcel235.log
/data/public/cdebes/soft/samtools-1.9/samtools view -@ 32 -bS  /data/public/cdebes/workspace/scripts/ssl1/$filea.sam > /data/public/cdebes/workspace/scripts/ssl1/$filea.WBcel235.bam
rm /data/public/cdebes/workspace/scripts/ssl1/$filea.WBcel235.sam
/data/public/cdebes/soft/samtools-1.9/samtools view -@ 32 -F 0x04 -f 0x02 -q 30 -b /data/public/cdebes/workspace/scripts/ssl1/$filea.WBcel235.bam > /data/public/cdebes/workspace/scripts/ssl1/$filea.WBcel235.paired.bam
/data/public/cdebes/soft/samtools-1.9/samtools sort -@ 32 /data/public/cdebes/workspace/scripts/ssl1/$filea.WBcel235.paired.bam -o /data/public/cdebes/workspace/scripts/ssl1/$filea.sorted.WBcel235.paired.bam
#bedtools genomecov -ibam /data/public/cdebes/workspace/scripts/ssl1/$filea.sorted.paired.bam -g $genome/m.genome -bga -pc -fs > /data/public/cdebes/workspace/scripts/ssl1/$filea.bedgraph

#path=/data/public/cdebes/workspace/scripts/ssl1_chip_seq/
#list=`ls $path*.bam | sed 's!.*/!!' | sed 's!.bam!!' | uniq`
#for f in ${list[@]};
#do
#    samtools sort -@ 12 $path$f.bam -o $path$f.sorted.bam
#    samtools index $path$f.sorted.bam
#done

#/data/public/cdebes/soft/MACS-master/bin/macs2 peakcall --nomodel --shift 37 --extsize 73  -i /data/public/cdebes/workspace/scripts/ssl1_chip_seq/$filea.bedgraph -o /data/public/cdebes/workspace/scripts/ssl1_chip_seq/$filea.peaks
#bedtools bamtobed -i /data/public/cdebes/workspace/scripts/ssl1_chip_seq/$filea.bam > /data/public/cdebes/workspace/scripts/ssl1_chip_seq/$filea.bed
#python /data/public/cdebes/soft/danpos-2.2.2/danpos.py dpos  /data/public/cdebes/workspace/scripts/$spe/$filea.bed -m 

#source /data/public/cdebes/myEnv/bin/activate
#python /data/public/cdebes/soft/danpos-2.2.2/danpos.py dpos  /data/public/cdebes/workspace/scripts/$spe/$filea.bed -m 1 -o $filea
#sort -k 1,1 -k2,2n /data/public/cdebes/workspace/scripts/$spe/$filea.bed > /data/public/cdebes/workspace/scripts/$spe/$filea.sorted.bed

#R1-RNA

# macs2 callpeak -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53943.bam /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53967.bam /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53991.bam -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53935.bam /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53959.bam /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53983.bam -f BAM -g ce --outdir /data/public/cdebes/workspace/scripts/ssl1_chip_seq/ -n RNA-WT1 --bdg  >& /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1.log &
# macs2 callpeak -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53945.bam /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53969.bam /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53993.bam -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53937.bam /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53961.bam /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53985.bam -f BAM -g ce --outdir /data/public/cdebes/workspace/scripts/ssl1_chip_seq/ -n RNA-SSL1 --bdg --nomodel --shift 37 --extsize 73 >& /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1.log &
# macs2 callpeak -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53947.bam /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53971.bam /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53995.bam -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53939.bam /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53963.bam /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53987.bam -f BAM -g ce --outdir /data/public/cdebes/workspace/scripts/ssl1_chip_seq/ -n RNA-WT1-UV --bdg --nomodel --shift 37 --extsize 73 >& /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-UV.log &
# macs2 callpeak -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53949.bam /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53973.bam /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53997.bam -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53941.bam /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53965.bam /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53989.bam -f BAM -g ce --outdir /data/public/cdebes/workspace/scripts/ssl1_chip_seq/ -n RNA-SSL1-UV --bdg --nomodel --shift 37 --extsize 73 >& /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-UV.log &

# macs2 callpeak -t /mnt/cellnet/SSL1/data/chip-seq/48_N2HTZC.bam -c /mnt/cellnet/SSL1/data/chip-seq/56_N2InputC.bam -f BAM -g ce --outdir /data/public/cdebes/workspace/scripts/ssl1_chip_seq/ -n HTZ-WT2 --bdg --nomodel --shift 37 --extsize 73 >& /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT2.log &

# macs2 callpeak -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53951.bam /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53975.bam /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53999.bam -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53935.bam /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53959.bam /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53983.bam -f BAM -g ce --outdir /data/public/cdebes/workspace/scripts/ssl1_chip_seq/ -n HTZ-WT1 --bdg --nomodel --shift 37 --extsize 73 >& /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1.log &

# macs2 callpeak -t /mnt/cellnet/SSL1/data/chip-seq/50_ssl1HTZC.bam -c /mnt/cellnet/SSL1/data/chip-seq/58_ssl1InputC.bam -f BAM -g ce --outdir /data/public/cdebes/workspace/scripts/ssl1_chip_seq/ -n HTZ-SSL2 --bdg --nomodel --shift 37 --extsize 73 >& /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL2.log &

# macs2 callpeak -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53953.bam /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53977.bam /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_54001.bam -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53937.bam /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53961.bam /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53985.bam -f BAM -g ce --outdir /data/public/cdebes/workspace/scripts/ssl1_chip_seq/ -n HTZ-SSL1 --bdg --nomodel --shift 37 --extsize 73 >& /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1.log &

# macs2 callpeak -t /mnt/cellnet/SSL1/data/chip-seq/49_N2HTZUV.bam  -c /mnt/cellnet/SSL1/data/chip-seq/57_N2InputUV.bam  -f BAM -g ce --outdir /data/public/cdebes/workspace/scripts/ssl1_chip_seq/ -n HTZ-WT1-UV2 --bdg --nomodel --shift 37 --extsize 73 >& /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV2.log &

# macs2 callpeak -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53955.bam /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53979.bam /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_54003.bam -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53939.bam /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53963.bam /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53987.bam -f BAM -g ce --outdir /data/public/cdebes/workspace/scripts/ssl1_chip_seq/ -n HTZ-WT1-UV --bdg --nomodel --shift 37 --extsize 73 >& /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV.log &

# macs2 callpeak -t /mnt/cellnet/SSL1/data/chip-seq/51_ssl1HTZUV.bam   -c /mnt/cellnet/SSL1/data/chip-seq/59_ssl1InputUV.bam  -f BAM -g ce --outdir /data/public/cdebes/workspace/scripts/ssl1_chip_seq/ -n HTZ-SSL1-UV2 --bdg --nomodel --shift 37 --extsize 73 >& /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV2.log &

# macs2 callpeak -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53957.bam /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53981.bam /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_54005.bam -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53941.bam /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53965.bam /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53989.bam -f BAM -g ce --outdir /data/public/cdebes/workspace/scripts/ssl1_chip_seq/ -n HTZ-SSL1-UV --bdg --nomodel --shift 37 --extsize 73 >& /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV.log &

# #R2
# macs2 callpeak -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53967.bam -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53959.bam -f BAM -g ce --outdir /data/public/cdebes/workspace/scripts/ssl1_chip_seq/ -n RNA-WT1-C2 --bdg --nomodel --shift 37 --extsize 73 >& /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-C2.log &
# macs2 callpeak -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53969.bam -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53961.bam -f BAM -g ce --outdir /data/public/cdebes/workspace/scripts/ssl1_chip_seq/ -n RNA-SSL1-C2 --bdg --nomodel --shift 37 --extsize 73 >& /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-C2.log &
# macs2 callpeak -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53971.bam -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53963.bam -f BAM -g ce --outdir /data/public/cdebes/workspace/scripts/ssl1_chip_seq/ -n RNA-WT1-UV2 --bdg --nomodel --shift 37 --extsize 73 >& /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-UV2.log &
# macs2 callpeak -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53973.bam -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53965.bam -f BAM -g ce --outdir /data/public/cdebes/workspace/scripts/ssl1_chip_seq/ -n RNA-SSL1-UV2 --bdg --nomodel --shift 37 --extsize 73 >& /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-UV2.log &

# macs2 callpeak -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53975.bam -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53959.bam -f BAM -g ce --outdir /data/public/cdebes/workspace/scripts/ssl1_chip_seq/ -n HTZ-WT1-C2 --bdg --nomodel --shift 37 --extsize 73 >& /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-C2.log &
# macs2 callpeak -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53977.bam -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53961.bam -f BAM -g ce --outdir /data/public/cdebes/workspace/scripts/ssl1_chip_seq/ -n HTZ-SSL1-C2 --bdg --nomodel --shift 37 --extsize 73 >& /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-C2.log &
# macs2 callpeak -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53979.bam -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53963.bam -f BAM -g ce --outdir /data/public/cdebes/workspace/scripts/ssl1_chip_seq/ -n HTZ-WT1-UV2 --bdg --nomodel --shift 37 --extsize 73 >& /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV2.log &
# macs2 callpeak -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53981.bam -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53965.bam -f BAM -g ce --outdir /data/public/cdebes/workspace/scripts/ssl1_chip_seq/ -n HTZ-SSL1-UV2 --bdg --nomodel --shift 37 --extsize 73 >& /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV2.log &

# #R3
# macs2 callpeak -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53991.bam -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53983.bam -f BAM -g ce --outdir /data/public/cdebes/workspace/scripts/ssl1_chip_seq/ -n RNA-WT1-C3 --bdg --nomodel --shift 37 --extsize 73 >& /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-C3.log &
# macs2 callpeak -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53993.bam -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53985.bam -f BAM -g ce --outdir /data/public/cdebes/workspace/scripts/ssl1_chip_seq/ -n RNA-SSL1-C3 --bdg --nomodel --shift 37 --extsize 73 >& /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-C3.log &
# macs2 callpeak -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53995.bam -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53987.bam -f BAM -g ce --outdir /data/public/cdebes/workspace/scripts/ssl1_chip_seq/ -n RNA-WT1-UV3 --bdg --nomodel --shift 37 --extsize 73 >& /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-UV3.log &
# macs2 callpeak -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53997.bam -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53989.bam -f BAM -g ce --outdir /data/public/cdebes/workspace/scripts/ssl1_chip_seq/ -n RNA-SSL1-UV3 --bdg --nomodel --shift 37 --extsize 73 >& /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-UV3.log &

# macs2 callpeak -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53999.bam -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53983.bam -f BAM -g ce --outdir /data/public/cdebes/workspace/scripts/ssl1_chip_seq/ -n HTZ-WT1-C3 --bdg --nomodel --shift 37 --extsize 73 >& /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-C3.log &
# macs2 callpeak -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_54001.bam -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53985.bam -f BAM -g ce --outdir /data/public/cdebes/workspace/scripts/ssl1_chip_seq/ -n HTZ-SSL1-C3 --bdg --nomodel --shift 37 --extsize 73 >& /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-C3.log &
# macs2 callpeak -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_54003.bam -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53987.bam -f BAM -g ce --outdir /data/public/cdebes/workspace/scripts/ssl1_chip_seq/ -n HTZ-WT1-UV3 --bdg --nomodel --shift 37 --extsize 73 >& /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV3.log &
# macs2 callpeak -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_54005.bam -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53989.bam -f BAM -g ce --outdir /data/public/cdebes/workspace/scripts/ssl1_chip_seq/ -n HTZ-SSL1-UV3 --bdg --nomodel --shift 37 --extsize 73 >& /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV3.log &


# #bdgcmp HTZ first sample
# macs2 bdgcmp -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT2_treat_pileup.bdg -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT2_control_lambda.bdg -m FE -o /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT2.bdg
# macs2 bdgcmp -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL2_treat_pileup.bdg -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL2_control_lambda.bdg -m FE -o /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL2.bdg
# macs2 bdgcmp -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV2_treat_pileup.bdg -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV2_control_lambda.bdg -m FE -o /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV2.bdg
# macs2 bdgcmp -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV2_treat_pileup.bdg -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV2_control_lambda.bdg -m FE -o /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV2.bdg

# #1
# macs2 bdgcmp -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1_treat_pileup.bdg -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1_control_lambda.bdg -m FE -o /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1.bdg
# macs2 bdgcmp -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1_treat_pileup.bdg -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1_control_lambda.bdg -m FE -o /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1.bdg
# macs2 bdgcmp -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV_treat_pileup.bdg -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV_control_lambda.bdg -m FE -o /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV.bdg
# macs2 bdgcmp -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV_treat_pileup.bdg -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV_control_lambda.bdg -m FE -o /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV.bdg

# macs2 bdgcmp -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-C2_treat_pileup.bdg -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-C2_control_lambda.bdg -m FE -o /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-C2.bdg
# macs2 bdgcmp -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-C2_treat_pileup.bdg -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-C2_control_lambda.bdg -m FE -o /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-C2.bdg
# macs2 bdgcmp -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV2_treat_pileup.bdg -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV2_control_lambda.bdg -m FE -o /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV2.bdg
# macs2 bdgcmp -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV2_treat_pileup.bdg -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV2_control_lambda.bdg -m FE -o /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV2.bdg

# macs2 bdgcmp -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-C3_treat_pileup.bdg -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-C3_control_lambda.bdg -m FE -o /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-C3.bdg
# macs2 bdgcmp -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-C3_treat_pileup.bdg -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-C3_control_lambda.bdg -m FE -o /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-C3.bdg
# macs2 bdgcmp -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV3_treat_pileup.bdg -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV3_control_lambda.bdg -m FE -o /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV3.bdg
# macs2 bdgcmp -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV3_treat_pileup.bdg -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV3_control_lambda.bdg -m FE -o /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV3.bdg

# #bdgcmp RNA
# macs2 bdgcmp -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1_treat_pileup.bdg -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1_control_lambda.bdg -m FE -o /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1.bdg
# macs2 bdgcmp -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1_treat_pileup.bdg -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1_control_lambda.bdg -m FE -o /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1.bdg
# macs2 bdgcmp -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-UV_treat_pileup.bdg -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-UV_control_lambda.bdg -m FE -o /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-UV.bdg
# macs2 bdgcmp -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-UV_treat_pileup.bdg -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-UV_control_lambda.bdg -m FE -o /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-UV.bdg

# macs2 bdgcmp -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-C2_treat_pileup.bdg -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-C2_control_lambda.bdg -m FE -o /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-C2.bdg
# macs2 bdgcmp -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-C2_treat_pileup.bdg -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-C2_control_lambda.bdg -m FE -o /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-C2.bdg
# macs2 bdgcmp -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-UV2_treat_pileup.bdg -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-UV2_control_lambda.bdg -m FE -o /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-UV2.bdg
# macs2 bdgcmp -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-UV2_treat_pileup.bdg -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-UV2_control_lambda.bdg -m FE -o /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-UV2.bdg

# macs2 bdgcmp -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-C3_treat_pileup.bdg -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-C3_control_lambda.bdg -m FE -o /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-C3.bdg
# macs2 bdgcmp -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-C3_treat_pileup.bdg -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-C3_control_lambda.bdg -m FE -o /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-C3.bdg
# macs2 bdgcmp -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-UV3_treat_pileup.bdg -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-UV3_control_lambda.bdg -m FE -o /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-UV3.bdg
# macs2 bdgcmp -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-UV3_treat_pileup.bdg -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-UV3_control_lambda.bdg -m FE -o /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-UV3.bdg


# macs2 bdgcmp -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT2_treat_pileup.bdg -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT2_control_lambda.bdg -m FE -o /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT2.bdg
# macs2 bdgcmp -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL2_treat_pileup.bdg -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL2_control_lambda.bdg -m FE -o /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL2.bdg
# macs2 bdgcmp -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV2_treat_pileup.bdg -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV2_control_lambda.bdg -m FE -o /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV2.bdg
# macs2 bdgcmp -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV2_treat_pileup.bdg -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV2_control_lambda.bdg -m FE -o /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV2.bdg
 
# sort -k1,1 -k2,2n /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT2.bdg > /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT2.sorted.bdg
# ~/Downloads/bedGraphToBigWig /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT2.sorted.bdg /data/public/cdebes/workspace/genomes/ce10.genome /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT2.bw
# sort -k1,1 -k2,2n /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL2.bdg > /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL2.sorted.bdg
# ~/Downloads/bedGraphToBigWig /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL2.sorted.bdg /data/public/cdebes/workspace/genomes/ce10.genome /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL2.bw
# sort -k1,1 -k2,2n /data/public/cdebes/workspace/scripts/ssl1_chip_seq/WT1-UV2.bdg > /data/public/cdebes/workspace/scripts/ssl1_chip_seq/WT1-UV2.sorted.bdg
# ~/Downloads/bedGraphToBigWig /data/public/cdebes/workspace/scripts/ssl1_chip_seq/WT1-UV2.sorted.bdg /data/public/cdebes/workspace/genomes/ce10.genome /data/public/cdebes/workspace/scripts/ssl1_chip_seq/WT1-UV2.bw
# sort -k1,1 -k2,2n /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV2.bdg > /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV2.sorted.bdg
# ~/Downloads/bedGraphToBigWig /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV2.sorted.bdg /data/public/cdebes/workspace/genomes/ce10.genome /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV2.bw
# ~/Downloads/wigCorrelate /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT2.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1.bw


# #wt HTZ
# sort -k1,1 -k2,2n /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1.bdg > /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1.sorted.bdg
# ~/Downloads/bedGraphToBigWig /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1.sorted.bdg /data/public/cdebes/workspace/genomes/ce10.genome /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1.bw
# sort -k1,1 -k2,2n /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-C2.bdg > /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-C2.sorted.bdg
# ~/Downloads/bedGraphToBigWig /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-C2.sorted.bdg /data/public/cdebes/workspace/genomes/ce10.genome /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-C2.bw
# sort -k1,1 -k2,2n /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-C3.bdg > /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-C3.sorted.bdg
# ~/Downloads/bedGraphToBigWig /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-C3.sorted.bdg /data/public/cdebes/workspace/genomes/ce10.genome /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-C3.bw

# ~/Downloads/wigCorrelate /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-C1.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-C2.bw

# #ssl1 HTZ
# sort -k1,1 -k2,2n /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1.bdg > /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1.sorted.bdg
# ~/Downloads/bedGraphToBigWig /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1.sorted.bdg /data/public/cdebes/workspace/genomes/ce10.genome /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1.bw
# sort -k1,1 -k2,2n /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-C2.bdg > /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-C2.sorted.bdg
# ~/Downloads/bedGraphToBigWig /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-C2.sorted.bdg /data/public/cdebes/workspace/genomes/ce10.genome /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-C2.bw
# sort -k1,1 -k2,2n /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-C3.bdg > /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-C3.sorted.bdg
# ~/Downloads/bedGraphToBigWig /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-C3.sorted.bdg /data/public/cdebes/workspace/genomes/ce10.genome /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-C3.bw

# #ssl1uv HTZ
# sort -k1,1 -k2,2n /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV.bdg > /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV.sorted.bdg
# ~/Downloads/bedGraphToBigWig /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV.sorted.bdg /data/public/cdebes/workspace/genomes/ce10.genome /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV.bw
# sort -k1,1 -k2,2n /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV2.bdg > /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV2.sorted.bdg
# ~/Downloads/bedGraphToBigWig /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV2.sorted.bdg /data/public/cdebes/workspace/genomes/ce10.genome /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV2.bw
# sort -k1,1 -k2,2n /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV3.bdg > /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV3.sorted.bdg
# ~/Downloads/bedGraphToBigWig /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV3.sorted.bdg /data/public/cdebes/workspace/genomes/ce10.genome /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV3.bw

# #wt1uv HTZ
# sort -k1,1 -k2,2n /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV.bdg > /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV.sorted.bdg
# ~/Downloads/bedGraphToBigWig /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV.sorted.bdg /data/public/cdebes/workspace/genomes/ce10.genome /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV.bw
# sort -k1,1 -k2,2n /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV2.bdg > /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV2.sorted.bdg
# ~/Downloads/bedGraphToBigWig /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV2.sorted.bdg /data/public/cdebes/workspace/genomes/ce10.genome /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV2.bw
# sort -k1,1 -k2,2n /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV3.bdg > /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV3.sorted.bdg
# ~/Downloads/bedGraphToBigWig /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV3.sorted.bdg /data/public/cdebes/workspace/genomes/ce10.genome /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV3.bw

# #wt RNA
# sort -k1,1 -k2,2n /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1.bdg > /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1.sorted.bdg
# ~/Downloads/bedGraphToBigWig /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1.sorted.bdg /data/public/cdebes/workspace/genomes/ce10.genome /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1.bw
# sort -k1,1 -k2,2n /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-C2.bdg > /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-C2.sorted.bdg
# ~/Downloads/bedGraphToBigWig /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-C2.sorted.bdg /data/public/cdebes/workspace/genomes/ce10.genome /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-C2.bw
# sort -k1,1 -k2,2n /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-C3.bdg > /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-C3.sorted.bdg
# ~/Downloads/bedGraphToBigWig /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-C3.sorted.bdg /data/public/cdebes/workspace/genomes/ce10.genome /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-C3.bw

# ~/Downloads/wigCorrelate /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-C1.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-C2.bw

# #ssl1 RNA
# sort -k1,1 -k2,2n /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1.bdg > /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1.sorted.bdg
# ~/Downloads/bedGraphToBigWig /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1.sorted.bdg /data/public/cdebes/workspace/genomes/ce10.genome /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1.bw
# sort -k1,1 -k2,2n /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-C2.bdg > /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-C2.sorted.bdg
# ~/Downloads/bedGraphToBigWig /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-C2.sorted.bdg /data/public/cdebes/workspace/genomes/ce10.genome /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-C2.bw
# sort -k1,1 -k2,2n /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-C3.bdg > /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-C3.sorted.bdg
# ~/Downloads/bedGraphToBigWig /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-C3.sorted.bdg /data/public/cdebes/workspace/genomes/ce10.genome /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-C3.bw

# #ssl1uv RNA
# sort -k1,1 -k2,2n /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-UV.bdg > /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-UV.sorted.bdg
# ~/Downloads/bedGraphToBigWig /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-UV.sorted.bdg /data/public/cdebes/workspace/genomes/ce10.genome /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-UV.bw
# sort -k1,1 -k2,2n /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-UV2.bdg > /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-UV2.sorted.bdg
# ~/Downloads/bedGraphToBigWig /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-UV2.sorted.bdg /data/public/cdebes/workspace/genomes/ce10.genome /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-UV2.bw
# sort -k1,1 -k2,2n /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-UV3.bdg > /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-UV3.sorted.bdg
# ~/Downloads/bedGraphToBigWig /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-UV3.sorted.bdg /data/public/cdebes/workspace/genomes/ce10.genome /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-UV3.bw

# #wt1uv RNA
# sort -k1,1 -k2,2n /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-UV.bdg > /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-UV.sorted.bdg
# ~/Downloads/bedGraphToBigWig /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-UV.sorted.bdg /data/public/cdebes/workspace/genomes/ce10.genome /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-UV.bw
# sort -k1,1 -k2,2n /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-UV2.bdg > /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-UV2.sorted.bdg
# ~/Downloads/bedGraphToBigWig /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-UV2.sorted.bdg /data/public/cdebes/workspace/genomes/ce10.genome /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-UV2.bw
# sort -k1,1 -k2,2n /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-UV3.bdg > /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-UV3.sorted.bdg
# ~/Downloads/bedGraphToBigWig /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-UV3.sorted.bdg /data/public/cdebes/workspace/genomes/ce10.genome /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-UV3.bw



# multiBigwigSummary bins -b /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-C1.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-C2.bw /data/public/cdebes//workspace/scripts/ssl1_chip_seq/RNA-WT1-C3.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-C1.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-C2.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-C3.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-UV1.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-UV2.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-UV3.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-UV1.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-UV2.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-UV3.bw -out RNA.npz&
# plotCorrelation -in RNA.npz -o RNA.png -c pearson -p heatmap --skipZeros --removeOutliers

# multiBigwigSummary bins -b /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-C1.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-C2.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-C3.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-C1.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-C2.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-C3.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV1.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV2.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV3.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV1.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV2.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV3.bw -out HTZ.npz&
# plotCorrelation -in HTZ.npz -o HTZ.png -c pearson -p heatmap --skipZeros --removeOutliers

# #5prim
# computeMatrix reference-point --referencePoint center -b 500 -a 500 -R /data/public/cdebes/workspace/genomes/celegans/5prim_ensembl_flank.bed /data/public/cdebes/workspace/genomes/celegans/intron_junction_ensembl_coverage_exon_intron_5norepeat.bed /data/public/cdebes/workspace/genomes/celegans/3prim_ensembl_flank.bed -S /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1.bw  /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV.bw  -o matrixHTZ5prim.gz&
# plotProfile --matrixFile matrixHTZ5prim.gz -out  profileHTZ5prim.png --perGroup

# computeMatrix reference-point --referencePoint center -b 500 -a 500 -R /data/public/cdebes/workspace/genomes/celegans/5prim_ensembl_flank.bed /data/public/cdebes/workspace/genomes/celegans/intron_junction_ensembl_coverage_exon_intron_5norepeat.bed /data/public/cdebes/workspace/genomes/celegans/3prim_ensembl_flank.bed -S /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1.bw  /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-UV.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1.bw  /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-UV.bw  -o matrixRNA5prim.gz&
# plotProfile --matrixFile matrixRNA5prim.gz -out  profileRNA5prim.png --perGroup

# #geneRNA
# computeMatrix scale-regions -b 200 -a 200 -R /data/public/cdebes/workspace/genomes/celegans/gene_ensembl.bed /data/public/cdebes/workspace/scripts/ssl1/NoRecIntronJunction.bed -S /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1.bw  /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-UV.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-UV.bw  -o matrixRNA_scale_gene.gz&
# plotProfile --matrixFile matrixRNA_scale_gene.gz -out  profileRNA_scale_gene.png --perGroup

# #geneHTZ
# computeMatrix scale-regions -b 200 -a 200 -R /data/public/cdebes/workspace/genomes/celegans/gene_ensembl.bed /data/public/cdebes/workspace/scripts/ssl1/NoRecIntronJunction.bed -S /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV.bw  -o matrixHTZ_scale_gene.gz&
# plotProfile --matrixFile matrixHTZ_scale_gene.gz -out  profileHTZ_scale_gene.png --perGroup
