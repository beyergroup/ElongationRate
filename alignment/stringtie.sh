#!/bin/bash -l

#SBATCH --account=cdebes
#SBATCH --error=/data/public/cdebes/stringtie.log
#SBATCH --job-name=stringtiehs
#SBATCH --partition=all
#SBATCH --ntasks=12

#rnorvegicus
list=(`find /data/public/cdebes/workspace/scripts/rnorvegicus/ -name *.sortedByCoord.out.bam`)

for i in ${list[*]}
do
    bname=`basename $i`
    /data/public/cdebes/soft/stringtie-1.3.3b.Linux_x86_64/stringtie $i -a 1 -c 1 -f 0.01 -G /data/public/cdebes/workspace/genomes/rnorvegicus/Rattus_norvegicus.Rnor_6.0.85.gtf  -p 12 -j 0 -a 5 -o /data/public/cdebes/workspace/scripts/rnorvegicus/${bname:0:10}.gtf -A /data/public/cdebes/workspace/scripts/rnorvegicus/${bname:0:10}.abundance.tab
done

/data/public/cdebes/soft/stringtie-1.3.3b.Linux_x86_64/stringtie --merge -G /data/public/cdebes/workspace/genomes/rnorvegicus/Rattus_norvegicus.Rnor_6.0.85.gtf /data/public/cdebes/workspace/scripts/rnorvegicus/SRR1873526.gtf /data/public/cdebes/workspace/scripts/rnorvegicus/SRR1873527.gtf /data/public/cdebes/workspace/scripts/rnorvegicus/SRR1873528.gtf /data/public/cdebes/workspace/scripts/rnorvegicus/SRR1873529.gtf /data/public/cdebes/workspace/scripts/rnorvegicus/SRR1873530.gtf /data/public/cdebes/workspace/scripts/rnorvegicus/SRR1873531.gtf -o /data/public/cdebes/workspace/scripts/rnorvegicus/combined_stringtie.gtf

for i in ${list[*]}
do
    bname=`basename $i`
    /data/public/cdebes/soft/stringtie-1.3.3b.Linux_x86_64/stringtie $i -v -e -G /data/public/cdebes/workspace/scripts/rnorvegicus/combined_stringtie.gtf -p 12 -o /data/public/cdebes/workspace/scripts/rnorvegicus/${bname:0:10}_filter1.gtf
done

#hsapiens
list=(`find /data/public/cdebes/workspace/scripts/hsapiens/ -name SN*.sortedByCoord.out.bam`)

for i in ${list[*]}
do
    bname=`basename $i`
    /data/public/cdebes/soft/stringtie-1.3.3b.Linux_x86_64/stringtie $i -a 1 -c 1 -f 0.01  -v -G /data/public/cdebes/workspace/genomes/hsapiens/Homo_sapiens.GRCh38.85.gtf -p 12 -o /data/public/cdebes/workspace/scripts/hsapiens/${bname:0:22}.gtf -A /data/public/cdebes/workspace/scripts/hsapiens/${bname:0:22}.abundance.tab
done

/data/public/cdebes/soft/stringtie-1.3.3b.Linux_x86_64/stringtie --merge -G /data/public/cdebes/workspace/genomes/hsapiens/Homo_sapiens.GRCh38.85.gtf /data/public/cdebes/workspace/scripts/hsapiens/SN7640229_17803_D1youn.gtf /data/public/cdebes/workspace/scripts/hsapiens/SN5370216_20212_D2youn.gtf /data/public/cdebes/workspace/scripts/hsapiens/SN5370216_20210_D3youn.gtf /data/public/cdebes/workspace/scripts/hsapiens/SN7640229_17804_D1oldA.gtf /data/public/cdebes/workspace/scripts/hsapiens/SN5370216_20213_D2oldA.gtf /data/public/cdebes/workspace/scripts/hsapiens/SN5370216_20211_D3oldA.gtf -o /data/public/cdebes/workspace/scripts/hsapiens/combined_stringtie.gtf

for i in ${list[*]}
do
    bname=`basename $i`
    /data/public/cdebes/soft/stringtie-1.3.3b.Linux_x86_64/stringtie $i -v -e -G /data/public/cdebes/workspace/scripts/hsapiens/combined_stringtie.gtf -p 12 -o /data/public/cdebes/workspace/scripts/hsapiens/${bname:0:22}.gtf -A /data/public/cdebes/workspace/scripts/hsapiens/${bname:0:22}.abundance.tab
done


#hsapiens_nascent
list=(`find /data/public/cdebes/workspace/scripts/hsapiens_nascent/ -name SN*.sortedByCoord.out.bam`)

for i in ${list[*]}
do
    bname=`basename $i`
    /data/public/cdebes/soft/stringtie-1.3.3b.Linux_x86_64/stringtie $i -a 1 -c 1 -f 0.01 -v -G /data/public/cdebes/workspace/genomes/hsapiens/Homo_sapiens.GRCh38.85.gtf -p 12 -o /data/public/cdebes/workspace/scripts/hsapiens/${bname:0:22}.gtf -A /data/public/cdebes/workspace/scripts/hsapiens/${bname:0:22}.abundance.tab
done

/data/public/cdebes/soft/stringtie-1.3.3b.Linux_x86_64/stringtie --merge -G /data/public/cdebes/workspace/genomes/hsapiens/Homo_sapiens.GRCh38.85.gtf /data/public/cdebes/workspace/scripts/hsapiens_nascent/SN7640329_38948_OldR1A.gtf /data/public/cdebes/workspace/scripts/hsapiens_nascent/SN7640329_38958_I79you.gtf /data/public/cdebes/workspace/scripts/hsapiens_nascent/SN7640329_38960_I10you.gtf /data/public/cdebes/workspace/scripts/hsapiens_nascent/SN7640329_38961_I10old.gtf /data/public/cdebes/workspace/scripts/hsapiens_nascent/SN7640329_38959_I79old.gtf /data/public/cdebes/workspace/scripts/hsapiens_nascent/SN7640329_38949_OldR2A.gtf /data/public/cdebes/workspace/scripts/hsapiens_nascent/SN7640329_38947_YoungR.gtf /data/public/cdebes/workspace/scripts/hsapiens_nascent/SN7640329_38946_YoungR.gtf -o /data/public/cdebes/workspace/scripts/hsapiens_nascent/combined_stringtie.gtf

for i in ${list[*]}
do
    bname=`basename $i`
    /data/public/cdebes/soft/stringtie-1.3.3b.Linux_x86_64/stringtie $i -v -e -G /data/public/cdebes/workspace/scripts/hsapiens_nascent/combined_stringtie.gtf -p 12 -o /data/public/cdebes/workspace/scripts/hsapiens_nascent/${bname:0:22}.gtf -A /data/public/cdebes/workspace/scripts/hsapiens_nascent/${bname:0:22}.abundance.tab
done

list=(`find /data/public/cdebes/workspace/scripts/hsapiens_blood/ -name *.sortedByCoord.out.bam`)
for i in ${list[*]}
do
    bname=`basename $i`
    /data/public/cdebes/soft/stringtie-1.3.3b.Linux_x86_64/stringtie $i -a 1 -c 1 -f 0.01  -v -G /data/public/cdebes/workspace/genomes/hsapiens/Homo_sapiens.GRCh38.85.gtf -p 12 -o /data/public/cdebes/workspace/scripts/hsapiens_blood/${bname:0:22}.gtf -A /data/public/cdebes/workspace/scripts/hsapiens_blood/${bname:0:22}.abundance.tab
done


