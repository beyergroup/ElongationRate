macs2 callpeak -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53943_bwa.bam /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53967_bwa.bam /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53991_bwa.bam -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53935_bwa.bam /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53959_bwa.bam /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53983_bwa.bam -f BAM -g ce --outdir /data/public/cdebes/workspace/scripts/ssl1_chip_seq/ -n RNA-WT1 --bdg  >& /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1.log &
macs2 callpeak -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53945_bwa.bam /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53969_bwa.bam /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53993_bwa.bam -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53937_bwa.bam /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53961_bwa.bam /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53985_bwa.bam -f BAM -g ce --outdir /data/public/cdebes/workspace/scripts/ssl1_chip_seq/ -n RNA-SSL1 --bdg --nomodel  >& /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1.log &
macs2 callpeak -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53947_bwa.bam /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53971_bwa.bam /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53995_bwa.bam -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53939_bwa.bam /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53963_bwa.bam /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53987_bwa.bam -f BAM -g ce --outdir /data/public/cdebes/workspace/scripts/ssl1_chip_seq/ -n RNA-WT1-UV --bdg --nomodel  >& /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-UV.log &
macs2 callpeak -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53949_bwa.bam /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53973_bwa.bam /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53997_bwa.bam -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53941_bwa.bam /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53965_bwa.bam /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53989_bwa.bam -f BAM -g ce --outdir /data/public/cdebes/workspace/scripts/ssl1_chip_seq/ -n RNA-SSL1-UV --bdg --nomodel  >& /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-UV.log &

macs2 callpeak -t /cellnet/SSL1/data/chip-seq/48_N2HTZC.bam -c /cellnet/SSL1/data/chip-seq/56_N2InputC.bam -f BAM -g ce --outdir /data/public/cdebes/workspace/scripts/ssl1_chip_seq/ -n HTZ-WT2 --bdg --nomodel --shift 37 --extsize 73 >& /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT2.log &

macs2 callpeak -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53951_bwa.bam /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53975_bwa.bam /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53999_bwa.bam -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53935_bwa.bam /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53959_bwa.bam /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53983_bwa.bam -f BAM -g ce --outdir /data/public/cdebes/workspace/scripts/ssl1_chip_seq/ -n HTZ-WT1 --bdg --nomodel --shift 37 --extsize 73 >& /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1.log &

macs2 callpeak -t /cellnet/SSL1/data/chip-seq/50_ssl1HTZC.bam -c /cellnet/SSL1/data/chip-seq/58_ssl1InputC.bam -f BAM -g ce --outdir /data/public/cdebes/workspace/scripts/ssl1_chip_seq/ -n HTZ-SSL2 --bdg --nomodel --shift 37 --extsize 73 >& /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL2.log &

macs2 callpeak -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53953_bwa.bam /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53977_bwa.bam /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_54001_bwa.bam -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53937_bwa.bam /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53961_bwa.bam /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53985_bwa.bam -f BAM -g ce --outdir /data/public/cdebes/workspace/scripts/ssl1_chip_seq/ -n HTZ-SSL1 --bdg --nomodel --shift 37 --extsize 73 >& /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1.log &

macs2 callpeak -t /cellnet/SSL1/data/chip-seq/49_N2HTZUV.bam  -c /cellnet/SSL1/data/chip-seq/57_N2InputUV.bam  -f BAM -g ce --outdir /data/public/cdebes/workspace/scripts/ssl1_chip_seq/ -n HTZ-WT1-UV2 --bdg --nomodel --shift 37 --extsize 73 >& /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV2.log &

macs2 callpeak -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53955_bwa.bam /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53979_bwa.bam /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_54003_bwa.bam -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53939_bwa.bam /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53963_bwa.bam /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53987_bwa.bam -f BAM -g ce --outdir /data/public/cdebes/workspace/scripts/ssl1_chip_seq/ -n HTZ-WT1-UV --bdg --nomodel --shift 37 --extsize 73 >& /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV.log &

macs2 callpeak -t /cellnet/SSL1/data/chip-seq/51_ssl1HTZUV.bam   -c /cellnet/SSL1/data/chip-seq/59_ssl1InputUV.bam  -f BAM -g ce --outdir /data/public/cdebes/workspace/scripts/ssl1_chip_seq/ -n HTZ-SSL1-UV2 --bdg --nomodel --shift 37 --extsize 73 >& /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV2.log &

macs2 callpeak -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53957_bwa.bam /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53981_bwa.bam /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_54005_bwa.bam -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53941_bwa.bam /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53965_bwa.bam /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53989_bwa.bam -f BAM -g ce --outdir /data/public/cdebes/workspace/scripts/ssl1_chip_seq/ -n HTZ-SSL1-UV --bdg --nomodel --shift 37 --extsize 73 >& /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV.log &

# #R2
# macs2 callpeak -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53967_bwa.bam -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53959_bwa.bam -f BAM -g ce --outdir /data/public/cdebes/workspace/scripts/ssl1_chip_seq/ -n RNA-WT1-C2 --bdg --nomodel --shift 37 --extsize 73 >& /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-C2.log &
# macs2 callpeak -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53969_bwa.bam -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53961_bwa.bam -f BAM -g ce --outdir /data/public/cdebes/workspace/scripts/ssl1_chip_seq/ -n RNA-SSL1-C2 --bdg --nomodel --shift 37 --extsize 73 >& /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-C2.log &
# macs2 callpeak -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53971_bwa.bam -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53963_bwa.bam -f BAM -g ce --outdir /data/public/cdebes/workspace/scripts/ssl1_chip_seq/ -n RNA-WT1-UV2 --bdg --nomodel --shift 37 --extsize 73 >& /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-UV2.log &
# macs2 callpeak -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53973_bwa.bam -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53965_bwa.bam -f BAM -g ce --outdir /data/public/cdebes/workspace/scripts/ssl1_chip_seq/ -n RNA-SSL1-UV2 --bdg --nomodel --shift 37 --extsize 73 >& /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-UV2.log &

# macs2 callpeak -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53975_bwa.bam -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53959_bwa.bam -f BAM -g ce --outdir /data/public/cdebes/workspace/scripts/ssl1_chip_seq/ -n HTZ-WT1-C2 --bdg --nomodel --shift 37 --extsize 73 >& /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-C2.log &
# macs2 callpeak -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53977_bwa.bam -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53961_bwa.bam -f BAM -g ce --outdir /data/public/cdebes/workspace/scripts/ssl1_chip_seq/ -n HTZ-SSL1-C2 --bdg --nomodel --shift 37 --extsize 73 >& /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-C2.log &
# macs2 callpeak -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53979_bwa.bam -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53963_bwa.bam -f BAM -g ce --outdir /data/public/cdebes/workspace/scripts/ssl1_chip_seq/ -n HTZ-WT1-UV2 --bdg --nomodel --shift 37 --extsize 73 >& /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV2.log &
# macs2 callpeak -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53981_bwa.bam -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53965_bwa.bam -f BAM -g ce --outdir /data/public/cdebes/workspace/scripts/ssl1_chip_seq/ -n HTZ-SSL1-UV2 --bdg --nomodel --shift 37 --extsize 73 >& /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV2.log &

# #R3
# macs2 callpeak -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53991_bwa.bam -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53983_bwa.bam -f BAM -g ce --outdir /data/public/cdebes/workspace/scripts/ssl1_chip_seq/ -n RNA-WT1-C3 --bdg --nomodel --shift 37 --extsize 73 >& /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-C3.log &
# macs2 callpeak -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53993_bwa.bam -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53985_bwa.bam -f BAM -g ce --outdir /data/public/cdebes/workspace/scripts/ssl1_chip_seq/ -n RNA-SSL1-C3 --bdg --nomodel --shift 37 --extsize 73 >& /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-C3.log &
# macs2 callpeak -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53995_bwa.bam -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53987_bwa.bam -f BAM -g ce --outdir /data/public/cdebes/workspace/scripts/ssl1_chip_seq/ -n RNA-WT1-UV3 --bdg --nomodel --shift 37 --extsize 73 >& /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-UV3.log &
# macs2 callpeak -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53997_bwa.bam -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53989_bwa.bam -f BAM -g ce --outdir /data/public/cdebes/workspace/scripts/ssl1_chip_seq/ -n RNA-SSL1-UV3 --bdg --nomodel --shift 37 --extsize 73 >& /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-UV3.log &

# macs2 callpeak -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53999_bwa.bam -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53983_bwa.bam -f BAM -g ce --outdir /data/public/cdebes/workspace/scripts/ssl1_chip_seq/ -n HTZ-WT1-C3 --bdg --nomodel --shift 37 --extsize 73 >& /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-C3.log &
# macs2 callpeak -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_54001_bwa.bam -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53985_bwa.bam -f BAM -g ce --outdir /data/public/cdebes/workspace/scripts/ssl1_chip_seq/ -n HTZ-SSL1-C3 --bdg --nomodel --shift 37 --extsize 73 >& /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-C3.log &
# macs2 callpeak -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_54003_bwa.bam -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53987_bwa.bam -f BAM -g ce --outdir /data/public/cdebes/workspace/scripts/ssl1_chip_seq/ -n HTZ-WT1-UV3 --bdg --nomodel --shift 37 --extsize 73 >& /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV3.log &
# macs2 callpeak -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_54005_bwa.bam -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/K002000093_53989_bwa.bam -f BAM -g ce --outdir /data/public/cdebes/workspace/scripts/ssl1_chip_seq/ -n HTZ-SSL1-UV3 --bdg --nomodel --shift 37 --extsize 73 >& /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV3.log &

#bdgcmp HTZ first sample
macs2 bdgcmp -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT2_treat_pileup.bdg -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT2_control_lambda.bdg -m FE -o /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT2.bdg
macs2 bdgcmp -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL2_treat_pileup.bdg -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL2_control_lambda.bdg -m FE -o /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL2.bdg
macs2 bdgcmp -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV2_treat_pileup.bdg -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV2_control_lambda.bdg -m FE -o /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV2.bdg
macs2 bdgcmp -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV2_treat_pileup.bdg -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV2_control_lambda.bdg -m FE -o /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV2.bdg

#1
macs2 bdgcmp -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1_treat_pileup.bdg -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1_control_lambda.bdg -m FE -o /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1.bdg
macs2 bdgcmp -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1_treat_pileup.bdg -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1_control_lambda.bdg -m FE -o /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1.bdg
macs2 bdgcmp -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV_treat_pileup.bdg -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV_control_lambda.bdg -m FE -o /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV.bdg
macs2 bdgcmp -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV_treat_pileup.bdg -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV_control_lambda.bdg -m FE -o /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV.bdg

macs2 bdgcmp -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-C2_treat_pileup.bdg -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-C2_control_lambda.bdg -m FE -o /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-C2.bdg
macs2 bdgcmp -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-C2_treat_pileup.bdg -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-C2_control_lambda.bdg -m FE -o /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-C2.bdg
macs2 bdgcmp -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV2_treat_pileup.bdg -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV2_control_lambda.bdg -m FE -o /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV2.bdg
macs2 bdgcmp -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV2_treat_pileup.bdg -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV2_control_lambda.bdg -m FE -o /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV2.bdg

macs2 bdgcmp -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-C3_treat_pileup.bdg -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-C3_control_lambda.bdg -m FE -o /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-C3.bdg
macs2 bdgcmp -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-C3_treat_pileup.bdg -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-C3_control_lambda.bdg -m FE -o /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-C3.bdg
macs2 bdgcmp -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV3_treat_pileup.bdg -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV3_control_lambda.bdg -m FE -o /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV3.bdg
macs2 bdgcmp -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV3_treat_pileup.bdg -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV3_control_lambda.bdg -m FE -o /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV3.bdg

#bdgcmp RNA
macs2 bdgcmp -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1_treat_pileup.bdg -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1_control_lambda.bdg -m FE -o /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1.bdg
macs2 bdgcmp -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1_treat_pileup.bdg -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1_control_lambda.bdg -m FE -o /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1.bdg
macs2 bdgcmp -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-UV_treat_pileup.bdg -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-UV_control_lambda.bdg -m FE -o /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-UV.bdg
macs2 bdgcmp -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-UV_treat_pileup.bdg -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-UV_control_lambda.bdg -m FE -o /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-UV.bdg

macs2 bdgcmp -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-C2_treat_pileup.bdg -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-C2_control_lambda.bdg -m FE -o /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-C2.bdg
macs2 bdgcmp -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-C2_treat_pileup.bdg -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-C2_control_lambda.bdg -m FE -o /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-C2.bdg
macs2 bdgcmp -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-UV2_treat_pileup.bdg -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-UV2_control_lambda.bdg -m FE -o /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-UV2.bdg
macs2 bdgcmp -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-UV2_treat_pileup.bdg -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-UV2_control_lambda.bdg -m FE -o /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-UV2.bdg

macs2 bdgcmp -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-C3_treat_pileup.bdg -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-C3_control_lambda.bdg -m FE -o /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-C3.bdg
macs2 bdgcmp -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-C3_treat_pileup.bdg -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-C3_control_lambda.bdg -m FE -o /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-C3.bdg
macs2 bdgcmp -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-UV3_treat_pileup.bdg -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-UV3_control_lambda.bdg -m FE -o /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-UV3.bdg
macs2 bdgcmp -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-UV3_treat_pileup.bdg -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-UV3_control_lambda.bdg -m FE -o /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-UV3.bdg

macs2 bdgcmp -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT2_treat_pileup.bdg -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT2_control_lambda.bdg -m FE -o /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT2.bdg
macs2 bdgcmp -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL2_treat_pileup.bdg -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL2_control_lambda.bdg -m FE -o /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL2.bdg
macs2 bdgcmp -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV2_treat_pileup.bdg -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV2_control_lambda.bdg -m FE -o /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV2.bdg
macs2 bdgcmp -t /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV2_treat_pileup.bdg -c /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV2_control_lambda.bdg -m FE -o /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV2.bdg
 
sort -k1,1 -k2,2n /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT2.bdg > /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT2.sorted.bdg
~/Downloads/bedGraphToBigWig /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT2.sorted.bdg /data/public/cdebes/workspace/genomes/ce10.genome /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT2.bw
sort -k1,1 -k2,2n /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL2.bdg > /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL2.sorted.bdg
~/Downloads/bedGraphToBigWig /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL2.sorted.bdg /data/public/cdebes/workspace/genomes/ce10.genome /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL2.bw
sort -k1,1 -k2,2n /data/public/cdebes/workspace/scripts/ssl1_chip_seq/WT1-UV2.bdg > /data/public/cdebes/workspace/scripts/ssl1_chip_seq/WT1-UV2.sorted.bdg
~/Downloads/bedGraphToBigWig /data/public/cdebes/workspace/scripts/ssl1_chip_seq/WT1-UV2.sorted.bdg /data/public/cdebes/workspace/genomes/ce10.genome /data/public/cdebes/workspace/scripts/ssl1_chip_seq/WT1-UV2.bw
sort -k1,1 -k2,2n /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV2.bdg > /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV2.sorted.bdg
~/Downloads/bedGraphToBigWig /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV2.sorted.bdg /data/public/cdebes/workspace/genomes/ce10.genome /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV2.bw

~/Downloads/wigCorrelate /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT2.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1.bw
~/Downloads/wigCorrelate /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-UV2.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-UV1.bw
#wt HTZ
sort -k1,1 -k2,2n /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1.bdg > /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1.sorted.bdg
~/Downloads/bedGraphToBigWig /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1.sorted.bdg /data/public/cdebes/workspace/genomes/ce10.genome /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1.bw
sort -k1,1 -k2,2n /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-C2.bdg > /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-C2.sorted.bdg
~/Downloads/bedGraphToBigWig /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-C2.sorted.bdg /data/public/cdebes/workspace/genomes/ce10.genome /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-C2.bw
sort -k1,1 -k2,2n /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-C3.bdg > /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-C3.sorted.bdg
~/Downloads/bedGraphToBigWig /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-C3.sorted.bdg /data/public/cdebes/workspace/genomes/ce10.genome /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-C3.bw

~/Downloads/wigCorrelate /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1.bw

#ssl1 HTZ
sort -k1,1 -k2,2n /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1.bdg > /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1.sorted.bdg
~/Downloads/bedGraphToBigWig /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1.sorted.bdg /data/public/cdebes/workspace/genomes/ce10.genome /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1.bw
sort -k1,1 -k2,2n /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-C2.bdg > /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-C2.sorted.bdg
~/Downloads/bedGraphToBigWig /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-C2.sorted.bdg /data/public/cdebes/workspace/genomes/ce10.genome /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-C2.bw
sort -k1,1 -k2,2n /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-C3.bdg > /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-C3.sorted.bdg
~/Downloads/bedGraphToBigWig /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-C3.sorted.bdg /data/public/cdebes/workspace/genomes/ce10.genome /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-C3.bw

#ssl1uv HTZ
sort -k1,1 -k2,2n /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV.bdg > /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV.sorted.bdg
~/Downloads/bedGraphToBigWig /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV.sorted.bdg /data/public/cdebes/workspace/genomes/ce10.genome /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV.bw
sort -k1,1 -k2,2n /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV2.bdg > /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV2.sorted.bdg
~/Downloads/bedGraphToBigWig /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV2.sorted.bdg /data/public/cdebes/workspace/genomes/ce10.genome /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV2.bw
sort -k1,1 -k2,2n /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV3.bdg > /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV3.sorted.bdg
~/Downloads/bedGraphToBigWig /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV3.sorted.bdg /data/public/cdebes/workspace/genomes/ce10.genome /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV3.bw

#wt1uv HTZ
sort -k1,1 -k2,2n /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV.bdg > /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV.sorted.bdg
~/Downloads/bedGraphToBigWig /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV.sorted.bdg /data/public/cdebes/workspace/genomes/ce10.genome /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV.bw
sort -k1,1 -k2,2n /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV2.bdg > /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV2.sorted.bdg
~/Downloads/bedGraphToBigWig /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV2.sorted.bdg /data/public/cdebes/workspace/genomes/ce10.genome /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV2.bw
sort -k1,1 -k2,2n /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV3.bdg > /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV3.sorted.bdg
~/Downloads/bedGraphToBigWig /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV3.sorted.bdg /data/public/cdebes/workspace/genomes/ce10.genome /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV3.bw

#wt RNA
sort -k1,1 -k2,2n /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1.bdg > /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1.sorted.bdg
~/Downloads/bedGraphToBigWig /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1.sorted.bdg /data/public/cdebes/workspace/genomes/ce10.genome /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1.bw
#sort -k1,1 -k2,2n /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-C2.bdg > /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-C2.sorted.bdg
#~/Downloads/bedGraphToBigWig /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-C2.sorted.bdg /data/public/cdebes/workspace/genomes/ce10.genome /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-C2.bw
#sort -k1,1 -k2,2n /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-C3.bdg > /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-C3.sorted.bdg
#~/Downloads/bedGraphToBigWig /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-C3.sorted.bdg /data/public/cdebes/workspace/genomes/ce10.genome /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-C3.bw

#ssl1 RNA
sort -k1,1 -k2,2n /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1.bdg > /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1.sorted.bdg
~/Downloads/bedGraphToBigWig /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1.sorted.bdg /data/public/cdebes/workspace/genomes/ce10.genome /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1.bw
#sort -k1,1 -k2,2n /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-C2.bdg > /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-C2.sorted.bdg
#~/Downloads/bedGraphToBigWig /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-C2.sorted.bdg /data/public/cdebes/workspace/genomes/ce10.genome /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-C2.bw
#sort -k1,1 -k2,2n /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-C3.bdg > /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-C3.sorted.bdg
#~/Downloads/bedGraphToBigWig /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-C3.sorted.bdg /data/public/cdebes/workspace/genomes/ce10.genome /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-C3.bw

#ssl1uv RNA
sort -k1,1 -k2,2n /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-UV.bdg > /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-UV.sorted.bdg
~/Downloads/bedGraphToBigWig /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-UV.sorted.bdg /data/public/cdebes/workspace/genomes/ce10.genome /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-UV.bw
#sort -k1,1 -k2,2n /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-UV2.bdg > /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-UV2.sorted.bdg
#~/Downloads/bedGraphToBigWig /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-UV2.sorted.bdg /data/public/cdebes/workspace/genomes/ce10.genome /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-UV2.bw
#sort -k1,1 -k2,2n /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-UV3.bdg > /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-UV3.sorted.bdg
#~/Downloads/bedGraphToBigWig /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-UV3.sorted.bdg /data/public/cdebes/workspace/genomes/ce10.genome /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-UV3.bw

#wt1uv RNA
sort -k1,1 -k2,2n /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-UV.bdg > /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-UV.sorted.bdg
~/Downloads/bedGraphToBigWig /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-UV.sorted.bdg /data/public/cdebes/workspace/genomes/ce10.genome /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-UV.bw
#sort -k1,1 -k2,2n /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-UV2.bdg > /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-UV2.sorted.bdg
#~/Downloads/bedGraphToBigWig /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-UV2.sorted.bdg /data/public/cdebes/workspace/genomes/ce10.genome /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-UV2.bw
#sort -k1,1 -k2,2n /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-UV3.bdg > /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-UV3.sorted.bdg
#~/Downloads/bedGraphToBigWig /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-UV3.sorted.bdg /data/public/cdebes/workspace/genomes/ce10.genome /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-UV3.bw

multiBigwigSummary bins -b /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-C1.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-C2.bw /data/public/cdebes//workspace/scripts/ssl1_chip_seq/RNA-WT1-C3.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-C1.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-C2.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-C3.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-UV1.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-UV2.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-UV3.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-UV1.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-UV2.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-UV3.bw -out RNA.npz&
plotCorrelation -in RNA.npz -o RNA.png -c pearson -p heatmap --skipZeros --removeOutliers

multiBigwigSummary bins -b /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV2.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/WT1-UV2.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT2.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL2.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-C1.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-C2.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-C3.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-C1.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-C2.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-C3.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV1.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV2.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV3.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV1.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV2.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV3.bw -out HTZ.npz&
plotCorrelation -in HTZ.npz -o HTZ.png -c pearson -p heatmap --skipZeros --removeOutliers

#5prim
computeMatrix reference-point --referencePoint center -b 500 -a 500 -R /data/public/cdebes/workspace/genomes/celegans/5prim_ensembl_flank.bed  /data/public/cdebes/workspace/genomes/celegans/3prim_ensembl_flank.bed -S /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1.bw  /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT2.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV2.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL2.bw  /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV2.bw -o matrixHTZ5prim.gz&
plotProfile --matrixFile matrixHTZ5prim.gz -out  profileHTZ5prim.png --perGroup

computeMatrix reference-point --referencePoint center -b 500 -a 500 -R /data/public/cdebes/workspace/genomes/celegans/5prim_ensembl_flank.bed /data/public/cdebes/workspace/genomes/celegans/3prim_ensembl_flank.bed -S /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1.bw  /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-UV.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1.bw  /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-UV.bw  -o matrixRNA5prim.gz&
plotProfile --matrixFile matrixRNA5prim.gz -out  profileRNA5prim.png --perGroup

computeMatrix reference-point --referencePoint center -b 500 -a 500 -R /data/public/cdebes/workspace/genomes/celegans/intron_junction_ensembl_coverage_exon_intron_5norepeat.bed -S /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1.bw  /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV.bw  -o matrixHTZ5exonintron.gz&
plotProfile --matrixFile matrixHTZ5exonintron.gz -out  profileHTZ5exonintron.png --perGroup

computeMatrix reference-point --referencePoint center -b 500 -a 500 -R /data/public/cdebes/workspace/genomes/celegans/intron_junction_ensembl_coverage_exon_intron_5norepeat.bed -S /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1.bw  /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-UV.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1.bw  /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-UV.bw  -o matrixRNA5exonintron.gz&
plotProfile --matrixFile matrixRNA5exonintron.gz -out  profileRNA5exonintron.png --perGroup

#geneRNA
computeMatrix scale-regions -b 200 -a 200 -R /data/public/cdebes/workspace/genomes/celegans/gene_ensembl.bed -S /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1.bw  /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-UV.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-UV.bw  -o matrixRNA_scale_gene.gz&
plotProfile --matrixFile matrixRNA_scale_gene.gz -out  profileRNA_scale_gene.png --perGroup

computeMatrix scale-regions --regionBodyLength 2000 -R /data/public/cdebes/workspace/scripts/ssl1/NoRecIntronJunction.bed -S /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1.bw  /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-WT1-UV.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/RNA-SSL1-UV.bw  -o matrixRNA_scale_intron.gz&
plotProfile --matrixFile matrixRNA_scale_intron.gz -out  profileRNA_scale_intron.png --perGroup

#geneHTZ
computeMatrix scale-regions -b 200 -a 200 -R /data/public/cdebes/workspace/genomes/celegans/gene_ensembl.bed -S /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV.bw  -o matrixHTZ_scale_gene.gz&
plotProfile --matrixFile matrixHTZ_scale_gene.gz -out  profileHTZ_scale_gene.png --perGroup

computeMatrix scale-regions --regionBodyLength 2000 -R /data/public/cdebes/workspace/scripts/ssl1/NoRecIntronJunction.bed -S /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1.bw  /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-WT1-UV.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1.bw /data/public/cdebes/workspace/scripts/ssl1_chip_seq/HTZ-SSL1-UV.bw  -o matrixHTZ_scale_intron.gz&
plotProfile --matrixFile matrixHTZ_scale_intron.gz -out  profileHTZ_scale_intron.png --perGroup
