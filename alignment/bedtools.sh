#Script to get filtered exonic and intronic annotation

#For introns:

cat ../genomes/$spe/intron_$anot.bed | bedtools sort |  bedtools subtract -s -a stdin -b ../genomes/$spe/exon_$anot\_c.bed | bedtools subtract -s -a stdin -b ../genomes/$spe/3prim_$anot.bed | bedtools subtract -s -a stdin -b ../genomes/$spe/5prim_$anot.bed | bedtools sort | bedtools merge -i stdin -s -c 4,5,6 -o distinct > ../genomes/$spe/intron_junction_$anot.bed

#For exons:

cat ../genomes/$spe/exon_$anot.bed | bedtools sort |  bedtools subtract -s -a stdin -b ../genomes/$spe/intron_$anot\_c.bed | bedtools subtract -s -a stdin -b ../genomes/$spe/3prim_$anot.bed | bedtools subtract -s -a stdin -b ../genomes/$spe/5prim_$anot.bed | bedtools sort | bedtools merge -i stdin -s -c 4,5,6 -o distinct > ../genomes/$spe/exon_junction_$anot.bed
