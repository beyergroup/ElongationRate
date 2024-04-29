#Script to get filtered exonic and intronic annotation. Use the same GTF file as in the alignment of the BAM files.

#Set variable names
gtf_input="Homo_sapiens.GRCh38.90.gtf" #This is the gtf file used for the alignment of the human samples
bed_gene="${gtf_input/gtf/gene.bed}"
bed_exon="${gtf_input/gtf/exon.bed}"
bed_5prim="${gtf_input/gtf/5prim.bed}"
bed_3prim="${gtf_input/gtf/3prim.bed}"
bed_intron="${gtf_input/gtf/intron_junction.bed}"

#Extract positions of different genomic elements

cat $gtf_input | awk '{OFS="\t"}{if ($3 == "gene") {print $1,$4,$5,"intron",$10,$7}}' | bedtools sort > $bed_gene #Bed file of genes
cat $gtf_input | awk '{OFS="\t"}{if (($3 == "exon") && ($5-$4 > 0)){print $1,$4,$5,"exon",$10,$7}}'  | bedtools sort > $bed_exon #Bed file of exons
cat $gtf_input| awk '{OFS="\t"}{if ($3 == "five_prime_utr") {print $1,$4,$5,"5prim",$10,$7}}' | bedtools sort > $bed_5prim #Bed file of  5′ untranslated regions
cat $gtf_input| awk '{OFS="\t"}{if ($3 == "three_prime_utr") {print $1,$4,$5,"3prim",$10,$7}}' | bedtools sort > $bed_3prim #Bed file of 3′ untranslated regions

#To find the positions of the introns, we substract the exons and UTR positions from the gene coordinates. Overlapping introns are merged to remove duplicated regions from the analysis.
cat $bed_gene  | bedtools subtract -s -a stdin -b $bed_exon | bedtools subtract -s -a stdin -b $bed_5prim | bedtools subtract -s -a stdin -b $bed_3prim" | bedtools sort | bedtools merge -i stdin -s -c 4,5,6 -o distinct > $bed intron #Bed file of introns

