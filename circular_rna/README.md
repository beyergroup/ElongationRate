For circular RNA, we aligned the reads using STAR version 2.5.1b with the following parameters: 

--chimSegmentMin 15 
--outSJfilterOverhangMin 15 15 15 15 
--alignSJoverhangMin 15
--alignSJDBoverhangMin 15 
--seedSearchStartLmax 30 
--outFilterMultimapNmax 20
--outFilterScoreMin 1 
--outFilterMatchNmin 1
--outFilterMismatchNmax 2
--chimScoreMin 15
--chimScoreSeparation 10 
--chimJunctionOverhangMin 15

We then extracted back spliced reads from the STAR chimeric output file and normalized the number of back spliced reads by the sum of back spliced (BS<sub>i</sub>) and spliced reads from linear transcripts (S1<sub>i</sub>, S2<sub>i</sub>) for an exon i, as explained in the paper. Thus, this score quantifies the percent of transcripts from this locus that resulted in circular RNA. Finally, we quantified the significance of the average change in circular RNA formation between two conditions using the Wilcoxon rank test.

