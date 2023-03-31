# Mismatch detection

Mismatch detection was performed using the tool rnaseqmut (https://github.com/davidliwei/rnaseqmut ), which detects mutations from the NM tag of BAM files. To avoid detection of RNA editing or DNA damage-based events, we only considered genomic positions with only one mismatch detected (that is, occurring in only one single read). Reads with indels were excluded and only mismatches with a distance of more than four from the beginning and the end of the read were considered. A coverage-level filter was applied so that only bases covered by at least 100 reads were kept. 

A substantial number of mismatches may result from technical sequencing errors. However, as young and old samples were always handled together in the same batch, we can exclude that consistent differences in the number of mismatches are due to technical biases. The fraction of RNA editing events is generally relatively low and not expected to globally increase with age.
