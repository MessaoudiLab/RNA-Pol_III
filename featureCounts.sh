#!/bin/bash
  
featureCounts -T 10 \
--countReadPairs \
-t CDS,exon,gene,transcript,start_codon,stop_codon,five_prime_utr,three_prime_utr \
-p \
-a /rhome/evance/bigdata/rna_seq/scripts/nextflow/reference/homo_sapiens_grch38/Homo_sapiens.GRCh38.104.gtf \
-o test.txt \
/rhome/evance/bigdata/rna_seq/scripts/nextflow/temp_test/samtools_sortedByName/aligned_sortedByName_bams/test-READ_GRCh38_Aligned.sortedByCoord.out.bam
