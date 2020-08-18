# MyRNAPipe
This pipeline includes the following contents:
1. Mapping from fastq file to Transcripts counts by Salmon
2. Import Salmon transcripts counts and merge them into gene counts by tximport
3. Normalization by DESeq2
4. Basic PCA analysis to check abnormalities
5. DEG analysis
6. GSEA analysis
7. Permutation based gene set enrichment analysis
