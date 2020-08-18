# MyRNAPipe
This pipeline includes the following contents:
1. Mapping from fastq file to Transcripts counts by Salmon
2. Import Salmon transcripts counts and merge them into gene counts by tximport
3. Normalization by DESeq2
4. Basic PCA analysis to check abnormalities using vst transformed counts from DESeq2
5. DEG analysis by DESeq2
6. GSEA analysis by ClusterProfile
7. Permutation based gene set enrichment analysis ClusterProfile
8. View Pathway map using Pathview
