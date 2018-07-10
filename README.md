# Rio
Scripts and command lines used to compare the new Rio (sweet sorghum) genome to the original sorghum reference as well as to identify changes in gene expression associated with sugar accumulation.

## Aligning the Genomes and Filtering Output
Step 1: Run the Nucmer alignment option from MUMmer:
```
nucmer --maxmatch -c 200 -l 100 -b 200 -g 500 BTx623.fa Rio.fa
```
Step 2: Run Blat to create pairwise alignments between transcripts:
```
ref1=SbicolorRiov2.1.primaryTrs.fa
ref2=Sbicolor_454_v3.1.1.transcript_primaryTranscriptOnly.fa 

blat -t=dna -q=dna $ref1 $ref2  Rio_Btx_v3.1.1_Blat.psl
```

