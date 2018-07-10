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
Step 3: Use the R code in `mummerParse.R` (which in turn relies on the functions in `pan_transcriptome_fxns.R` to perform local alignments of genes within the Mummer blocks identified in the nucmer output from Step 1.  Alignments are created using a Needleman-Wunsch algorithm, with scores calculated from the Blat .psl file created in Step 2.

Step 4: Run [Assemblytics](assemblytics.com) web analysis tool to process the delta file from nucmer. Then run `parseAssemblytics.R` to combine info with the Mummer blocks, and to determine if there are actually colinear gene alignments with regions called as SVs by Assemblytics.

_Note that some further error checking may be needed (i.e. done by hand) for particularly difficult/repetitive regions._
  
