# Antisense transcription can induce expression memory via stable promoter repression
Verena Mutzel, Till Schwämmle, Svearike Oeverdieck, Lucia Librenjak, Benedikt Boesen, Melissa Bothe, Rutger A.F. Gjaltema, Ilona Dunkel, Gemma Noviello and Edda G. Schulz 

Data and code used to perform computational analyses in Mutzel et al., 2024. (https://www.biorxiv.org/content/10.1101/2024.03.06.583761v1)


## Abstract
The capacity of cells to retain a memory of previous signals enables them to adopt unique cell fates and adjust to their surrounding environment. The underlying gene expression memory can arise from mutual repression of two genes, forming a toggle switch. Such mutual repression may occur at antisense loci, where two convergently oriented genes repress each other in cis. Under which conditions antisense transcription can generate expression memory remains poorly understood. To address this question, we combine mathematical modeling, genomics and a synthetic biology approach. Through simulations we show that stable memory can emerge, if both genes in an antisense pair transcribe through the convergent promoter and induce a stable repressive chromatin state. Genome-wide analysis of nascent transcription further supports antisense-mediated promoter repression with promoter-overlapping antisense gene pairs exhibiting mutually exclusive expression. Through constructing a synthetic antisense locus in mouse embryonic stem cells (mESCs) we then show that such a locus architecture can indeed maintain a memory of a transient stimulus. Mutual repression and the capacity for memory formation are elevated, when mESCs differentiate, showing that epigenetic memory is a cell type-specific property. Our finding that stem cells adapt their ability to remember stimuli as they differentiate might help to elucidate how stemness is maintained.


## Description
Align_methylation: Contains instructions, code and files to reproduce WGBS data alignment (GSE253792). (Code by Melissa Bothe)

Analyze_overlaps: Contains instructions, code and files to reproduce genome-wide analysis of antisense transcription for figures 3-4 and supplemental figures 1-2. (Code by Till Schwämmle)

Nascent_transcriptome: Contains instructions, code and files to reproduce the generation of a nascent transcriptome based on TT-seq data from Gjaltema et al., 2022 (GSE167358). (Code by Till Schwämmle, with help from Melissa Bothe)

AS_Model: Contains all code that was used to perform and analyze the simulations of the antisense model (Code by Verena Mutzel)

synthetic_AS_locus: Contains the data and code of the flow cytometry measurements of the synthetic antisense locus



