# Antisense transcription can induce expression memory via stable promoter repression
Verena Mutzel, Till Schwämmle, Melissa Bothe, Rutger A.F. Gjaltema, Benedikt Boesen, Lucia Librenjak, Svearike Oeverdiek, Ilona Dunkel, Gemma Noviello and Edda G. Schulz 

Data and code used to perform computational analyses in Mutzel et al., 2024. (LINK/TO/PAPER)


## Abstract
Antisense transcription is widespread across genomes. Despite this, a systematic understanding of its biological function is lacking. Here, we set out to characterize the regulatory responses that can be elicited by antisense transcription. To address this question, we stably integrate a synthetic reporter construct with an inducible antisense promoter into mouse embryonic stem cells (mESCs) and measure stimulus-response curves and hysteresis. We developed a modeling framework describing transcription from initiation to transcript degradation to systematically characterize the regulatory responses elicited by antisense transcription. The model predicts that ultrasensitivity and hysteresis in the dose-response curve can arise for antisense pairs with the ability to stably repress the convergent promoter through the recruitment of epigenetic modifications. When we adapt our experimental system to fit these conditions, it indeed displays expression memory over timescales of several days. 

To corroborate our findings in a comprehensive manner, we identify antisense transcription across the genome using nascent RNA-sequencing data. We find that globally repression via antisense transcription is induced by promoter repression, but not primarily via transcriptional interference. Additionally, we find that high levels of DNA methylation at the promoter, but not of repressive histone marks, such as H3K27me3 or H3K9me3, correlate with repression via antisense transcription.  

Our unbiased systematic analysis establishes a potential general functional role of antisense transcription in converting graded input signals into non-linear hysteretic transcriptional responses. This allows genes to respond in an all-or-nothing manner, and to remember transient stimuli for extended periods and lock in alternative expression states allowing cells to establish and maintain molecularly distinct cell fates.


## Description
Align_methylation: Contains instructions, code and files to reproduce WGBS data alignment (GSEXXX). (Code by Melissa Bothe)


Analyze_overlaps: Contains instructions, code and files to reproduce genome-wide analysis of antisense transcription for figures 3-4 and supplemental figures 1-2. (Code by Till Schwämmle)

Nascent_transcriptome: Contains instructions, code and files to reproduce the generation of a nascent transcriptome based on TT-seq data from Gjaltema et al., 2022 (GSE167358). (Code by Till Schwämmle and Melissa Bothe)


