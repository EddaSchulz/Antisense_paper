#!/usr/bin/bash

#Changed designate TSS from +1 to +2 for comp
#Intersect TSS with complete gene body instead only with tss
#

#Specify a working directory
work_dir=/project/agsgpa/Till/verena_antisense/pipeline_test/final_pipeline/

#Specify path to ChromHMM software
CHROMHMM=/project/agsgpa/NGS_Tools/ChromHMM/ChromHMM.jar

#Select path to directory containing input BAM files
bam_dir=${work_dir}bam_files/

#Select directories containing input files and scripts
script_dir=${work_dir}scripts/
files_dir=${work_dir}files/

#Create additional required directories
mkdir -p ${work_dir}data/
data_dir=${work_dir}data/

mkdir -p ${work_dir}out_chromHMM/
chromHMM_dir=${work_dir}out_chromHMM/

mkdir -p ${work_dir}binarize_chromHMM/
binarize_dir=${work_dir}binarize_chromHMM/

mkdir -p ${work_dir}raw_assembly/
raw_dir=${work_dir}raw_assembly/

mkdir -p ${work_dir}binned_assembly/
binned_dir=${work_dir}binned_assembly/

mkdir -p ${work_dir}transcript_source/
source_dir=${work_dir}transcript_source/

mkdir -p ${work_dir}final_assembly/
final_dir=${work_dir}final_assembly/

#Removing reads in blacklisted regions from TT-seq file
for day in d0 d2 d4
do
	echo -e "Removing Blacklisted regions from $day bam files"
	bedtools intersect -v -a ${bam_dir}XX_TT_$day\_plus.bam -b ${files_dir}mm10-blacklist.v2.bed > ${bam_dir}XX_TT_bl_$day\_plus.bam
	samtools index ${bam_dir}XX_TT_bl_$day\_plus.bam

	bedtools intersect -v -a ${bam_dir}XX_TT_$day\_minus.bam -b ${files_dir}mm10-blacklist.v2.bed > ${bam_dir}XX_TT_bl_$day\_minus.bam
	samtools index ${bam_dir}XX_TT_bl_$day\_minus.bam
done


#Uses deeptools to create count tables for the TT-seq data
for day in d0 d2 d4
do
  echo -e "Counting TTseq Reads with multiBamSummary for $day"
  prun python3 multiBamSummary bins -b ${bam_dir}XX_TT_bl_$day\_plus.bam ${bam_dir}XX_TT_bl_$day\_minus.bam --genomeChunkSize 2407883318 \
  -l TT_XX_$day\_plus TT_XX_$day\_minus -bs 200 -o ${data_dir}XX_$day\_binned_counts.npz --outRawCounts ${data_dir}XX_$day\_binned_counts.txt -p 30 --centerReads
  sed -i 's/[#'\'']//g' ${data_dir}XX_$day\_binned_counts.txt
  sed '1d' ${data_dir}XX_$day\_binned_counts.txt | sort -k 1,1 -k2,2n - > ${data_dir}XX_$day\_binned_counts_sorted.txt
  rm ${data_dir}XX_$day\_binned_counts.npz
done

#Uses GENOSTAN package to assign transcription states based on the TTseq counts
for day in d0 d2 d4
do
  echo -e "Running Genostan for $day"
  Rscript ${script_dir}TT_genostan.R $day $data_dir
done


echo -e "Binarizing bam files for ChromHMM"
java -mx4000M -jar $CHROMHMM BinarizeBam -paired ${files_dir}mm10_chrom_sizes.txt $bam_dir ${files_dir}chromHMM_table.txt \
  ${work_dir}binarize_chromHMM/

echo -e "Running ChromHMM for CUT&Tag"
java -mx10000M -jar $CHROMHMM LearnModel -p 40 -s 222 $binarize_dir $chromHMM_dir 6 mm10


echo -e 'Reducing ChromHMM tracks'
${script_dir}reduce_chromHMM.py $chromHMM_dir $data_dir

echo -e 'Designating candidate transcription start sites according to TT-seq'
${script_dir}designate_tss.py $data_dir

echo -e 'Creating BED files of candidate TSS regions overlapping to CnT enhancer/promoter, CAGE peak or annotated TSS'
for day in d0 d2 d4
do
  grep -P '\t1\t0\t' ${data_dir}$day\_XX_chromHMM.bed | cut -f 1,2,3 | bedtools merge -i - | \
    bedtools sort -i - | sed 's/$/\tpromoter/' - > ${data_dir}$day\_XX_promoter_chromHMM.bed
  grep -P '\t2\t0\t' ${data_dir}$day\_XX_chromHMM.bed | cut -f 1,2,3 | bedtools merge -i - | \
    bedtools sort -i - | sed 's/$/\tenhancer/' - > ${data_dir}$day\_XX_enhancer_chromHMM.bed
  cat ${data_dir}$day\_XX_promoter_chromHMM.bed ${data_dir}$day\_XX_enhancer_chromHMM.bed | sort -k 1,1 -k2,2n - > ${data_dir}$day\_XX_tss_chromHMM.bed
  for strand in plus minus
  do
    bedtools intersect -wo -a ${data_dir}XX_$day\_pot_tss_$strand\.bed -b ${data_dir}$day\_XX_tss_chromHMM.bed | cut -f 1,2,3,7 \
      > ${data_dir}$day\_XX_pot_tss_$strand\_chromHMM_overlap.bed
    bedtools intersect -wo -a ${data_dir}XX_$day\_pot_tss_$strand\.bed -b ${files_dir}FANTOM5_CAGE_peaks_mm10_$strand\.bed | cut -f 1,2,3 | \
      sed 's/$/\tCAGE/' - > ${data_dir}$day\_XX_pot_tss_$strand\_CAGE_overlap.bed
    bedtools intersect -wo -a ${data_dir}XX_$day\_pot_tss_$strand\.bed -b ${files_dir}GENCODE_tss_$strand\.bed | cut -f 1,2,3 | \
     sed 's/$/\tGENCODE/' - > ${data_dir}$day\_XX_pot_tss_$strand\_GENCODE_overlap.bed
    cat ${data_dir}$day\_XX_pot_tss_$strand\_chromHMM_overlap.bed ${data_dir}$day\_XX_pot_tss_$strand\_CAGE_overlap.bed \
      ${data_dir}$day\_XX_pot_tss_$strand\_GENCODE_overlap.bed | bedtools sort -i - | bedtools merge -i - | sed 's/$/\ttss/' - \
      > ${data_dir}$day\_XX_tss_$strand\.bed
  done
done


echo -e 'Combining TSS beds with transcribed regions'
for day in d0 d2 d4
do
  grep -P 'plus' ${data_dir}Genostan_XX_$day\_plus_reduced.bed | bedtools merge -i - | \
   bedtools sort -i - | bedtools subtract -a - -b ${data_dir}$day\_XX_tss_plus.bed | sed 's/$/\tplus/' - > ${data_dir}$day\_XX_plus_trans.bed
  grep -P 'minus' ${data_dir}Genostan_XX_$day\_minus_reduced.bed | bedtools merge -i - | \
   bedtools sort -i - | bedtools subtract -a - -b ${data_dir}$day\_XX_tss_minus.bed | sed 's/$/\tminus/' - > ${data_dir}$day\_XX_minus_trans.bed
  cat ${data_dir}$day\_XX_tss_plus.bed ${data_dir}$day\_XX_plus_trans.bed | sort -k 1,1 -k2,2n - > ${data_dir}$day\_XX_plus_raw.bed
  cat ${data_dir}$day\_XX_tss_minus.bed ${data_dir}$day\_XX_minus_trans.bed | sort -k 1,1 -k2,2n - > ${data_dir}$day\_XX_minus_raw.bed
done

echo -e 'Creating RGB beds of raw tracks'
${script_dir}create_raw_rgb.py $data_dir $raw_dir

echo -e 'Marking track file with annotated genes'
for day in d0 d2 d4
do
  for strand in plus minus
  do
    grep 'tss'  ${data_dir}$day\_XX_$strand\_raw.bed | bedtools intersect -wo -a - -b ${files_dir}GENCODE_tss_$strand\.bed | cut -f 1,2,3,4,8,9,10 |
      sort -k1,1 -k2,2n -k5,5 - | uniq > ${data_dir}$day\_XX_$strand\_tss_annot.bed
    grep "$strand" ${data_dir}$day\_XX_$strand\_raw.bed| bedtools intersect -wo -a - -b ${files_dir}GENCODE_gene_$strand\.bed | cut -f 1,2,3,4,8,9,10 |
      sort -k1,1 -k2,2n -k5,5 - | uniq > ${data_dir}$day\_XX_$strand\_gene_annot.bed
  done
  grep 'tss' ${data_dir}$day\_XX_plus_raw.bed | bedtools intersect -v -a - -b ${files_dir}GENCODE_gene_plus.bed | sed 's/$/\t+\t.\t./' - \
    > ${data_dir}$day\_XX_plus_tss_no_annot.bed
  grep 'plus' ${data_dir}$day\_XX_plus_raw.bed | bedtools intersect -v -a - -b ${files_dir}GENCODE_gene_plus.bed | sed 's/$/\t+\t.\t./' - \
    > ${data_dir}$day\_XX_plus_gene_no_annot.bed
  grep 'tss' ${data_dir}$day\_XX_minus_raw.bed | bedtools intersect -v -a - -b ${files_dir}GENCODE_gene_minus.bed | sed 's/$/\t-\t.\t./' - \
    > ${data_dir}$day\_XX_minus_tss_no_annot.bed
  grep 'minus' ${data_dir}$day\_XX_minus_raw.bed | bedtools intersect -v -a - -b ${files_dir}GENCODE_gene_minus.bed | sed 's/$/\t-\t.\t./' - \
    > ${data_dir}$day\_XX_minus_gene_no_annot.bed
  for strand in plus minus
  do
  cat ${data_dir}$day\_XX_$strand\_tss_annot.bed ${data_dir}$day\_XX_$strand\_gene_annot.bed ${data_dir}$day\_XX_$strand\_tss_no_annot.bed \
    ${data_dir}$day\_XX_$strand\_gene_no_annot.bed | sort -k1,1 -k2,2n -k5,5 - > ${data_dir}$day\_XX_$strand\_full_annot.bed
  done
done

echo -e 'Removing sporadic transcription and merging bins in annotated transcripts'
for day in d0 d2 d4
do
  ${script_dir}remove_sporadic_transcripts.py $day $data_dir $binned_dir
  ${script_dir}make_tss_bins.py $day $binned_dir $data_dir
  Rscript ${script_dir}find_tss_source.R $day $data_dir $source_dir
done

echo -e 'Counting reads in TSS and TTS bins'
for day in d0 d2 d4
do
  for strand in plus minus
  do
    prun python3 multiBamSummary BED-file --BED ${data_dir}XX_$day\_$strand\_tss_bins.bed -b ${bam_dir}XX_TT_bl_$day\_$strand\.bam \
      -l TT_$strand -o ${data_dir}XX_$day\_tss_$strand\_binned_counts.npz --outRawCounts ${data_dir}XX_$day\_tss_$strand\_binned_counts.txt \
      -p 20 --centerReads
      prun python3 multiBamSummary BED-file --BED ${data_dir}XX_$day\_$strand\_tts_bins.bed -b ${bam_dir}XX_TT_bl_$day\_$strand\.bam \
        -l TT_$strand -o ${data_dir}XX_$day\_tts_$strand\_binned_counts.npz --outRawCounts ${data_dir}XX_$day\_tts_$strand\_binned_counts.txt \
        -p 20 --centerReads
    sed -i 's/[#'\'']//g' ${data_dir}XX_$day\_tss_$strand\_binned_counts.txt
    sed -i 's/[#'\'']//g' ${data_dir}XX_$day\_tts_$strand\_binned_counts.txt
    rm  ${data_dir}XX_$day\_tss_$strand\_binned_counts.npz
    rm  ${data_dir}XX_$day\_tts_$strand\_binned_counts.npz
  done
done

echo -e 'Trimming transcript ends'
for day in d0 d2 d4
do
  Rscript ${script_dir}trim_transcript_ends.R $day $data_dir $final_dir
done
