#!/usr/bin/bash
#Master script to analyze antisense transcription for Mutzel et al., 2024

bam_dir=''
work_dir=$(pwd)'/'
path=''

help() {
   echo "Generates nascent transcriptome."
   echo
   echo "Syntax: ./master_generate_Analyze_overlaps.sh [-b]"
   echo "options:"
   echo "b     Provide directory containing TT-seq/CUT&Tag BAM files for XXdXic cells (Gjaltema et al., 2022). [mandatory]"
   echo "d     Provides working directory (Standard is current directory)."
   echo "h     Prints this help."
   echo "m     Provide directory containing BS_seq data (from /Antisense_paper/Process_BSseq/) [mandatory]"
   echo "p     Provide path to /Antisense_paper/Analyze_overlaps/. [mandatory]"
   echo "w     Provide directory containing TT-seq BIGWIG files for XXdXic cells (Gjaltema et al., 2022). [mandatory]"
   echo
}

parse_args() {
    case "$1" in
        -b)
            bam_dir="$2"
            ;;
        -d)
            work_dir="$2"
            ;;
        -m)
            bs_dir="$2"
            ;;
        -p)
            path="$2"
            ;;
        -w)
            bigwig_dir="$2"
            ;;
        -h)
            help
            exit 0
            ;;
        *)
            echo "Unknown or badly placed parameter '$1'." 1>&2
            exit 1
            ;;
    esac
}

while [[ "$#" -ge 1 ]]; do
    parse_args "$1" "$2"
    shift; shift
done

if [[ $path == '' ]]
then
	echo -e "Please provide the path to /Antisense_paper/Nascent_transcriptome/ with -p"
  exit 1
fi

if [[ $bam_dir == '' ]]
then
	echo -e "Please provide the path to a directory containing TT-seq, ATAC-seq and CUT&Tag BAM files with -b"
  exit 1
fi

if [[ $bs_dir == '' ]]
then
	echo -e "Please provide the path to a directory containing BS_seq data with -m (from /Antisense_paper/Process_BSseq/)"
  exit 1
fi

if [[ $bigwig_dir == '' ]]
then
	echo -e "Please provide the path to a directory containing TT-seq BIGWIG files with -w"
  exit 1
fi


bigwig_dir=$(realpath $bigwig_dir)'/'
bam_dir=$(realpath $bam_dir)'/'
bs_dir=$(realpath $bs_dir)'/'
work_dir=$(realpath $work_dir)'/'
path=$(realpath $path)'/'


#Select directories containing input files and scripts
script_dir=${path}scripts/
files_dir=${path}files/

mkdir -p ${work_dir}overlaps/
overlaps_dir=${work_dir}overlaps/

mkdir -p ${work_dir}data/
data_dir=${work_dir}data/

mkdir -p ${work_dir}count_tables/
count_dir=${work_dir}count_tables/

mkdir -p ${work_dir}figures/
fig_dir=${work_dir}figures/

echo -e "Plotting TSS sources and number"
Rscript ${script_dir}plot_tss_source.R $files_dir $fig_dir

echo -e "Finding overlaps and plotting basic statistics"
for day in d0 d2 d4
do
  Rscript ${script_dir}find_overlaps.R $day $files_dir $overlaps_dir
  Rscript ${script_dir}desc_overlaps.R $day $bam_dir $overlaps_dir $files_dir $count_dir
done

Rscript ${script_dir}plot_overlap_stats.R $overlaps_dir $count_dir $fig_dir

echo -e "Create control sets, make computeMatrix BEDs and plot stats"
for day in d0 d4
do
  Rscript ${script_dir}make_computeMatrix_bed.R $files_dir $overlaps_dir $count_dir $data_dir $day
  Rscript ${script_dir}plot_control_stats.R $data_dir $files_dir $count_dir $fig_dir $day
done

echo -e "Creating linegraphs"
for day in d0 d4
do
  ${script_dir}linegraph_computematrix.sh $bigwig_dir $data_dir $count_dir $day
  ${script_dir}linegraph_computematrix_antisense.sh $bigwig_dir $data_dir $count_dir $day
  Rscript ${script_dir}plot_linegraph_3prime_overlap.R $data_dir $count_dir $fig_dir $day
  Rscript ${script_dir}plot_linegraph_3prime_overlap_antisense.R $data_dir $count_dir $fig_dir $day
done

echo -e "Calculate CnT enrichment around promoters"
Rscript ${script_dir}cnt_lfc_heatmaps.R $data_dir $bam_dir $fig_dir


echo -e "Calculate and plot methylation around promoters"
Rscript ${script_dir}plot_BSseq.R $bs_dir $data_dir $fig_dir


echo -e "Plotting composite overlaps"
Rscript ${script_dir}create_composite_overlap_set.R $files_dir $overlaps_dir $bam_dir
Rscript ${script_dir}plot_ratio_area.R $overlaps_dir $bam_dir $fig_dir
