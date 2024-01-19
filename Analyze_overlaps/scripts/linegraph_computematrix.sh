#!/usr/bin/bash

tt_dir=$1
bed_dir=$2
count_dir=$3
day=$4

cd $work_dir

echo -e "Counting 3'overlaps and controls using computeMatrix"
prun python3 computeMatrix scale-regions -S ${tt_dir}XX_TT_$day\_plus.bw -R ${bed_dir}XX_$day\_3prime_case_free_plus.bed -m 120 \
-bs 10 --outFileNameMatrix ${count_dir}XX_$day\_3prime_case_free_plus.txt -o ${count_dir}trash.txt -p 20
prun python3 computeMatrix scale-regions -S ${tt_dir}XX_TT_$day\_minus.bw -R ${bed_dir}XX_$day\_3prime_case_free_minus.bed -m 120 \
-bs 10 --outFileNameMatrix ${count_dir}XX_$day\_3prime_case_free_minus.txt -o ${count_dir}trash.txt -p 20

prun python3 computeMatrix scale-regions -S ${tt_dir}XX_TT_$day\_plus.bw -R ${bed_dir}XX_$day\_3prime_controls_free_plus.bed -m 120 \
-bs 10 --outFileNameMatrix ${count_dir}XX_$day\_3prime_controls_free_plus.txt -o ${count_dir}trash.txt -p 20
prun python3 computeMatrix scale-regions -S ${tt_dir}XX_TT_$day\_minus.bw -R ${bed_dir}XX_$day\_3prime_controls_free_minus.bed -m 120 \
-bs 10 --outFileNameMatrix ${count_dir}XX_$day\_3prime_controls_free_minus.txt -o ${count_dir}trash.txt -p 20


prun python3 computeMatrix scale-regions -S ${tt_dir}XX_TT_$day\_plus.bw -R ${bed_dir}XX_$day\_3prime_case_overlap_plus.bed -m 80 \
-bs 10 --outFileNameMatrix ${count_dir}XX_$day\_3prime_case_overlap_plus.txt -o ${count_dir}trash.txt -p 20
prun python3 computeMatrix scale-regions -S ${tt_dir}XX_TT_$day\_minus.bw -R ${bed_dir}XX_$day\_3prime_case_overlap_minus.bed -m 80 \
-bs 10 --outFileNameMatrix ${count_dir}XX_$day\_3prime_case_overlap_minus.txt -o ${count_dir}trash.txt -p 20

prun python3 computeMatrix scale-regions -S ${tt_dir}XX_TT_$day\_plus.bw -R ${bed_dir}XX_$day\_3prime_controls_overlap_plus.bed -m 80 \
-bs 10 --outFileNameMatrix ${count_dir}XX_$day\_3prime_controls_overlap_plus.txt -o ${count_dir}trash.txt -p 20
prun python3 computeMatrix scale-regions -S ${tt_dir}XX_TT_$day\_minus.bw -R ${bed_dir}XX_$day\_3prime_controls_overlap_minus.bed -m 80 \
-bs 10 --outFileNameMatrix ${count_dir}XX_$day\_3prime_controls_overlap_minus.txt -o ${count_dir}trash.txt -p 20


rm ${count_dir}trash.txt
