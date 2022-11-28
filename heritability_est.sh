#!/bin/bash
# Loop over phenotypes and perform heritability estimation 

while getopts "m:" OPTION
do 
	case "$OPTION" in 
		m) mode=$OPTARG
	esac
done

if [ -z "$mode" ]; then mode="bash"; fi # default mode is bash 

if [ "$mode" != "bash" ] && [ "$mode" != "slurm" ]; then
	echo "ERROR: -m flag must be set to \"bash\" or \"slurm\"" 
	exit 
fi 

mkdir -p logs

phenotypes=("Total_AUC" "HarvWeight" "RetroFat_norm" "FastGluc" \
	"FastIns" "SWIM" "CLIMB" "FLOAT" "REST_EPISODE_COUNT_5" \
	"MOVEMENT_EPISODE_COUNT_5" "VERTICAL_EPISODE_COUNT_5" \
	"EpiFat_norm" "OmentalFat_norm" "ClosedJunc" "OpenJunc")

sexes=("male" "female") 

for sex in ${sexes[@]}; do
	datapath="derived_data/transformed_data_${sex}-G3package.csv"
	matpath="derived_data/relatedness_matrix_MetBehThesis_${sex}.txt"
	for pheno in ${phenotypes[@]}; do
		if [ $mode == "slurm" ]; then 
			logfile="logs/${sex}-${pheno}.out"
			sbatch --mem=30G -t 24:00:00 --output=${logfile} --wrap="module add r; Rscript pheno_heritability_est.R --args --datapath=${datapath} --matpath=${matpath} --sex=${sex} --pheno=${pheno}"
		fi
		if [ $mode == "bash" ]; then 
			logfile="logs/${sex}-${pheno}.log"
			nohup Rscript pheno_heritability_est.R --args --datapath=${datapath} --matpath=${matpath} --sex=${sex} --pheno=${pheno} > ${logfile} 2>&1 < /dev/null &
		fi
	done
done