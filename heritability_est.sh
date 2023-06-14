#!/bin/bash
# Loop over phenotypes and perform heritability estimation 

while getopts "m:" OPTION
do 
	case "$OPTION" in 
		m) mode=$OPTARG
	esac
done

if [ -z "$mode" ]; then mode="bash"; fi # default mode is bash 

if [ "$mode" != "bash" ] && [ "$mode" != "slurm" ] && [ "$mode" != "dryrun" ]; then
	echo "ERROR: -m flag must be set to \"bash\", \"slurm\", or \"dryrun\""
	exit 
fi 

mkdir -p logs

phenotypes=("HarvWeight" "RetroFat_norm" "EpiFat_norm" \
	"OmentalFat_norm" "Total_AUC" "FastGluc" "FastIns" \
	"REST_EPISODE_COUNT_5" "MOVEMENT_EPISODE_COUNT_5" \
	"VERTICAL_EPISODE_COUNT_5" "FLOAT" "SWIM" "CLIMB" \
	"OpenJunc", "ClosedJunc" "Factor1" "Factor2" "Factor3" \
	"Factor4")

sexes=("Male" "Female") 

for sex in ${sexes[@]}; do
	datapath="derived_data/transformed_data_${sex}-G3package.csv"
	matpath="derived_data/relatedness_matrix_${sex}_PilotThesis.txt"
	for pheno in ${phenotypes[@]}; do
            cmd="Rscript pheno_heritability_est.R --args --datapath=${datapath} --matpath=${matpath} --sex=${sex} --pheno=${pheno} --num_samples=100000 --num_chains=3 --burnin=1000 --thin=10"
		if [ $mode == "slurm" ]; then 
			logfile="logs/${sex}-${pheno}.out"
			sbatch --mem=30G -t 72:00:00 --output=${logfile} --wrap="module add r; ${cmd}"
		fi
		if [ $mode == "bash" ]; then 
			logfile="logs/${sex}-${pheno}.log"
			eval "nohup ${cmd} > ${logfile} 2>&1 < /dev/null &"
		fi
		if [ $mode == "dryrun" ]; then 
			logfile="logs/${sex}-${pheno}.log"
			echo "${cmd} > ${logfile} 2>&1"
		fi
	done
done
