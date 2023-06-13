#!/bin/bash
#SBATCH -p q_cgpu
#SBATCH -N 1
#SBATCH -G 1
#SBATCH -n 10
#SBATCH --job-name=RoseTTAFold
# usage: sbatch run_e2e_slurm.sh -f <path/to/fasta/file> -o <path/to/output/dir>
while getopts f:o: flag
do
    case "${flag}" in
        f) fasta=${OPTARG};;
        o) output_dir=${OPTARG};;
    esac
done

# set base repository path
REPO_BASE=/data/users/rolz/repos/RoseTTAFold_mod

# make the script stop when error (non-true exit code) is occured
set -e

############################################################
# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
# __conda_setup="$('conda' 'shell.bash' 'hook' 2> /dev/null)"
# eval "$__conda_setup"
# unset __conda_setup
# <<< conda initialize <<<
############################################################

#SCRIPT=`realpath -s $0`
#export PIPEDIR=`dirname $SCRIPT`
export PIPEDIR=$REPO_BASE

CPU="$SLURM_JOB_CPUS_PER_NODE" # number of CPUs to use
MEM="64" # max memory (in GB)

# Inputs:
IN="$fasta"                # input.fasta
WDIR=`realpath -s $output_dir`  # working/outputs folder


LEN=`tail -n1 $IN | wc -m`

mkdir -p $WDIR/log

source /data/users/rolz/anaconda3/bin/activate /data/users/rolz/anaconda3/envs/RoseTTAFold
conda activate RoseTTAFold
############################################################
# 1. generate MSAs
############################################################
if [ ! -s $WDIR/t000_.msa0.a3m ]
then
    echo "Running HHblits"
    $PIPEDIR/input_prep/make_msa.sh $IN $WDIR $CPU $MEM > $WDIR/log/make_msa.stdout 2> $WDIR/log/make_msa.stderr
fi


############################################################
# 2. predict secondary structure for HHsearch run
############################################################
if [ ! -s $WDIR/t000_.ss2 ]
then
    echo "Running PSIPRED"
    $PIPEDIR/input_prep/make_ss.sh $WDIR/t000_.msa0.a3m $WDIR/t000_.ss2 > $WDIR/log/make_ss.stdout 2> $WDIR/log/make_ss.stderr
fi


############################################################
# 3. search for templates
############################################################
DB="$PIPEDIR/pdb100_2021Mar03/pdb100_2021Mar03"
if [ ! -s $WDIR/t000_.hhr ]
then
    echo "Running hhsearch"
    HH="hhsearch -b 50 -B 500 -z 50 -Z 500 -mact 0.05 -cpu $CPU -maxmem $MEM -aliw 100000 -e 100 -p 5.0 -d $DB"
    cat $WDIR/t000_.ss2 $WDIR/t000_.msa0.a3m > $WDIR/t000_.msa0.ss2.a3m
    $HH -i $WDIR/t000_.msa0.ss2.a3m -o $WDIR/t000_.hhr -atab $WDIR/t000_.atab -v 0 > $WDIR/log/hhsearch.stdout 2> $WDIR/log/hhsearch.stderr
fi


############################################################
# 4. end-to-end prediction
############################################################
if [ ! -s $WDIR/t000_.e2e.pdb ]
then
    echo "Running end-to-end prediction"
    python $PIPEDIR/network/predict_e2e.py \
        -m $PIPEDIR/weights \
        -i $WDIR/t000_.msa0.a3m \
        -o $WDIR/t000_.e2e \
        --hhr $WDIR/t000_.hhr \
        --atab $WDIR/t000_.atab \
        --db $DB 1> $WDIR/log/network.stdout 2> $WDIR/log/network.stderr
fi


source /data/users/rolz/anaconda3/bin/activate /data/users/rolz/anaconda3/envs/folding
############################################################
# 5. build side-chains by running FastRelax
############################################################
if [ ! -s $WDIR/t000_.e2e_relaxed.pdb ]
then
    echo "Running fast relax"
    python $PIPEDIR/folding/postprocess_structure.py \
        -i $WDIR/t000_.e2e.pdb \
        -o $WDIR/t000_.e2e_relaxed.pdb \
        -p $CPU \
        -n 10 1> $WDIR/log/postprocess.stdout 2> $WDIR/log/postprocess.stderr
fi
echo "Done"
