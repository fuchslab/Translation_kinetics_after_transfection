#!/bin/bash

#SBATCH -o slurm_messages/output_%j.txt
#SBATCH -e slurm_messages/error_%j.txt
#SBATCH -J simulating_data
#SBATCH -p icb_interactive
#SBATCH -w icb-iris
#SBATCH --cpus-per-task=1
#SBATCH --mem=1000
#SBATCH -t 0-00:10:00
#SBATCH --nice=10000

# This bash script can be submitted to SLURM with 
# `sbatch bash_scripts/simulating_data_submit.sh <argument>`
# and one of the arguments [dataset1_no_error, dataset1_with_error, 
# dataset2_no_error, dataset2_with_error]

# the config file defines the following environment variables:
# WORK_DIR, PATH_IMG_EXT, CONTAINER_NAME, TMP_DIR, PATH_IMG_COMP
# and the function: start_container
source ./project_config.sh

echo job name: $SLURM_JOB_NAME
echo dataset: $1
echo
echo job id: $SLURM_JOB_ID
echo
echo node name: $SLURMD_NODENAME

if [ ! -d ${PATH_IMG_EXT}/${CONTAINER_NAME} ] ; then
    # create a local folder for your own user
    # -p: don't report error if it already exists
    mkdir -p $TMP_DIR
    
    # extract the container image from a tarball
    ch-tar2dir $PATH_IMG_COMP $PATH_IMG_EXT
fi

# start a charliecloud container;
# the user home directory is mounted automatically; however,
# note that the path to the "user home" folder does not contain "/icb" 
start_container ${WORK_DIR}/bash_scripts/simulating_data_run.sh $1

# move and rename output and error message files
output_message_file="slurm_messages/simulating_data/output_${1}.txt"
error_message_file="slurm_messages/simulating_data/error_${1}.txt"
mkdir -p "slurm_messages/simulating_data"
mv "slurm_messages/output_${SLURM_JOB_ID}.txt" $output_message_file
mv "slurm_messages/error_${SLURM_JOB_ID}.txt" $error_message_file

# once done, remove directory again
#rm -rf ${PATH_IMG_EXT}/${CONTAINER_NAME}
