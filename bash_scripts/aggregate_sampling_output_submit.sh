#!/bin/bash

#SBATCH --output=slurm_messages/output_%j.txt
#SBATCH --error=slurm_messages/error_%j.txt
#SBATCH --job-name=aggregate_sampling_output
#SBATCH -p cpu_p
#SBATCH -x ibis-ceph-[002-006,008-009,017-018],ibis216-010-[001-004,007,020-023,025,027-029,031,033-037,064,068-071],ibis216-224-011,icb-rsrv[05-06],ibis216-010-051,ibis216-224-010,ibis-ceph-[010-016,019],ibis216-010-[011-012,024,026,030,032],icb-rsrv08 # remaining: icb-neu*
#SBATCH --cpus-per-task=2
#SBATCH --mem=1gb
#SBATCH --time=0-08:00:00
#SBATCH --nice=10000

# This bash script can be submitted to SLURM with 
# `sbatch bash_scripts/aggregate_sampling_output_submit.sh <argument1> <argument1>`
# and one of the arguments [experimental_data_eGFP,experimental_data_d2eGFP]
# and [SDE,ODE]

# the config file defines the following environment variables:
# WORK_DIR, PATH_IMG_EXT, CONTAINER_NAME, TMP_DIR, PATH_IMG_COMP
# and the function: start_container
source ./project_config.sh

echo job name: $SLURM_JOB_NAME
echo dataset: $1
echo
echo model type: $2
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

# only start the Stan sampling if there is no stanfit_object for the corresponding index yet
if [ ! -f "intermediate_output_files/aggregate_sampling_output/${1}/aggregated_output_${2}.Rdata" ] ; then
    start=`date +%s`

    # start the charliecloud container;
    start_container ${WORK_DIR}/bash_scripts/aggregate_sampling_output_run.sh $@

    # move and rename output and error message files
    output_message_file="slurm_messages/aggregate_sampling_output/output_${1}_${2}.txt"
    error_message_file="slurm_messages/aggregate_sampling_output/error_${1}_${2}.txt"
    mv "slurm_messages/output_${SLURM_JOB_ID}.txt" $output_message_file
    mv "slurm_messages/error_${SLURM_JOB_ID}.txt" $error_message_file

    end=`date +%s`
    echo
    echo
    echo Duration in seconds: $(echo "$end - $start" | bc -l)
else
    rm "slurm_messages/output_${SLURM_JOB_ID}.txt"
    rm "slurm_messages/error_${SLURM_JOB_ID}.txt"
fi

#rm -rf ${PATH_IMG_EXT}/${CONTAINER_NAME}
