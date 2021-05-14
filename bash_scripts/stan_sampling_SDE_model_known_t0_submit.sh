#!/bin/bash

#SBATCH --output=slurm_messages/output_%A_%a.txt
#SBATCH --error=slurm_messages/error_%A_%a.txt
#SBATCH --job-name=stan_sampling_SDE_model_known_t0
#SBATCH --partition=icb_cpu
#SBATCH --exclude=ibis-ceph-[002-006,008-009,017-018],ibis216-010-[001-004,007,020-023,025,027-029,031,033-037,064,068-071],ibis216-224-011,icb-rsrv[05-06],ibis216-010-051,ibis216-224-010,ibis-ceph-[010-016,019],ibis216-010-[011-012,024,026,030,032],icb-rsrv08 # remaining: icb-neu*
#SBATCH --cpus-per-task=16
#SBATCH --mem=8gb
#SBATCH --time=2-00:00:00
#SBATCH --nice=10000
#SBATCH --array=1-100

# This bash script can be submitted to SLURM with 
# `sbatch bash_scripts/stan_sampling_SDE_model_mult_error_known_t0_submit.sh <argument>`
# and one of the arguments [experimental_data_eGFP,experimental_data_d2eGFP,
# simulated_data_dataset1_with_error,simulated_data_dataset1_no_error,
# simulated_data_dataset2_with_error,simulated_data_dataset2_no_error]

# the config file defines the following environment variables:
# WORK_DIR, PATH_IMG_EXT, CONTAINER_NAME, TMP_DIR, PATH_IMG_COMP
# and the function: start_container
source ./project_config.sh

echo job name: $SLURM_JOB_NAME
echo dataset: $1
echo
echo job id: $SLURM_JOB_ID
echo array job id: $SLURM_ARRAY_JOB_ID
echo array task id: $SLURM_ARRAY_TASK_ID
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
if [ ! -f "intermediate_output_files/stanfit_objects/${1}_SDE/stanfit_object_$SLURM_ARRAY_TASK_ID.rds" ] ; then
  # only start the Stan sampling for SDE model if the stanfit_object for the corresponding index of the ODE model exists
  if [ -f "intermediate_output_files/stanfit_objects/${1}_ODE/stanfit_object_$SLURM_ARRAY_TASK_ID.rds" ] ; then
      start=`date +%s`
  
      # start the charliecloud container;
      start_container ${WORK_DIR}/bash_scripts/stan_sampling_SDE_model_known_t0_run.sh $1 $SLURM_ARRAY_TASK_ID
  
      # move and rename output and error message files
      mkdir -p "slurm_messages/$SLURM_JOB_NAME/${1}_SDE/"
      output_message_file="slurm_messages/$SLURM_JOB_NAME/${1}_SDE/output_$SLURM_ARRAY_TASK_ID.txt"
      error_message_file="slurm_messages/$SLURM_JOB_NAME/${1}_SDE/error_$SLURM_ARRAY_TASK_ID.txt"
      mv "slurm_messages/output_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.txt" $output_message_file
      mv "slurm_messages/error_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.txt" $error_message_file
  
      end=`date +%s`
      echo
      echo
      echo Duration in seconds: $(echo "$end - $start" | bc -l)
  else
      echo "intermediate_output_files/stanfit_objects/${1}_ODE/stanfit_object_$SLURM_ARRAY_TASK_ID.rds does not exist"
  fi
else
    rm "slurm_messages/output_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.txt"
    rm "slurm_messages/error_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.txt"
fi



