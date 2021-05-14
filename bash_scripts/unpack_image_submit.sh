#!/bin/bash

#SBATCH --output=slurm_messages/unpacked_%J.txt
#SBATCH --job-name=unpack_container_image
#SBATCH --partition=icb_cpu
#SBATCH --nodelist=icb-neu-002
#SBATCH --cpus-per-task=1
#SBATCH --mem=1000
#SBATCH --time=00-00:10:00
#SBATCH --nice=0


echo job name: $SLURM_JOB_NAME
echo
echo job id: $SLURM_JOB_ID
echo
echo node name: $SLURMD_NODENAME

if [ ! -d "/localscratch/${USER}/r3.6.2_rstan_rmd" ] ; then
    # create a local folder for your own user
    # -p: don't report error if it already exists
    mkdir -p /localscratch/${USER}/tmp
    
    # extract the container image from a tarball
    ch-tar2dir /storage/groups/biostat01/projects/mRNA_transfection_project/container_image/r3.6.2_rstan_rmd.tar.gz /localscratch/${USER}
else
    echo was already unpacked
fi


output_message_file="slurm_messages/unpacked_on_${SLURMD_NODENAME}.txt"
mv "slurm_messages/unpacked_${SLURM_JOB_ID}.txt" $output_message_file
