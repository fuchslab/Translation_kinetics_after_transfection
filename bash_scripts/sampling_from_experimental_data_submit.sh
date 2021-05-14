#!/bin/bash

#SBATCH -o slurm_messages/%x/output_sampling_from_experimental_data.txt
#SBATCH -e slurm_messages/%x/error_sampling_from_experimental_data.txt
#SBATCH -J sampling_from_experimental_data
#SBATCH -p icb_interactive
#SBATCH -w icb-iris
#SBATCH --cpus-per-task=1
#SBATCH --mem=1000
#SBATCH -t 0-00:10:00
#SBATCH --nice=10000


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
fi

# start a charliecloud container;
# the user home directory is mounted automatically; however,
# note that the path to the "user home" folder does not contain "/icb" 
ch-run -b /localscratch:/localscratch/  -b /storage/groups/:/storage/groups  /localscratch/${USER}/r3.6.2_rstan_rmd/ -- /bin/bash /storage/groups/biostat01/projects/mRNA_transfection_project/bash_scripts/sampling_from_experimental_data_run.sh

# once done, remove directory again
# rm -rf /localscratch/${USER}/
