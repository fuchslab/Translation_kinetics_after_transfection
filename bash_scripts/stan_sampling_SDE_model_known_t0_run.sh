#!/bin/bash

# redefine path for a tmp directory (the default tmp directory is mounted as non-executable on the new cluster)
export TMP=$TMP_DIR && export TMPDIR=$TMP && export TEMP=$TMP

echo 
echo output from ${SLURM_JOB_NAME}_run.sh starts here.
echo
echo TMP: $TMP

cd $WORK_DIR

echo current path: $PWD

Rscript --vanilla R_code/stan_sampling_SDE_model_known_t0.R $@


echo 
echo output from ${SLURM_JOB_NAME}_run.sh ends here.
echo
