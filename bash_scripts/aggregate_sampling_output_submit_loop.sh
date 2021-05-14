#!/bin/bash

v_datasets=(experimental_data_eGFP experimental_data_d2eGFP
            simulated_data_dataset1_no_error simulated_data_dataset1_with_error 
            simulated_data_dataset2_no_error simulated_data_dataset2_with_error)
            
v_model_types=(SDE ODE)

for dataset in "${v_datasets[@]}"
do
    for model_type in "${v_model_types[@]}"
    do
    
    sbatch bash_scripts/aggregate_sampling_output_submit.sh $dataset $model_type
    
    done
done