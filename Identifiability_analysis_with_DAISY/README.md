## Structural identifiability analysis using DAISY (Differential Algebra for Identifiability of SYstems)

The results in this folder where generated using the `DAISY` software which can be downloaded from https://daisy.dei.unipd.it where also instructions for installing `REDUCE` and a tutorial how to use `DAISY` are provided. Here we use the  PSL version of `Reduce` which provides a command line interface.

Once `DAISY` and `Reduce` are installed:

* Run REDUCE PSL by opening the Terminal and entering the command: `path_to_REDUCE_folder/Reduce/psl/redpsl`
* Write the instruction: `load daisy$`

* To perform the identifability analysis you can then run either

  * `IN "path_to_DAISY_model_files\example_model.txt"$` 
  
    and the results will be printed to the command line; or
  * you first specify an output file by 
  
    `OUT "path_to_DAISY_model_files\example_model_result.txt"$` 
    
    `IN "path_to_DAISY_model_files\example_model.txt"$` 
    
    `SHUT "path_to_DAISY_model_files\example_model_result.txt"$` 
* Before starting a new identifiability study, run:

  `CLEAR ALL$`
  
* Exit `REDUCE` by writing

  `BYE$`  or   `QUIT$`
  
  
### Contained files
* `mRNA_transformed_model_1M.txt` 

    -- definition of the ODE model of the translation kinetics (1st order moments)
    
* `mRNA_transformed_model_1M_result.txt` 

    -- `DAISY` output for the ODE model
    
* `mRNA_transformed_model_2M.txt`

    -- definition of the surrogate model for the SDE model of the translation kinetics (2nd order moments)
    
* `mRNA_transformed_model_2M_result.txt` 
    
    -- `DAISY` output for the surrogate model of the SDE model

* `mRNA_transformed_model_reparametrized_1M.txt` 

    -- definition of a reparametrized ODE model of the translation kinetics (1st order moments)
    
* `mRNA_transformed_model_reparametrized_1M_result.txt` 

    -- `DAISY` output for the reparametrized ODE model showing that the degradation rates are locally identifiable
