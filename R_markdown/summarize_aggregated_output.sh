# R_markdown/summarize_aggregated_output.sh

dataset=(experimental_data_d2eGFP experimental_data_eGFP \
simulated_data_dataset1_with_error simulated_data_dataset2_with_error \
simulated_data_dataset1_no_error simulated_data_dataset2_no_error)

model_type=(ODE SDE)

for ds in "${dataset[@]}"
do
  :
  for mt in "${model_type[@]}"
  do
    :
    echo
    echo dataset: $ds 
    echo model_type: $mt
    Rscript --vanilla -e 'pars <- commandArgs(trailingOnly=TRUE);
      rmarkdown::render(input = "R_markdown/summarize_aggregated_output.Rmd",
      output_format = "pdf_document", knit_root_dir = getwd(),
      output_file = paste0("pdf_files/summarize_aggregated_output_", 
      pars[1], "_", pars[2], ".pdf"),
      params = list(dataset = pars[1], model_type = pars[2]))' $ds $mt
    done
done