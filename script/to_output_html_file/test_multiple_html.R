options(java.parameters = c("-XX:+UseConcMarkSweepGC", "-Xmx8192m"))
gc()

library(stringr)
library(purrr)

##########################################################################################################################################################################################################################
# Output HTML files for files ranked by VIOLA score

input_path = "/Users/justine_labory/Desktop/github/VIOLA/test_filter_after_VAE/res_viola/raw_vcf"
list_files <- list.files(input_path, pattern = "_viola_score.csv", full.names = T, recursive = T)
list_files <- list_files[-21]


for (i in 1:length(list_files)) {

  id_patient <-  str_split(basename(list_files[i]), "_")[[1]] [1]

  rmarkdown::render(
    input  = "/Users/justine_labory/Desktop/These_ANR/All_Projects/Projet_VIOLA/script_to_output_html_file/version2_080125/template_viola_report.Rmd",
    params = list(
      directory = paste0(input_path, "/", id_patient, "/"),
      file      = basename(list_files[i])
    ),
    ## To output all the variants of the files
    #output_file = paste0("/Users/justine_labory/Desktop/github/VIOLA/test_filter_after_VAE/res_viola/raw_vcf/viola_report_080125/viola_score/", id_patient,"_VS_report.html")
    
    ## To output only the top 20 variants of the files
    output_file = paste0("/Users/justine_labory/Desktop/github/VIOLA/test_filter_after_VAE/res_viola/raw_vcf/viola_report_080125/viola_score_top20/", id_patient,"_VS_report.html")
  )
}

## Test on 1 file to debug
# rmarkdown::render(
#   input  = "/Users/justine_labory/Desktop/These_ANR/All_Projects/Projet_VIOLA/script_to_output_html_file/version2_080125/template_viola_report.Rmd",
#   params = list(
#     directory = "/Users/justine_labory/Desktop/github/VIOLA/test_filter_after_VAE/res_viola/raw_vcf/1-VT/",
#     file      = "1-VT_viola_score.csv"
#   ),
#   output_file = paste0("/Users/justine_labory/Desktop/github/VIOLA/test_filter_after_VAE/res_viola/raw_vcf/viola_report_080125/1-VT_test_report.html")
# )

##########################################################################################################################################################################################################################
# Output HTML files for files ranked by VIOLA score genotype

input_path = "/Users/justine_labory/Desktop/github/VIOLA/test_filter_after_VAE/res_viola/raw_vcf"
list_files <- list.files(input_path, pattern = "viola_score_genotype.csv", full.names = T, recursive = T)
list_files <- list_files[-21]

for (i in 1:length(list_files)) {
  
  id_patient <-  str_split(basename(list_files[i]), "_")[[1]] [1]
  
  rmarkdown::render(
    input  = "/Users/justine_labory/Desktop/These_ANR/All_Projects/Projet_VIOLA/script_to_output_html_file/version2_080125/template_viola_report_VS_genotype.Rmd",
    params = list(
      directory = paste0(input_path, "/", id_patient, "/"),
      file      = basename(list_files[i])
    ),
    ## To output all the variants of the files
    # output_file = paste0("/Users/justine_labory/Desktop/github/VIOLA/test_filter_after_VAE/res_viola/raw_vcf/viola_report_080125/viola_score_genotype/", id_patient,"_VS_genotype_report.html")
    
    ## To output only the top 20 variants of the files
    output_file = paste0("/Users/justine_labory/Desktop/github/VIOLA/test_filter_after_VAE/res_viola/raw_vcf/viola_report_080125/viola_score_genotype_top20/", id_patient,"_VS_genotype_report.html")
  )
}

## Test on 1 file to debug
# rmarkdown::render(
#   input  = "/Users/justine_labory/Desktop/These_ANR/All_Projects/Projet_VIOLA/script_to_output_html_file/version2_080125/template_viola_report_VS_genotype.Rmd",
#   params = list(
#     directory = "/Users/justine_labory/Desktop/github/VIOLA/test_filter_after_VAE/res_viola/raw_vcf/1-VT/",
#     file      = "1-VT_viola_score_genotype.csv"
#   ),
#   output_file = paste0("/Users/justine_labory/Desktop/github/VIOLA/test_filter_after_VAE/res_viola/raw_vcf/viola_report_080125/1-VT_test_report.html")
# )
