Pre-processing of positive WES
================

- <a href="#introduction" id="toc-introduction">Introduction</a>
- <a href="#step-1-annotation-with-vep"
  id="toc-step-1-annotation-with-vep">Step 1: Annotation with VEP</a>
- <a href="#step-2-keep-rare-variants"
  id="toc-step-2-keep-rare-variants">Step 2: Keep rare variants</a>
- <a href="#step-3-transform-vcf-to-tsv"
  id="toc-step-3-transform-vcf-to-tsv">Step 3: Transform vcf to tsv</a>
- <a href="#step-4-annotation-with-cadd"
  id="toc-step-4-annotation-with-cadd">Step 4: Annotation with CADD</a>
- <a
  href="#step-5-merge-vep-and-cadd-files-to-get-rare-variants-annotated-by-cadd"
  id="toc-step-5-merge-vep-and-cadd-files-to-get-rare-variants-annotated-by-cadd">Step
  5: Merge VEP and CADD files to get rare variants annotated by CADD</a>
  - <a href="#load-data" id="toc-load-data">Load data</a>
  - <a href="#loop-on-all-files" id="toc-loop-on-all-files">Loop on all
    files</a>
  - <a href="#loop-on-all-files-without-removing-00-variants"
    id="toc-loop-on-all-files-without-removing-00-variants">Loop on all
    files without removing 0/0 variants</a>
  - <a href="#look-at-differences-between-the-two-loops"
    id="toc-look-at-differences-between-the-two-loops">Look at differences
    between the two loops</a>

Load libraries

``` r
library(xlsx)
library(dplyr)
```

# Introduction

We have a cohort composed of **28 patients** in which there are several
families. Whole Exome Sequencing (WES) has been performed and for all
patients, the responsible variant has been identified. We have “raw” vcf
files and we need to annotate them with VEP to keep only rare variants
and with CADD to run VIOLA.

Several steps has been executed on **genotoul** server. On **genotoul**
server, we can ask for the installation of almost all tools but we have
to use *sbatch* command to run a script or a tool.

# Step 1: Annotation with VEP

This step has been performed on **genotoul** server.

First, we add to create the *slurm* file.

``` bash
#!/bin/bash
#SBATCH -p workq
#SBATCH -t  01-00:00:00 #Acceptable time formats include "minutes", "minutes:seconds", "hours:minutes:seconds", "days-hours", "days-hours:minutes" and "days-hours:minutes:seconds". 
#SBATCH --mail-user=justine.labory@inrae.fr   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END

#Load modules
module load bioinfo/Ensembl-VEP/111.0

# no --species option because default is home_sapiens
# no --assembly option because it takes the higher available, and here it was GRCH38 at this time
# --cache_version specified because module load was for 110 version, and only 106 was installed at this time

vep -i $1
    -o $2
    --verbose
    --cache
    --dir_cache=/home/jlabory/save/tool/vep/
    --offline
    --assembly GRCh37
    --vcf
    --check_existing
    --af
    --af_1kg
    --af_gnomad
```

Then, we can run VEP.

``` bash
for file in $(ls /home/jlabory/save/data/positive_wes_vcf/*.vcf) ;
do new_name=$(echo $file | xargs -n1 basename  | cut -d '_' -f1) ;
sbatch slurm_loop_vep.sh $file "/home/jlabory/work/tool/vep/res/positive_wes/"$new_name"_output_vep.vcf" ;
done
```

# Step 2: Keep rare variants

This step has been performed on **genotoul** server.

First, we add to create the *slurm* file.

``` bash
#!/bin/bash
#SBATCH -p workq
#SBATCH -t  01-00:00:00 #Acceptable time formats include "minutes", "minutes:seconds", "hours:minutes:seconds", "days-hours", "days-hours:minutes" and "days-hours:minutes:seconds". 
#SBATCH --mail-user=justine.labory@inrae.fr   # email address
#SBATCH --mail-type=END

#Load modules
module load bioinfo/Ensembl-VEP/111.0


filter_vep -i $1 -o $2 -filter "SYMBOL and ((AF < 0.01 or gnomAD_AF < 0.01) or (not AF and not gnomAD_AF and not EUR_AF))"
```

Then, we can run *filter_vep*.

``` bash
for file in $(ls /home/jlabory/work/tool/vep/res/positive_wes/*.vcf) ; do new_name=$(echo $file | xargs -n1 basename | cut -d  '_' -f1) ;
sbatch slurm_loop_filter-vep.sh $file "/home/jlabory/work/tool/filter_vep/res/positive_wes/"$new_name"_rare_variant.vcf" ;
done
```

# Step 3: Transform vcf to tsv

Since, we will use *bcftools* to transform vcf files into tsv files,
first, we have to transform 1,2 and 3 for chromosomes into chr1, chr2
and chr3…

``` bash
for file in $(ls /home/jlabory/work/tool/filter_vep/res/positive_wes/*.vcf) ;
do new_name=$(echo $file | xargs -n1 basename | cut -d  '_' -f1) ;
awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' $file > "/home/jlabory/work/tool/filter_vep/res/positive_wes/"$new_name"_rare_variant_with_chr.vcf" ;
done
```

Then, we have to create a *slurm* file.

``` bash
#!/bin/bash
#SBATCH -p workq
#SBATCH -t  01-00:00:00 #Acceptable time formats include "minutes", "minutes:seconds", "hours:minutes:seconds", "days-hours", "days-hours:minutes" and "days-hours:minutes:seconds". 
#SBATCH --mail-user=justine.labory@inrae.fr   # email address
#SBATCH --mail-type=END

#Load module
module load bioinfo/Bcftools/1.17

echo -e "CHROM\tPOS\tREF\tALT\tGT\t$(bcftools +split-vep -l $1 | cut -f 2 | tr '\n' '\t' | sed 's/\t$//')" > $2 ;
bcftools +split-vep -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\t%CSQ\n' -d -A tab $1 >> $2
```

And now, we can run *bcftools +split-vep*.

``` bash
for file in $(ls /home/jlabory/work/tool/filter_vep/res/positive_wes/*with_chr.vcf) ;
do new_name=$(echo $file | xargs -n1 basename | cut -d  '_' -f1) ;
sbatch slurm_loop_split-vep.sh $file "/home/jlabory/work/tool/bcftools/split_vep/res/positive_wes/"$new_name"_rare_variant.tsv" ;
done
```

# Step 4: Annotation with CADD

This step was executed on **mdlab server**.

# Step 5: Merge VEP and CADD files to get rare variants annotated by CADD

## Load data

``` r
path_vep = "/Users/justine_labory/Desktop/These_ANR/data/res_splitvep/positive_wes"
path_cadd = "/Users/justine_labory/Desktop/These_ANR/data/res_cadd/res_positive_wes"

# Load table of responsible gene information
df_resp_var <- read.xlsx("/Users/justine_labory/Desktop/These_ANR/positive_variant_WES_MD.xlsx", sheetIndex = 1)
df_resp_var <- df_resp_var[,1:4]
head(df_resp_var)
```

    ##   File.code   Gene Chromosome   variant
    ## 1    01-WES   UFM1         13  38923901
    ## 2    02-WES   UFM1         13  38923901
    ## 3    03-WES  EPCAM          2  47601109
    ## 4    03-WES SLC3A1          2  44547590
    ## 5    04-WES  SYNE1          6 152652562
    ## 6    05-WES  SYNE1          6 152652562

## Loop on all files

- Step 1: Keep rare variants from VEP file
- Step 2: Merge CADD and VEP file
- Step 3: Filter out 0/0 variants

But be careful, for patient 24, the *GT* (genotype) information is not
available so we can not filter with *GT* column so we have to remove
this step for the pre-processing of this file.

``` r
for (i in c(sprintf("%02d",1:6),10:27)) {
  
  id_patient = paste0(i,"-WES")
  print(id_patient)
  
  # Read VEP file
  pat_vep <- read.table(paste0(path_vep, "/",id_patient,"_rare_variant.tsv"), header = T, sep = "\t")
  
  # Filter out variants with AF > 0.01 in any population
  pat_vep_f <- pat_vep %>%
    mutate_at(vars(AF:gnomADe_SAS_AF), ~ ifelse(. == ".", 0, as.numeric(.))) %>%
    filter_at(vars(AF:gnomADe_SAS_AF), all_vars((. < 0.01)))
  
  # Modify CHROM column in order to have 1,2,3 instead of chr1, chr2, chr3.
  pat_vep_f$CHROM <- str_replace_all(pat_vep_f$CHROM, "chr","")
  # Create an unique id with chromosome number, variant position, the gene name and the transcript. 
  pat_vep_f$unique_id = paste(pat_vep_f$CHROM, pat_vep_f$POS, pat_vep_f$SYMBOL, pat_vep_f$Feature, sep = "_")
  
  
  # Read CADD file for the same patient
  pat_cadd <- read.delim(paste0(path_cadd, "/",id_patient,"_cadd_output.tsv"), skip = 1, header = T)
  
  # Create an unique id with chromosome number, variant position, the gene name and the transcript.
  pat_cadd$unique_id <- paste(pat_cadd$X.Chrom, pat_cadd$Pos, pat_cadd$GeneName, pat_cadd$FeatureID, sep = "_")
  
  
  df_resp_var_pat <- df_resp_var %>%
    filter(File.code == id_patient)
  
  if (nrow(df_resp_var_pat) == 1) {
    
    pattern_resp_var_1 = paste(df_resp_var_pat$Chromosome, df_resp_var_pat$variant, df_resp_var_pat$Gene, sep = "_")
    pattern_resp_var_2 = "NOTAPATTERN"
    
  } else {
    
    pattern_resp_var_1 = paste(df_resp_var_pat$Chromosome[1], df_resp_var_pat$variant[1], df_resp_var_pat$Gene[1], sep = "_")
    pattern_resp_var_2 = paste(df_resp_var_pat$Chromosome[2], df_resp_var_pat$variant[2], df_resp_var_pat$Gene[2], sep = "_")
    
  }
  
  
  
  
  df_merge <- pat_cadd %>%
    # Remove variants which are not associated to a gene
    filter(! is.na(GeneName)) %>%
    # Merge CADD file with VEP file to get only rare variants
    inner_join(pat_vep_f, by=c("unique_id" = "unique_id",
                                "X.Chrom" = "CHROM",
                                "Pos" = "POS",
                                "Ref"="REF" ,
                                "Alt" = "ALT",
                                "GeneName" = "SYMBOL",
                                "FeatureID" = "Feature")) %>%
    # Remove variants with a genotype 0/0 which means that they are homozygous for the reference allele
    filter(GT != "0/0") %>%
    # Put unique_id column at the first position
    relocate(unique_id) %>%
    mutate(responsible_variant = case_when(
      str_starts(unique_id, pattern_resp_var_1) ~ 1,
      str_starts(unique_id, pattern_resp_var_2) ~ 1,
      TRUE ~ 0
    ))
  
  print(table(df_merge$responsible_variant))
  
  write.csv(df_merge, paste0("/Users/justine_labory/Desktop/These_ANR/All_Projects/Projet_VIOLA/positive_wes/input_viola/",id_patient,"_input_viola.csv"), row.names = F)
}
```

## Loop on all files without removing 0/0 variants

- Step 1: Keep rare variants from VEP file
- Step 2: Merge CADD and VEP file

``` r
for (i in c(sprintf("%02d",1:6),10:27)) {
  
  id_patient = paste0(i,"-WES")
  print(id_patient)
  
  # Read VEP file
  pat_vep <- read.table(paste0(path_vep, "/",id_patient,"_rare_variant.tsv"), header = T, sep = "\t")
  
  # Filter out variants with AF > 0.01 in any population
  pat_vep_f <- pat_vep %>%
    mutate_at(vars(AF:gnomADe_SAS_AF), ~ ifelse(. == ".", 0, as.numeric(.))) %>%
    filter_at(vars(AF:gnomADe_SAS_AF), all_vars((. < 0.01)))
  
  # Modify CHROM column in order to have 1,2,3 instead of chr1, chr2, chr3.
  pat_vep_f$CHROM <- str_replace_all(pat_vep_f$CHROM, "chr","")
  # Create an unique id with chromosome number, variant position, the gene name and the transcript. 
  pat_vep_f$unique_id = paste(pat_vep_f$CHROM, pat_vep_f$POS, pat_vep_f$SYMBOL, pat_vep_f$Feature, sep = "_")
  
  
  # Read CADD file for the same patient
  pat_cadd <- read.delim(paste0(path_cadd, "/",id_patient,"_cadd_output.tsv"), skip = 1, header = T)
  
  # Create an unique id with chromosome number, variant position, the gene name and the transcript.
  pat_cadd$unique_id <- paste(pat_cadd$X.Chrom, pat_cadd$Pos, pat_cadd$GeneName, pat_cadd$FeatureID, sep = "_")
  
  
  df_resp_var_pat <- df_resp_var %>%
    filter(File.code == id_patient)
  
  if (nrow(df_resp_var_pat) == 1) {
    
    pattern_resp_var_1 = paste(df_resp_var_pat$Chromosome, df_resp_var_pat$variant, df_resp_var_pat$Gene, sep = "_")
    pattern_resp_var_2 = "NOTAPATTERN"
    
  } else {
    
    pattern_resp_var_1 = paste(df_resp_var_pat$Chromosome[1], df_resp_var_pat$variant[1], df_resp_var_pat$Gene[1], sep = "_")
    pattern_resp_var_2 = paste(df_resp_var_pat$Chromosome[2], df_resp_var_pat$variant[2], df_resp_var_pat$Gene[2], sep = "_")
    
  }
  
  
  
  
  df_merge <- pat_cadd %>%
    # Remove variants which are not associated to a gene
    filter(! is.na(GeneName)) %>%
    # Merge CADD file with VEP file to get only rare variants
    inner_join(pat_vep_f, by=c("unique_id" = "unique_id",
                                "X.Chrom" = "CHROM",
                                "Pos" = "POS",
                                "Ref"="REF" ,
                                "Alt" = "ALT",
                                "GeneName" = "SYMBOL",
                                "FeatureID" = "Feature")) %>%
    # Remove variants with a genotype 0/0 which means that they are homozygous for the reference allele
    # filter(GT != "0/0") %>%
    # Put unique_id column at the first position
    relocate(unique_id) %>%
    mutate(responsible_variant = case_when(
      str_starts(unique_id, pattern_resp_var_1) ~ 1,
      str_starts(unique_id, pattern_resp_var_2) ~ 1,
      TRUE ~ 0
    ))
  
  print(table(df_merge$responsible_variant))
  
  write.csv(df_merge, paste0("/Users/justine_labory/Desktop/These_ANR/All_Projects/Projet_VIOLA/positive_wes/input_viola_allvar/",id_patient,"_input_viola_allvar.csv"), row.names = F)
}
```

## Look at differences between the two loops

``` r
path1 <- "/Users/justine_labory/Desktop/These_ANR/All_Projects/Projet_VIOLA/positive_wes/input_viola"
path2 <- "/Users/justine_labory/Desktop/These_ANR/All_Projects/Projet_VIOLA/positive_wes/input_viola_allvar"

list1 <- list.files(path = path1, patter="csv", full.names = T)
list2 <- list.files(path = path2, patter="csv", full.names = T)

df_res <- as.data.frame(matrix(nrow = length(list1), ncol = 5))
names(df_res) <- c("id_patient", "nb_all_variant", "nb_variant_filtered", "nb_remove_var", "perc_remove_var")

for (i in 1:length(list1)) {
  
  df1 = read.csv(list1[i])
  df2 = read.csv(list2[i])
  
  df_res$id_patient[i] = substr(basename(list1[i]), 1, 6)
  df_res$nb_all_variant[i] = nrow(df2)
  df_res$nb_variant_filtered[i] = nrow(df1)
  df_res$nb_remove_var[i] = nrow(df2) - nrow(df1)
  df_res$perc_remove_var[i] = round(( (nrow(df2) - nrow(df1)) / nrow(df2) ) *100, 2)
  
}

df_res

# write.csv(df_res, "/Users/justine_labory/Desktop/These_ANR/All_Projects/Projet_VIOLA/positive_wes/table_preprocessing_input_viola_file.csv", row.names = F)
```

    ##    id_patient nb_all_variant nb_variant_filtered nb_remove_var perc_remove_var
    ## 1      01-WES          82881               37550         45331           54.69
    ## 2      02-WES           2909                2906             3            0.10
    ## 3      03-WES          82881               39353         43528           52.52
    ## 4      04-WES          82881               39552         43329           52.28
    ## 5      05-WES          82881               39703         43178           52.10
    ## 6      06-WES           2631                2631             0            0.00
    ## 7      10-WES           2500                2500             0            0.00
    ## 8      11-WES           2407                2407             0            0.00
    ## 9      12-WES           8530                8530             0            0.00
    ## 10     13-WES           8469                8469             0            0.00
    ## 11     14-WES           8049                8049             0            0.00
    ## 12     15-WES           8972                8972             0            0.00
    ## 13     16-WES           8332                8332             0            0.00
    ## 14     17-WES           8159                8159             0            0.00
    ## 15     18-WES           8062                8062             0            0.00
    ## 16     19-WES           8652                8652             0            0.00
    ## 17     20-WES           7311                7311             0            0.00
    ## 18     21-WES           7296                7296             0            0.00
    ## 19     22-WES           7285                7285             0            0.00
    ## 20     23-WES           7633                7633             0            0.00
    ## 21     24-WES           3054                3054             0            0.00
    ## 22     25-WES           3420                3417             3            0.09
    ## 23     26-WES           2919                2911             8            0.27
    ## 24     27-WES          18171               18055           116            0.64
