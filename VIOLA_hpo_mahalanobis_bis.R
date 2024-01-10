#!/usr/bin/Rscript

library(dplyr)
library("data.table")
library(ontologyIndex)
library(tidyr)
library(stringr)
library(optparse)


# Import functions

jaccard <- function(a, b) {
  
  ###""" This function compute the jaccard index between 2 vectors.
  ### It takes in input 2 vectors a and b.
  ### And it returns a float that is the jaccard index.""" 
  
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}

jaccard_for_hpo <- function(hpo_patient,hpo_gene) {
  
  ###""" This function compute the jaccard index between 2 HPO terms of all combinations of HPO terms of a patient and HPO terms of a gene.
  ### It takes in input 2 vectors, one with hpo terms of the patient and the other one with hpo terms of a gene.
  ### And it returns a dataframe in which there are all the possible combinations of HPO with the corresponding jaccard index.""" 
  
  combinaison <- CJ(hpo_patient, hpo_gene, unique = TRUE)
  
  jaccard_index <- c()
  for (i in 1:nrow(combinaison)) {
    
    v=get_ancestors(hpo,combinaison$hpo_patient[i])
    u=get_ancestors(hpo,combinaison$hpo_gene[i])
    
    if (length(v) > 0 & length(v) > 0) {
      
      index <- jaccard(v, u)
      jaccard_index <- c(jaccard_index, index)
      
    } else {
      jaccard_index <- c(jaccard_index, NA)
    }
    
  }
  combinaison$jaccard_index <- jaccard_index
  
  return(combinaison)
  
}

###############################################################################################################################################################################

select_anomaly_variant_with_hpo <- function(dbscan_res, hpo_patient, jaccard_th) {
  
  ###""" This function sort the anomaly variants (output of VIOLA).
  ### It takes in input the output of VIOLA (dataframe), the hpo of the patient (vector) and the threshold of jaccad index (float value) to sort anomaly variants.
  ### And it returns 2 objects, the first one is a dataframe in which there are sorted anomaly variants and the second one is a dataframe of combinations of HPO terms that have a jaccard index superior to the fixed thereshold.""" 
  
  jaccard_th <- as.numeric(jaccard_th)
  
  #anomaly_variant = viola_res
  
  # Get gene name of these variants
  anomaly_unique_gene <- unique(dbscan_res$GeneName)
  # Remove empty strings of the vector
  anomaly_unique_gene <- anomaly_unique_gene[anomaly_unique_gene != ""]
  
  hpo_simili_score_2_gene <- as.data.frame(matrix(nrow=length(anomaly_unique_gene), ncol=2))
  names(hpo_simili_score_2_gene) <- c("GeneName","jaccard_index")
  hpo_simili_score_2_gene$GeneName <- anomaly_unique_gene
  
  for (i in 1:length(anomaly_unique_gene)) {
    
    if (anomaly_unique_gene[i] %in% hpo_term_md_related_to_a_gene$gene_symbol) {
      
      hpo_gene <- hpo_term_md_related_to_a_gene$hpo_id[hpo_term_md_related_to_a_gene$gene_symbol == anomaly_unique_gene[i]]
      
      simili_gene <- jaccard_for_hpo(hpo_patient, hpo_gene)
      hpo_simili_score_2_gene$jaccard_index[i] <- max(simili_gene$jaccard_index, na.rm = T)
      
    } else {
      hpo_simili_score_2_gene$jaccard_index[i] = 0
    }
    
  }
  
  # Keep only genes with a jaccard superior to the threshold given by the user
  filter_jaccard <- hpo_simili_score_2_gene %>%
    dplyr::filter(jaccard_index >= jaccard_th)
  
  if (nrow(filter_jaccard) > 0) {
    
    #anomaly_variant_with_jaccard_score <- merge(anomaly_variant, filter_jaccard, by="GeneName")
    anomaly_variant_with_jaccard_score <- dbscan_res %>%
      inner_join(filter_jaccard)

    
    return(anomaly_variant_with_jaccard_score)
    
  } else {
    
    return(paste0("No gene has a jaccard index > ", jaccard_th))
  }
  
}

##########################################
############## MAIN ######################
##########################################

#################
## ARGUMENTS ####
#################

# create an option parser
option_list <- list(
  make_option(c("-f", "--input"), dest="input_file",
              help="path to input file"),
  make_option(c("-o", "--output_path"), dest="output_path",
              help="path to output directory"),
  make_option(c("-t", "--hpo_table"), dest="hpo_table_path",
              help="path to table containing patient hpo terms")
)
parser <- OptionParser(option_list=option_list)

# parse the command line arguments
args <- parse_args(parser)

# use the parsed arguments
input_file <- args$input_file
output_path <- args$output_path
hpo_table_path <- args$hpo_table_path

# Load databases
hpo_term_md_related_to_a_gene <- read.delim("/home/justine/hpo_db/hpoterms_md_to_gene.txt")

# Load HPO database
hpo <- get_ontology("https://raw.githubusercontent.com/obophenotype/human-phenotype-ontology/master/hp.obo",
                    extract_tags = "everything")

hpo_tibble <- simplify2array(hpo) %>% #Convert to array
  as_tibble()

# Load patient hpo terms
#table_hpo <- read.csv("/home/justine/projet_VIOLA/clinical_positive_exome_hpoterms_for_viola.csv")
table_hpo <- read.csv(hpo_table_path)

#################
## HPO FILTER ###
#################

message(date(), " : Apply HPO filter")

viola_res = read.csv(input_file)
file_prefix <- strsplit(basename(input_file), '_')[[1]][1]
hpo_patient <- table_hpo$hpo_id[table_hpo$patient_id == file_prefix]

res_HPO = select_anomaly_variant_with_hpo(viola_res, hpo_patient, 0.8)
write.csv(res_HPO, file = paste0(output_path,file_prefix,'_res_HPO.csv'), row.names = FALSE)

########################
## MAHALANOBIS SCORE ###
########################
  
message(date(), " : Compute mahalanobis score")

# Features hg19
features = c('ConsScore', 'SIFTval', 'PolyPhenVal', 'priPhCons', 'mamPhCons', 'verPhCons', 'priPhyloP', 'mamPhyloP',
             'verPhyloP', 'EncH3K27Ac', 'EncH3K4Me1', 'EncH3K4Me3', 'EncExp', 'EncNucleo', 'EncOCC', 'EncOCCombPVal',
             'EncOCDNasePVal', 'EncOCFairePVal', 'EncOCpolIIPVal', 'EncOCctcfPVal', 'EncOCmycPVal', 'EncOCDNaseSig',
             'EncOCFaireSig', 'EncOCpolIISig', 'EncOCctcfSig', 'EncOCmycSig', 'SpliceAI.acc.gain', 'SpliceAI.acc.loss',
             'SpliceAI.don.gain', 'SpliceAI.don.loss', 'MMSp_acceptorIntron', 'MMSp_acceptor', 'MMSp_exon', 'MMSp_donor',
             'MMSp_donorIntron','Freq100bp', 'Rare100bp', 'Sngl100bp', 'Freq1000bp', 'Rare1000bp', 'Sngl1000bp', 'Freq10000bp',
             'Rare10000bp','Sngl10000bp', 'RawScore' , 'PHRED')


# # Feature hg38
# features = c('priPhCons', 'mamPhCons', 'verPhCons', 'priPhyloP', 'mamPhyloP',
#              'verPhyloP', 'EncodeH3K27ac.max', 'EncodeH3K27ac.sum','EncodeH3K4me1.max', 'EncodeH3K4me1.sum','EncodeH3K4me3.sum',
#              'EncodeH3K4me3.max', 'EncodeDNase.sum','EncodeDNase.max', 'MMSp_acceptorIntron', 'MMSp_acceptor', 'MMSp_exon', 'MMSp_donor',
#              'MMSp_donorIntron', 'Freq100bp', 'Rare100bp', 'Sngl100bp', 'Freq1000bp', 'Rare1000bp', 'Sngl1000bp', 'Freq10000bp',
#              'Rare10000bp', 'Sngl10000bp', 'RawScore', 'PHRED')



data <- res_HPO[, features]
#data <- data %>% mutate(across(all_of(features), ~replace_na(.,round(median(., na.rm=TRUE)))))
data[features] <- sapply(data[features],as.double)
data <- data %>%
  mutate(across(all_of(features),
                ~replace_na(., median(., na.rm=TRUE))))


data$mahalnobis<- mahalanobis(data, colMeans(data), cov(data), tol=1e-25)
data$pvalue <- pchisq(data$mahalnobis, df=3, lower.tail=FALSE)
outliers <- data[data$pvalue < 0.001, "pvalue", drop = FALSE]
#dim(outliers)

res = cbind(res_HPO, data[, c("mahalnobis", "pvalue")])


# Select rows with pvalue < 0.001, Cluster = -1, and responsible_variant = 1
subset_res <- res[res$pvalue < 0.001, ]
subset_res <- subset_res %>%
  arrange(pvalue)

subset_res$rank <- rank(subset_res$pvalue, ties.method = "min")

write.csv(subset_res, file = paste0(output_path,file_prefix,'_mahalanobis_res.csv'), row.names = FALSE)

