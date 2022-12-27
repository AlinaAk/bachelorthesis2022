#This is code to replicate the analyses and figures from the Bachelorthesis
#of Alina Akerlund in 2022. Code developed by Alina Akerlund, some functions
#co-authored by Markus Gumbel

#### Libraries ####

#install.packages("dplyr")
library(dplyr)
#install.packages("tibble")
library(tibble)
#install.packages("seqinr")
library(seqinr)
#install.packages("/home/mgumbel/devel/gcatbase_0.6.0_R_x86_64-pc-linux-gnu.tar.gz")
library(gcatbase)
#install.packages("rlist")
library(rlist)
#install.packages("devtools")
library(devtools)
#install_github("VariantEffect/hgvsParseR")
library(hgvsParseR)
#install.packages("ggplot2")
library(ggplot2)
#install.packages("tidyr")
library(tidyr)
#install.packages("patchwork")
library(patchwork)

#### Functions ####

#' Function to filter cosmic_cds for all CDS variants of relevant genes
#' @param fasta List of CDS-sequences
#' @param gene_names Data frame with column GENE_NAME
#' @return List with all CDS variants of the searched genes
cosmic_genes <- function(fasta, gene_names) {
  new_list <- list()
  for(i in 1:nrow(gene_names)) {
    elements <- fasta[grep(gene_names[i, "GENE_NAME"], names(fasta))]
    new_list <- append(new_list, elements)
  }
  return(new_list)
}

#' Function to map a tissue to another
#' @param data Data frame in which the tissue should be edited
#' @param tissue Name of the tissue to be edited
#' @param new_tissue New name of the tissue
#' @return Data frame with edited tissue names
map_tissue <- function(data, tissue, new_tissue) {
  data$Primary.site[data$Primary.site == tissue] <- new_tissue
  return(data)
}

#' Function to calculate the mean of TPM values for each gene per each tissue
#' @param data Data frame in which the new values should be saved
#' @param data_tissues Data frame with column Primary.site and all values
#' found in the first parameter
#' @param data_genes Data frame with column GENE_NAME and all values
#' found in the first parameter
#' @return Data frame with newly calculated TPM values
mean_tpm <- function(data, data_tissues, data_genes) {
  #filter for one gene at a time
  for(i in 1:nrow(data_genes)) {
    print(i)
    temp_data_gene <-  data %>% filter(GENE_NAME %in% 
      data_genes[i, "GENE_NAME"])

    #go through tissues and see if there is more than one TPM value
    for(j in 1:nrow(data_tissues)) {
      temp_data_tissue <- temp_data_gene %>% filter(Primary.site %in% 
        data_tissues[j, "Primary.site"])

      #calculate mean if there is more than one value
      if(nrow(temp_data_tissue) > 1) {
        tpm_mean <- mean(temp_data_tissue$TPM)

        #change old different TPM values with new mean TPM value
        data[data$GENE_NAME == data_genes[i, "GENE_NAME"] & 
          data$Primary.site == data_tissues[j, "Primary.site"], "TPM"] <- 
          tpm_mean
      }
    }
  }
  return(data)
}

#' Function to calculate the mean of Z-score values for each gene
#' per each tissue
#' @param data Data frame in which the new values should be saved
#' @param data_tissues Data frame with column Primary.site and all values
#' found in the first parameter
#' @param data_genes Data frame with column GENE_NAME and all values
#' found in the first parameter
#' @return Data frame with newly calculated Z-score values
mean_z_score <- function(data, data_tissues, data_genes) {
  #filter for one gene at a time
  for(i in 1:nrow(data_genes)) {
    print(i)
    temp_data_gene <-  data %>% filter(GENE_NAME %in% 
      data_genes[i, "GENE_NAME"])

    #go through tissues and see if there is more than one Z-score value
    for(j in 1:nrow(data_tissues)) {
      temp_data_tissue <- temp_data_gene %>% filter(Primary.site %in% 
        data_tissues[j, "Primary.site"])

      #calculate mean if there is more than one value
      if(nrow(temp_data_tissue) > 1) {
        mean_z <- mean(temp_data_tissue$Z_SCORE)

        #change old different Z-score values with new mean Z-score value
        data[data$GENE_NAME == data_genes[i, "GENE_NAME"] &
          data$Primary.site == data_tissues[j, "Primary.site"], "Z_SCORE"] <-
          mean_z
      }
    }
  }
  return(data)
}

#' Function to create a table for median and sd part of Z-score calculation
#' @param data Data frame with all TPM-values
#' @param data_tissues Data frame with column Primary.site and all values
#' found in the first parameter
#' @return Data frame with columns: Primary.site, MEDIAN_TPM and SD_TPM
get_tpm_table <- function(data, data_tissues) {
  #new data frame
  tab <- data.frame(Primary.site = character(), MEDIAN_TPM = double(),
    SD_TPM = double())
  
  #calculate median- and sd-part of Z-Score calculation for each tissue
  for (i in 1:nrow(data_tissues)) {
    print(i)
    temp_tissue <- data %>% filter(Primary.site %in% 
      data_tissues[i, "Primary.site"]) %>% distinct()
    temp_median <- median(log2((temp_tissue$TPM) + 1))
    temp_sd <- sd(log2((temp_tissue$TPM) + 1))
    tab[nrow(tab) + 1, ] <- c(data_tissues[i, "Primary.site"], temp_median, temp_sd)
  }
  tab <- tab %>% mutate(MEDIAN_TPM = as.numeric(MEDIAN_TPM),
    SD_TPM = as.numeric(SD_TPM))
  return(tab)
}

#' Function to calculate Z-score with the TPM-value (median and sd 
#' over one tissue)
#' @param data Data frame with all TPM values
#' @param median_sd_table Data frame with columns GENE_NAME, MEDIAN_TPM
#' and SD_TPM
#' @return Data frame with calculated Z-score values in the new column
#' Z_SCORE
get_zscore_values <- function(data, median_sd_table) {
  #add column to data frame
  data <- data %>% add_column(Z_SCORE = "")
  data$Z_SCORE <- as.numeric(data$Z_SCORE)

  #calculate and save the Z-Score in data$Z_SCORE
  for(i in 1:nrow(tpm_table)) {
    print(i)
    data$Z_SCORE[data$Primary.site == tpm_table[i, "Primary.site"]] <-
      (log2(data$TPM[data$Primary.site == tpm_table[i, "Primary.site"]] + 1) -
      tpm_table[i, "MEDIAN_TPM"]) / tpm_table[i, "SD_TPM"]
  }
  return(data)
}

#' Function to weight the scaled Z-score on a scale from 0-1
#' @param data Data frame with the column SCALED_Z_SCORE
#' @return Data frame with new column WEIGHTED_Z_SCORE
weighted_z_score <- function(data) {
  #add column to data frame
  data <- data %>% add_column(WEIGHTED_Z_SCORE = "")
  data$WEIGHTED_Z_SCORE <- as.numeric(data$WEIGHTED_Z_SCORE)

  #calculate the weighted Z-score for each scaled Z-score
  for(i in 1:nrow(data)) {
    print(i)
    temp_w <- 1 / (1 + exp(-data[i, "SCALED_Z_SCORE"]))
    data[i, "WEIGHTED_Z_SCORE"] <- temp_w
  }
  return(data)
}

#' Function to get CDS-length and weighting of the length in one table
#' @param census_genes_cds List with all sequences
#' @param gene_names Data frame with the column GENE_NAMES
cosmic_genes_length <- function(census_genes_cds, gene_names) {
  #create a new data frame
  tab <- data.frame(GENE_NAME = character(), CDS_LENGTH = numeric(),
    WEIGHTED_CDS_LENGTH = double())
  
  #count the length for each gene's CDS and save it in the new data frame
  for (i in 1:nrow(gene_names)) {
    print(i)
    gene <- census_genes_cds[[gene_names[i, "GENE_NAME"]]]
    seq <- getSequence(gene, as.string = TRUE)
    seq_length <- nchar(seq)
    tab[nrow(tab) + 1, ] <- c(gene_names[i, "GENE_NAME"], seq_length/3, NA)
  }
  tab <- tab %>% mutate(CDS_LENGTH = as.numeric(CDS_LENGTH),
    WEIGHTED_CDS_LENGTH = as.numeric(WEIGHTED_CDS_LENGTH))
  
  #calculate the sum of all lengths together
  sum_cds <- sum(tab$CDS_LENGTH)

  #calculate the weighting of each gene's CDS-length
  for(j in 1:nrow(tab)) {
    print(j)
    w_length <- tab[j, "CDS_LENGTH"] / sum_cds
    tab[j, "WEIGHTED_CDS_LENGTH"] <- w_length
  }
  return(tab)
}

#' Function to parse HGVSC code of each row in exp_mutants
#' @param data Data frame with column HGVSC
#' @param fasta List of all sequences
#' @param circ_code Circular code to search for in the sequences
#' @return Data frame with new column CODE_USAGE
get_mutation_code_usage <- function(data, fasta, circ_code) {
  #add column to data frame
  data <- data %>% add_column(CODE_USAGE = "")
  
  #iterate through data and parse each HGVSC code
  for(j in 1:nrow(data)) {
    print(j)
    temp_data <- data[j, ]
    result <- parseHGVS(temp_data$Mutation.CDS)
    #combine temp_data and result
    temp_data <- cbind(temp_data, result[1, ])
    #mutate cds depending on the type of mutation
    if(temp_data[1, "type"] == "insertion" |
      temp_data[1, "type"] == "singledeletion" |
      temp_data[1, "type"] == "deletion") {
        seq <- mutated_indels(temp_data[1, ], fasta)
    } else if(temp_data[1, "type"] == "substitution") {
      seq <- mutated_codon(temp_data[1, ], fasta)
    } else if(temp_data[1, "type"] == "invalid") {
      seq <- NULL
      cu <- NULL
    }

    #calc code usage of mutated sequence and add to data
    if(!is.null(seq)) {
      seq_norm <- normalize(seq, lowercase = FALSE)
      cu <- code_usage(seq_norm, circ_code)
    } else if(is.null(seq)) {
      cu <- NULL
    }
      
    if(!is.null(cu)) {
      data[j, "CODE_USAGE"] <- cu
    } else if(is.null(cu)) {
      data[j, "CODE_USAGE"] <- NA
    }
  }
  return(data)
}

#' Mutated codons.
#' @param mut_entry A mutation entry (from COSMIC table)
#' @param cds_ref Reference CDS.
#' @return List with three attributes:
#' "from", "to", "length"
mutated_codon <- function(mut_entry, cds_ref) {
  cds_name <- mut_entry$GENE_NAME
  cds <- normalize(getSequence(cds_ref[[cds_name]],
    as.string = TRUE)) # reference sequence

  if (mut_entry$type != "substitution") {
    stop("Internal error: not a substitution!")
  }
  start <- mut_entry$start
  from_base <- mut_entry$ancestral
  to_base <- mut_entry$variant
  cds_mut <- modify_subst(cds, start, from_base, to_base)
  cds_mut
}

#' Modify a string sequence.
#' @param seq The sequence as a string.
#' @param idx The index position of the base to modify.
#' @param from_base Old base
#' @param to_base New base
#' @return Modified sequence
modify_subst <- function(seq, idx, from_base, to_base) {
  if (is.na(idx)) {
    print("Invalid index.")
    return(NULL)
  }
  s <- strsplit(seq, "")[[1]]

  if (idx > length(s)) {
    warning(paste0("Mutation index ", idx, " greater than sequence length ",
      length(s)))
    NULL
  } else if (s[idx] != from_base) {
    env <- substr(seq, idx - 1, idx + 1)
    warning(paste0(
      "Position ", idx, ": ancestral base ", s[idx], " should be ", from_base,
      " (range ", env, ")"
    ))
    NULL
    # stop("Illegal substitution.")
  } else {
    s[idx] <- to_base
    paste(s, collapse = "")
  }
}

#' New sequence after insert or deletions.
#'
#' `mut_entry` contains the attributes `start`, `end` and `variant`
#' which are relevant here. The values are according to:
#' http://varnomen.hgvs.org/ or
#' den Dunnen, J.T. et al. (2016) �HGVS Recommendations for the Description of Sequence Variants: 2016 Update�,
#' Human Mutation, 37(6), pp. 564�569. doi:10.1002/humu.22981.
#' `start` and `end` range from 1 to n where n is the
#' number of nucleotides in `cds_ref`. Their meaning depends on the kind of
#' modification:
#' _insert_
#' The new nucleotides are in `variant` and are inserted
#' _after_ the `start` position. I.e. there is one exception:
#' if nucleotides should be inserted at the beginning of the sequence, then
#' `start` is 0.
#' Explanation from http://varnomen.hgvs.org/recommendations/RNA/variant/insertion/
#' positions_flanking: position two nucleotides flanking insertion site, e.g. 123_124
#' _singledeletion_
#' Only the `start`-element is used. It represents the position of the
#' nucleotide to be deleted.
#' _deletion_
#' `start` and `end` define the range of the positions of the nucleotides
#' to be deleted.
#' @param mut_entry A mutation entry (from COSMIC table)
#' @param cds_ref Reference CDS.
#' @return Altered sequence
mutated_indels <- function(mut_entry, cds_ref) {
  cds_name <- mut_entry$GENE_NAME
  cds <- getSequence(cds_ref[[cds_name]], as.string = TRUE) # reference sequence

  if (!(mut_entry$type == "insertion" |
    mut_entry$type == "singledeletion" |
    mut_entry$type == "deletion")) {
    stop("Internal error: not an indels!")
  }
  # Note: substr(string, from , to). to is inclusive.
  mod_cds <- if (mut_entry$type == "insertion") {
    after_pos <- mut_entry$start
    cds_insert <- mut_entry$variant
    # Build the new sequence:
    cds_b <- substr(cds, 1, after_pos) # begin to mutation
    cds_e <- substr(cds, after_pos + 1, nchar(cds)) # mutation to end
    paste0(cds_b, cds_insert, cds_e)
  } else if (mut_entry$type == "singledeletion") {
    pos <- mut_entry$start
    # Build the new sequence:
    cds_b <- substr(cds, 1, pos - 1) # begin to mutation
    cds_e <- substr(cds, pos + 1, nchar(cds)) # mutation to end
    paste0(cds_b, cds_e)
  } else { # deletion of two or more nucs.
    start <- mut_entry$start
    end <- mut_entry$end
    # Build the new sequence:
    cds_b <- substr(cds, 1, start - 1) # begin to mutation
    cds_e <- substr(cds, end + 1, nchar(cds)) # mutation to end
    paste0(cds_b, cds_e)
  }

  mod_cds
}

#' Function to calculate the code usage of codons
#' relatively to the other codons
#' @param seq Sequence to be analyzed
#' @param x Set of codons to calculate the code usage for
code_usage <- function(seq, x) {
 h = classify(seq, X)
 u = length(h[h == 1]) / length(h)
} 
  
#' Function to calculate code usage per gene
#' @param data Data frame with the column GENE_NAMES
#' @param fasta List of all sequences
#' @param gene_names Data frame with the column GENE_NAMES
#' @param circ_code Circular code to search for in the sequences
#' @return Data frame with new column CODE_USAGE
get_code_usage <- function(data, fasta, gene_names, circ_code) {
  #new data frame
  tab <- data.frame(GENE_NAME = character(), CODE_USAGE = double())
  #calculate the code usage for each gene
  for (i in 1:nrow(gene_names)) {
    print(i)
    gene <- fasta[[gene_names[i, "GENE_NAME"]]]
    seq <- getSequence(gene, as.string = TRUE)
    seq_norm <- normalize(seq, lowercase = FALSE)
    cu <- code_usage(seq_norm, circ_code)
    tab[nrow(tab) + 1, ] <- c(gene_names[i, "GENE_NAME"], cu)
  }
  tab <- tab %>% mutate(CODE_USAGE = as.numeric(CODE_USAGE))

  #save the code usage for each gene in the main data frame
  data <- data %>% inner_join(tab, by = "GENE_NAME")
  return(data)
}

#' Function to calc the mean of code usage for one gene in one tissue
#' for carcinogenic sequences
#' @param data Data frame with column CODE_USAGE
#' @param data_tissues Data frame with column Primary.site and all values
#' found in the first parameter
#' @param data_genes Data frame with column GENE_NAME and all values
#' found in the first parameter
mean_code_usage <- function(data, data_tissues, data_genes) {
  for(i in 1:nrow(data_genes)) {
    print(i)
    temp_data_gene <-  data %>% filter(GENE_NAME %in% 
      data_genes[i, "GENE_NAME"])

    #go through tissues and see if there is more than one Z-score value
    for(j in 1:nrow(data_tissues)) {
      temp_data_tissue <- temp_data_gene %>% filter(Primary.site %in% 
        data_tissues[j, "Primary.site"])

      #calc mean if there is more than one value
      if(nrow(temp_data_tissue) > 1) {
        mean_cu <- mean(temp_data_tissue$CODE_USAGE)
        #change old different Z-score values with new mean Z-score value
        data[data$GENE_NAME == data_genes[i, "GENE_NAME"] &
          data$Primary.site == data_tissues[j, "Primary.site"], "CODE_USAGE"] <-
          mean_cu
      }
    }
  }
  return(data)
}

#' Function to calculate the weighted arith mean for each code usage per
#' tissue
#' @param data Data frame with column GENE_NAME, Primary.site,
#' WEIGHTED_Z_SCORE, WEIGHTED_CDS_LENGTH
#' @param data_tissues Data frame with column Primary.site and all values
#' found in the first parameter
#' @return Data frame with columns Primary.site, CODE_USAGE and
#' WEIGHTED_CODE_USAGE
weighted_code_usage_per_tissue <- function(data, data_tissues) {
  #new data frame
  tab <- data.frame(Primary.site = character(),
    WEIGHTED_CODE_USAGE = double(), CODE_USAGE = double())
  
  #calculate the mean for each tissue and save it in tab
  for (i in 1:nrow(data_tissues)) {
    print(i)
    temp_tissue <- data %>% filter(Primary.site ==
      data_tissues[i, "Primary.site"]) %>% distinct()
    upper_sum <- sum(temp_tissue$WEIGHTED_Z_SCORE *
      temp_tissue$WEIGHTED_CDS_LENGTH *
      temp_tissue$CODE_USAGE, na.rm = TRUE)
    lower_sum <- sum(temp_tissue$WEIGHTED_Z_SCORE *
      temp_tissue$WEIGHTED_CDS_LENGTH, na.rm = TRUE)
    weighted_cu <- upper_sum/lower_sum
    tab[nrow(tab) + 1, ] <- c(data_tissues[i, "Primary.site"], weighted_cu,
      mean(temp_tissue$CODE_USAGE, na.rm = TRUE))
  }
  tab <- tab %>% mutate(WEIGHTED_CODE_USAGE = as.double(WEIGHTED_CODE_USAGE),
    CODE_USAGE = as.double(CODE_USAGE))
  return(tab)
}

#### Chapter 3.2.2: Reading the data ####

## mutants ##
mutants <- read.table(file = "~/home2/CosmicMutantExportCensus.tsv.gz", #change location if needed
                      header = TRUE,
                      sep = "\t",
                      quote = "",
                      dec = ".",
                      blank.lines.skip = FALSE)
saveRDS(mutants, file = "~/home2/mutants.rds") #change location if needed

## expression ##
expression <- read.table(file = '~/home2/CosmicCompleteGeneExpression.tsv.gz', #change location if needed
                         header = TRUE,
                         sep = "",
                         quote = "",
                         dec = ".",
                         blank.lines.skip = FALSE)
saveRDS(expression, file = "~/home2/expression.rds") #change location if needed

## hpa_data ##
hpa_data <- read.table(file = "~/home2/rna_tissue_hpa.tsv", #change location if needed
                      header = TRUE,
                      sep = "\t",
                      quote = "",
                      dec = ".",
                      blank.lines.skip = FALSE)
saveRDS(hpa_data, file = "~/home2/hpa_data.rds") #change location if needed

## cosmic_cds ##
cosmic_cds <- read.fasta(file = "~/home2/All_COSMIC_Genes.fasta.gz") #change location if needed
saveRDS(cosmic_cds, file = "~/home2/cosmic_cds.rds") #change location if needed

#### Chapter 3.2.3: Preparation of the data ####

#read the objects faster if they aren't already in your workspace
expression <- readRDS("~/home2/expression.rds")
mutants <- readRDS("~/home2/mutants.rds")
hpa_data <- readRDS("~/home2/hpa_data.rds")
cosmic_cds <- readRDS("~/home2/cosmic_cds.rds")

### Rename and delete columns ###
## mutants ##
mutants <- mutants %>% rename(GENE_NAME = Gene.name, SAMPLE_NAME = Sample.name,
  SAMPLE_ID = ID_sample) %>% select(-c("Gene.CDS.length"))

## hpa_data ##
hpa_data <- hpa_data %>% rename(GENE_NAME = Gene.name,
  Primary.site = Tissue) %>% select(-c("pTPM", "nTPM"))

### Filter for Census Genes ###
#list of census genes
mutants_genes <- mutants %>% select(GENE_NAME) %>% distinct()

## expression ##
expression <- expression %>% filter(GENE_NAME %in% mutants_genes$GENE_NAME)

## hpa_data ##
hpa_data <- hpa_data %>% filter(GENE_NAME %in% mutants_genes$GENE_NAME)

## cosmic_cds ##
cosmic_cds <- cosmic_genes(cosmic_cds, mutants_genes)

### Sort out irrelevant data ###
mutants <- subset(mutants, Mutation.Description != "Unknown")

### Mapping the tissues of mutants and hpa_data ###
## mutants ##
mutants <- map_tissue(mutants, "thyroid", "thyroid_gland")
mutants <- map_tissue(mutants, "prostate", "male_genital_tract")
mutants <- map_tissue(mutants, "meninges", "central_nervous_system")
mutants <- map_tissue(mutants, "parathyroid", "parathyroid_gland")
mutants <- map_tissue(mutants, "penis", "male_genital_tract")
mutants <- map_tissue(mutants, "pituitary", "central_nervous_system")
mutants <- map_tissue(mutants, "pleura", "lung")
mutants <- map_tissue(mutants, "perineum", "perineum")
mutants <- map_tissue(mutants, "autonomic_ganglia", "autonomic_nervous_system")
mutants <- map_tissue(mutants, "testis", "male_genital_tract")
mutants <- map_tissue(mutants, "placenta", "female_genital_tract")
mutants <- map_tissue(mutants, "vulva", "female_genital_tract")
mutants <- map_tissue(mutants, "female_genital_tract_(site_indeterminate)",
  "female_genital_tract")
mutants <- map_tissue(mutants, "vagina", "female_genital_tract")
mutants <- map_tissue(mutants, "uterine_adnexa", "female_genital_tract")
mutants <- map_tissue(mutants, "fallopian_tube", "female_genital_tract")
mutants <- map_tissue(mutants, "pericardium", "heart")
mutants <- map_tissue(mutants, "paratesticular_tissues", "male_genital_tract")
mutants <- map_tissue(mutants, "retroperitoneum", "stomach")

## hpa_data ##
hpa_data <- map_tissue(hpa_data, "adipose tissue", "soft_tissue")
hpa_data <- map_tissue(hpa_data, "adrenal gland", "adrenal_gland")
hpa_data <- map_tissue(hpa_data, "amygdala", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "angular gyrus", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "angular gyrus (white matter)", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "anterior cingulate cortex, supragenual-dorsal", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "anterior cingulate cortex, supragenual-ventral", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "anterior cingulate gyrus, pregenual-dorsal", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "anterior cingulate gyrus, pregenual-ventral", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "anterior cingulate gyrus, subgenual", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "anterior cochlear nucleus, ventral", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "anterior funiculus (white matter)", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "anterior insular cortex, dorsal", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "anterior insular cortex, ventral", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "anterior thalamic nucleus, ventral", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "appendix", "haematopoietic_and_lymphoid_tissue")
hpa_data <- map_tissue(hpa_data, "arcuate nucleus", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "area parastriata, inferior", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "area parastriata, parietal", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "area parastriata, superior", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "area parastriata, temporal", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "area postrema", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "area striata", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "area striata (white matter)", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "basal amygdala", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "bone marrow", "bone")
hpa_data <- map_tissue(hpa_data, "caudate nucleus", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "central amygdala", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "centromedial thalamic nucleus", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "cerebellar cortex", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "cerebellar nuclei", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "cerebellar white matter", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "cervical spinal cord, central gray", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "cervical spinal cord, dorsal horn", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "cervical spinal cord, ventral horn", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "cervix", "female_genital_tract")
hpa_data <- map_tissue(hpa_data, "choroid plexus", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "claustrum", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "colon", "large_intestine")
hpa_data <- map_tissue(hpa_data, "corpus callosum, genu (white matter)", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "corpus callosum, splenium (white matter)", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "corticomedial amygdala", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "cuneate nuclei", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "cuneiform nucleus", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "dentate gyrus", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "dorsal cochlear nucleus", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "dorsal funiculus (white matter)", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "dorsal medullary reticular nucleus", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "dorsal motor vagal complex", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "dorsal raphe nucleus", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "dorsal tegmental nucleus", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "dorsolateral prefrontal cortex", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "dorsolateral tegmental area", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "dorsomedial nucleus", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "dorsomedial prefrontal cortex", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "duodenum", "small_intestine")
hpa_data <- map_tissue(hpa_data, "entorhinal gyrus", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "epididymis", "male_genital_tract")
hpa_data <- map_tissue(hpa_data, "esophagus", "oesophagus")
hpa_data <- map_tissue(hpa_data, "fallopian tube", "female_genital_tract")
hpa_data <- map_tissue(hpa_data, "flocculonodular lobe", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "frontal cortex", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "frontal eye field", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "frontal operculum", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "frontomarginal cortex", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "frontopolar cortex", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "fusiform gyrus", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "gallbladder", "biliary_tract")
hpa_data <- map_tissue(hpa_data, "gigantocellular reticular nuclei, medullary", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "gigantocellular reticular nuclei, pars", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "globus pallidus, externus", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "globus pallidus, internus", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "gyrus rectus", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "habenula", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "heart muscle", "heart")
hpa_data <- map_tissue(hpa_data, "hippocampus", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "hippocampus, CA1", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "hippocampus, CA2", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "hippocampus, CA3", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "inferior colliculus", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "inferior frontal gyrus, opercular", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "inferior frontal gyrus, orbital", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "inferior frontal gyrus, triangular", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "inferior olive", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "inferior parietal lobule", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "inferior temporal gyrus", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "intraparietal deep sulcus", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "Kolliker-Fuse nucleus", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "lateral amygdala", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "lateral funiculus (white matter)", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "lateral geniculate body", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "lateral hypothalamic area", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "lateral lemniscal nuclei", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "lateral medullary reticular nucleus", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "lateral parabrachial nucleus", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "lateral thalamic nuclei", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "lateral vestibular nucleus", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "lingual gyrus", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "locus coeruleus", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "lymph node", "haematopoietic_and_lymphoid_tissue")
hpa_data <- map_tissue(hpa_data, "mamillary body", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "medial dorsal thalamic nucleus", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "medial geniculate body", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "medial olivary nucleus", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "medial parabrachial nucleus", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "medial periolivary nuclei", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "medial vestibular nucleus", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "median raphe nucleus", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "middle cingulate cortex", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "middle frontal gyrus", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "middle temporal gyrus", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "motor facial nucleus", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "motor hypoglossal nucleus", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "motor trigeminal nucleus", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "nuclei of the trapezoid body", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "nucleus accumbens", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "nucleus ambiguus", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "nucleus basalis of Meynert", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "nucleus cuneatus", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "nucleus gracilis", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "nucleus of the diagonal band", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "nucleus raphe magnus", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "nucleus raphe obscurus", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "nucleus raphe pallidus", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "nucleus reuniens", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "nucleus rhomboideus", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "nucleus tractus solitarii", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "occipital cortex", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "occipital cortex (white matter)", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "olfactory area", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "olfactory tubercle", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "orbitofrontal gyrus, anterior", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "orbitofrontal gyrus, lateral", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "orbitofrontal gyrus, medial", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "orbitofrontal gyrus, posterior", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "paracentral lobule, anterior", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "paracentral lobule, posterior", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "parahippocampal cortex", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "parahippocampal cortex (white matter)", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "paramedian reticular nucleus", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "parathyroid gland", "parathyroid_gland")
hpa_data <- map_tissue(hpa_data, "paraventricular nucleus", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "parietal cortex", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "parietal operculum", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "parieto-insular cortex", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "parieto-occipital transitional area", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "parieto-temporal junction (white matter)", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "parvicellular reticular nuclei", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "parvicellular reticular nuclei, medullary", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "pedunculopontine tegmental nucleus", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "periaquaductal grey, anterior", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "periaquaductal grey, dorsal", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "periaquaductal grey, lateral", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "perirhinal gyrus", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "peritrigeminal nucleus", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "piriform cortex", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "placenta", "female_genital_tract")
hpa_data <- map_tissue(hpa_data, "pontine nuclei", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "pontine raphe nucleus", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "postcentral gyrus", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "postcentral gyrus, dorsal", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "postcentral gyrus, middle", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "postcentral gyrus, ventral", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "posterior cingulate cortex", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "posterior cingulate cortex, dorsal", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "posterior cingulate cortex, ventral", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "posterior insular cortex", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "posterior thalamic nucleus", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "posteroventral cochlear nucleus", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "precentral gyrus", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "precentral gyrus, dorsal", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "precentral gyrus, middle", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "precentral gyrus, ventral", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "precuneus", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "premotor cortex, dorsal", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "premotor cortex, ventral", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "preoptic area", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "prepositus hypoglossal nucleus", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "pretectal area", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "principal sensory trigeminal nucleus", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "principal sensory trigeminal nucleus", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "prostate", "male_genital_tract")
hpa_data <- map_tissue(hpa_data, "pulvinar", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "putamen", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "rectum", "large_intestine")
hpa_data <- map_tissue(hpa_data, "red nucleus", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "reticular pontine nucleus, caudal", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "reticular pontine nucleus, oral", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "reticulotegmental nucleus", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "retrosplenial cortex", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "salivary gland", "salivary_gland")
hpa_data <- map_tissue(hpa_data, "seminal vesicle", "male_genital_tract")
hpa_data <- map_tissue(hpa_data, "septal nuclei", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "skeletal muscle", "skeletal_muscle")
hpa_data <- map_tissue(hpa_data, "small intestine", "small_intestine")
hpa_data <- map_tissue(hpa_data, "smooth muscle", "smooth_muscle")
hpa_data <- map_tissue(hpa_data, "somatomor cortex, precentral gyrus (white matter)", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "somatosensory cortex, postcentral gyrus, ventral (white matter)", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "spinal trigeminal nucleus, caudal", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "spinal trigeminal nucleus, interpolar", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "spinal trigeminal nucleus, oral", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "spinal vestibular nucleus", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "spleen", "haematopoietic_and_lymphoid_tissue")
hpa_data <- map_tissue(hpa_data, "stria terminalis, bed nucleus", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "subcallosal gyrus", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "subcentral gyrus, S2", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "subcoeruleus area", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "subiculum", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "substantia nigra", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "subthalamic nucleus", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "superior colliculus", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "superior frontal gyrus", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "superior olive", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "superior parietal lobule", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "superior temporal gyrus", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "superior vestibular nucleus", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "supplementary motor cortex", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "supramarginal gyrus", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "supraoptic nucleus", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "temporal cortex", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "temporal pole", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "temporal white matter", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "temporo-insular cortex, parainsular gyrus", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "temporo-occipital transitional zone", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "testis", "male_genital_tract")
hpa_data <- map_tissue(hpa_data, "thyroid gland", "thyroid_gland")
hpa_data <- map_tissue(hpa_data, "tongue", "upper_aerodigestive_tract")
hpa_data <- map_tissue(hpa_data, "tonsil", "haematopoietic_and_lymphoid_tissue")
hpa_data <- map_tissue(hpa_data, "transversal temporal gyrus, anterior", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "transversal temporal gyrus, posterior", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "urinary bladder", "urinary_tract")
hpa_data <- map_tissue(hpa_data, "ventral medullary reticular nucleus", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "ventral periolivary nuclei", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "ventral posterolateral thalamic nucleus", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "ventral posteromedial thalamic nucleus", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "ventral tegmental area", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "ventral thalamic nuclei", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "ventrolateral medulla, A1-C1 cell groups", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "ventrolateral prefrontal cortex", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "ventrolateral tegmental area, A5 NE cell group", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "ventromedial nucleus", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "ventromedial prefrontal cortex", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "vermis", "central_nervous_system")
hpa_data <- map_tissue(hpa_data, "zona incerta", "central_nervous_system")

### Join mutants and expression ###
exp_mutants <- expression %>%
  inner_join(mutants, by = c("SAMPLE_ID" = "SAMPLE_ID",
  "SAMPLE_NAME" = "SAMPLE_NAME", "GENE_NAME" = "GENE_NAME",
  "ID_STUDY" = "ID_STUDY"))
  #%>% select(-c("Site.subtype.1", "Site.subtype.2", "Site.subtype.3",
  #"Primary.histology", "Histology.subtype.1",
  #"Histology.subtype.2", "Histology.subtype.3",
  #"Genome.wide.screen", "Mutation.zygosity",
  #"LOH", "GRCh", "Mutation.strand",
  #"Resistance.Mutation", "Pubmed_PMID",
  #"Sample.Type", "Tumour.origin", "Age"))

#### Chapter 3.2.4: Calculations ####

### Calculate gene expressionlevel mean for each gene per each tissue ###
## hpa_data ##
#parameters
hpa_data_tissues <- hpa_data %>% select(Primary.site) %>% distinct()
hpa_data_genes <- hpa_data %>% select(GENE_NAME) %>% distinct()
#calculation
hpa_data <- mean_tpm(hpa_data, hpa_data_tissues, hpa_data_genes)
#delete duplicates
hpa_data <- hpa_data %>% distinct()

## exp_mutants ##
#parameters
exp_mutants_tissues <- exp_mutants %>% select(Primary.site) %>% distinct()
exp_mutants_genes <- exp_mutants %>% select(GENE_NAME) %>% distinct()
#calculation
exp_mutants <- mean_z_score(exp_mutants, exp_mutants_tissues,
  exp_mutants_genes)


### Calculate Z-scores with TPM values ###
#parameter
hpa_data_tissues <- hpa_data %>% select(Primary.site) %>% distinct()
#first, get a table with the median- and sd-part of the calculation
tpm_table <- get_tpm_table(hpa_data, hpa_data_tissues)
#cut out all TPM values that are 0, because it can't be projected on Z-scores
hpa_data <- hpa_data %>% filter(TPM != 0)
#calculate Z-scores
hpa_data <- get_zscore_values(hpa_data, tpm_table)


### Scale Z-scores ###
## exp_mutants ##
exp_mutants <- exp_mutants %>% mutate(SCALED_Z_SCORE = scale(Z_SCORE))
#relocate SCALED_Z_SCORE after Z_SCORE
exp_mutants <- exp_mutants %>% relocate(SCALED_Z_SCORE, .after = Z_SCORE)

## hpa_data ##
hpa_data <- hpa_data %>% mutate(SCALED_Z_SCORE = scale(Z_SCORE))


### Calculate weighted Z-scores ###
## exp_mutants ##
exp_mutants <- weighted_z_score(exp_mutants)
#relocate weighted Z-score after scaled Z-score
exp_mutants <- exp_mutants %>% relocate(WEIGHTED_Z_SCORE,
  .after = SCALED_Z_SCORE)

## hpa_data ##
hpa_data <- weighted_z_score(hpa_data)


### Calculate CDS-length and its weight ###
## exp_mutants ##
#parameter
exp_mutants_genes <- exp_mutants %>% select(GENE_NAME) %>% distinct()
#get for every gene: GENE_NAME, CDS_LENGTH, WEIGHTED_CDS_LENGTH
census_cds_length <- cosmic_genes_length(cosmic_cds, exp_mutants_genes)
#join CDS data with exp_mutants
exp_mutants <- exp_mutants %>% inner_join(census_cds_length, by = "GENE_NAME")
#relocate the CDS-length after the weighted Z-score and the weighted CDS-length
#after the CDS-length
exp_mutants <- exp_mutants %>% relocate(CDS_LENGTH,
  .after = WEIGHTED_Z_SCORE) %>% relocate(WEIGHTED_CDS_LENGTH,
  .after = CDS_LENGTH)

## hpa_data ##
#parameter
hpa_data_genes <- hpa_data %>% select(GENE_NAME) %>% distinct()
#get for every gene: GENE_NAME, CDS_LENGTH, WEIGHTED_CDS_LENGTH
census_cds_length <- cosmic_genes_length(cosmic_cds, hpa_data_genes)
#join CDS data with hpa_data
hpa_data <- hpa_data %>% inner_join(census_cds_length, by = "GENE_NAME")


### Calculate the code usage of x_0 ###
#parameter
x_0 <- code(tuples = c("AAC", "AAT", "ACC", "ATC", "ATT", "CAG", "CTC", "CTG",
  "GAA", "GAC", "GAG", "GAT", "GCC", "GGC", "GGT", "GTA", "GTC", "GTT",
  "TAC", "TTC"))
## exp_mutants ##
#calculate code usage for all mutations
exp_mutants <- get_mutation_code_usage(exp_mutants, cosmic_cds, x_0)
#relocate the code usage after the weighted CDS-length
exp_mutants <- exp_mutants %>% relocate(CODE_USAGE,
  .after = WEIGHTED_CDS_LENGTH) %>%
  mutate(CODE_USAGE = as.numeric(CODE_USAGE))

## hpa_data ##
hpa_data <- get_code_usage(hpa_data, cosmic_cds,
  mutants_genes, x_0)


### Calculate the mean of the code usage for each gene per tissue ###
## exp_mutants ##
#parameter
exp_mutants_tissues <- exp_mutants %>% select(Primary.site) %>% distinct()
exp_mutants_genes <- exp_mutants %>% select(GENE_NAME) %>% distinct()
#calculation
exp_mutants <- mean_code_usage(exp_mutants, exp_mutants_tissues,
  exp_mutants_genes)
#subset exp_mutants so that there is only one row per gene per tissue
exp_mutants <- exp_mutants %>% select(GENE_NAME, Primary.site, Z_SCORE,
  SCALED_Z_SCORE, WEIGHTED_Z_SCORE, CDS_LENGTH, WEIGHTED_CDS_LENGTH,
  CODE_USAGE) %>% distinct()


### Calculate the weighted code usage per tissue ###
weighted_cu_table_carc <- weighted_code_usage_per_tissue(exp_mutants,
  exp_mutants_tissues)


#### Graphics/plots ####

exp_mut_z_plot <- ggplot(exp_mutants, aes(Z_SCORE)) +
  geom_histogram(bins = 500) +
  scale_x_continuous(limits = c(-5, 5)) +
  scale_y_continuous(limits = c(0, 10500)) +
  geom_vline(xintercept = mean(exp_mutants$Z_SCORE), color = "red") +
  labs(x = "Z_SCORE", y = "Anzahl") +
  theme_gray(base_size = 16) +
  theme(axis.text.x = element_text(color = "#3a3a3a", size = 16),
    axis.text.y = element_text(color = "#3a3a3a", size = 16),
    plot.title = element_text(hjust = 0.5), ) +
  geom_text(aes(mean(exp_mutants$Z_SCORE), label = "0.1273344733",
    10000), color = "red", hjust = -0.1, size = 5)

exp_mut_scaled_z_plot <- ggplot(exp_mutants, aes(SCALED_Z_SCORE)) +
  geom_histogram(bins = 500) +
  scale_x_continuous(limits = c(-5, 5)) +
  scale_y_continuous(limits = c(0, 10500)) +
  geom_vline(xintercept = mean(exp_mutants$SCALED_Z_SCORE), color = "red") +
  labs(x = "SCALED_Z_SCORE", y = "Anzahl") +
  theme_gray(base_size = 16) +
  theme(axis.text.x = element_text(color = "#3a3a3a", size = 16),
    axis.text.y = element_text(color = "#3a3a3a", size = 16),
    plot.title = element_text(hjust = 0.5, size = 18)) +
  geom_text(aes(mean(exp_mutants$SCALED_Z_SCORE), label = "0",
    10000), color = "red", hjust = -1.5, size = 5)

hpa_data_z_plot <- ggplot(new_hpa_data, aes(Z_SCORE)) +
  geom_histogram(bins = 500) +
  scale_x_continuous(limits = c(-5, 5)) +
  scale_y_continuous(limits = c(0, 250)) +
  geom_vline(xintercept = mean(new_hpa_data$Z_SCORE), color = "red") +
  labs(x = "Z_SCORE", y = "Anzahl") +
  theme_gray(base_size = 16) +
  theme(axis.text.x = element_text(color = "#3a3a3a", size = 16),
    axis.text.y = element_text(color = "#3a3a3a", size = 16),
    plot.title = element_text(hjust = 0.5)) +
  geom_text(aes(mean(new_hpa_data$Z_SCORE), label = "0.033684001",
    250), color = "red", hjust = -0.1, size = 5)

hpa_data_scaled_z_plot <- ggplot(new_hpa_data, aes(SCALED_Z_SCORE)) +
  geom_histogram(bins = 500) +
  scale_x_continuous(limits = c(-5, 5)) +
  scale_y_continuous(limits = c(0, 250)) +
  geom_vline(xintercept = mean(new_hpa_data$SCALED_Z_SCORE), color = "red") +
  labs(x = "SCALED_Z_SCORE", y = "Anzahl") +
  theme_gray(base_size = 16) +
  theme(axis.text.x = element_text(color = "#3a3a3a", size = 16),
    axis.text.y = element_text(color = "#3a3a3a", size = 16),
    plot.title = element_text(hjust = 0.5, size = 18)) +
  geom_text(aes(mean(new_hpa_data$SCALED_Z_SCORE), label = "0",
    250), color = "red", hjust = -1.5, size = 5)


z_plots <- exp_mut_z_plot + exp_mut_scaled_z_plot + 
  hpa_data_z_plot + hpa_data_scaled_z_plot +
  plot_annotation(
    title = "Verteilung der Z-Scores in exp_mutants (oben) und hpa_data (unten)
      vor (links) und nach der Skalierung (rechts)",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 18))
  )

ggsave(filename = "z_plots.png",
  plot = z_plots,
  path = "~/home2")



hpa_weighted_z_plot <- ggplot(new_hpa_data, aes(x = SCALED_Z_SCORE,
  y = WEIGHTED_Z_SCORE)) +
  geom_vline(xintercept = 0, color = "lightgrey") +
  geom_hline(yintercept = 0.5, color = "lightgrey") +
  geom_line(colour = "red") +
  theme_classic() +
  theme(axis.text.x = element_text(color = "#3a3a3a", size = 16),
    axis.text.y = element_text(color = "#3a3a3a", size = 16))

exp_mut_weighted_z_plot <- ggplot(exp_mutants, aes(x = SCALED_Z_SCORE,
  y = WEIGHTED_Z_SCORE)) +
  geom_vline(xintercept = 0, color = "lightgrey") +
  geom_hline(yintercept = 0.5, color = "lightgrey") +
  geom_line(colour = "red") +
  theme_classic() +
  theme(axis.text.x = element_text(color = "#3a3a3a", size = 16),
    axis.text.y = element_text(color = "#3a3a3a", size = 16))

weighted_z_plots <- exp_mut_weighted_z_plot + hpa_weighted_z_plot +
  plot_annotation(
    title = "Verteilung der gewichteten Z-Scores in exp_mutants (links) und hpa_data (rechts)",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 18))
  )

ggsave(filename = "weighted_z_plots.png",
  plot = weighted_z_plots,
  path = "~/home2", width = 12)



mean_cu_healthy <- 0.4286
mean_cu_carc <- 0.4283

#prepare data to be plotted
same <- inner_join(weighted_cu_table_healthy, weighted_cu_table_carc,
  by = "Primary.site")
same <- same %>% rename(WEIGHTED_CODE_USAGE_HEALTHY = WEIGHTED_CODE_USAGE.x,
    CODE_USAGE_HEALTHY = CODE_USAGE.x,
    WEIGHTED_CODE_USAGE_CARC = WEIGHTED_CODE_USAGE.y,
    CODE_USAGE_CARC = CODE_USAGE.y)

sub_same <- same %>% select(!c(WEIGHTED_CODE_USAGE_HEALTHY,
  WEIGHTED_CODE_USAGE_CARC))
sub_same_tall <- sub_same %>% pivot_longer(cols =
  c(CODE_USAGE_HEALTHY, CODE_USAGE_CARC),
  names_to = "Gewebeart", values_to = "CODE_USAGE",
  names_prefix = "CODE_USAGE_")

sub_same_tall$Gewebeart[sub_same_tall$Gewebeart == "HEALTHY"] <- "normal"
sub_same_tall$Gewebeart[sub_same_tall$Gewebeart == "CARC"] <- "karzinogen"


code_usage_plot <- ggplot() +
  geom_col(data = sub_same_tall, aes(x = CODE_USAGE, y = Primary.site,
    fill = Gewebeart), position = "dodge") +
  geom_vline(xintercept = mean_cu_healthy,
    color = "black", linewidth = 0.5, linetype = "longdash") +
  geom_vline(xintercept = mean_cu_carc,
    color = "red", linewidth = 0.5, linetype = "longdash") +
  scale_fill_manual(values = c("#0F3277", "#A1D0E4")) +
  theme_classic() +
  scale_x_continuous(breaks = c(0.0, 0.1, 0.2, 0.3, 0.4,
    mean_cu_healthy)) +
  theme(axis.text.x = element_text(angle = 55, hjust = 1,
    color = "#3a3a3a"),
    axis.text.y = element_text(color = "#3a3a3a", size = 14),
    title = element_text(size = 16),
    legend.text = element_text(size = 14)) +
    labs(title = "Die ungewichtete Code Usage pro Gewebe",
      x = "Code Usage", y = "Gewebe")

ggsave(filename = "code_usage_plot.png",
  plot = code_usage_plot,
  path = "~/home2")


mean_w_cu_healthy <- 0.4231
mean_w_cu_carc <- 0.3605

w_sub_same <- same %>% select(!c(CODE_USAGE_HEALTHY, CODE_USAGE_CARC))
w_sub_same_tall <- w_sub_same %>% pivot_longer(cols =
  c(WEIGHTED_CODE_USAGE_HEALTHY, WEIGHTED_CODE_USAGE_CARC),
  names_to = "Gewebeart", values_to = "WEIGHTED_CODE_USAGE",
  names_prefix = "WEIGHTED_CODE_USAGE_")

w_sub_same_tall$Gewebeart[w_sub_same_tall$Gewebeart == "HEALTHY"] <- "normal"
w_sub_same_tall$Gewebeart[w_sub_same_tall$Gewebeart == "CARC"] <- "karzinogen"


weighted_code_usage_plot <- ggplot() +
  geom_col(data = w_sub_same_tall, aes(x = WEIGHTED_CODE_USAGE,
    y = Primary.site, fill = Gewebeart), position = "dodge") +
  geom_vline(xintercept = mean_w_cu_healthy,
    color = "black", linewidth = 0.5, linetype = "longdash") +
  geom_vline(xintercept = mean_w_cu_carc,
    color = "red", linewidth = 0.5, linetype = "longdash") +
  scale_fill_manual(values = c("#0F3277", "#A1D0E4")) +
  theme_classic() +
  scale_x_continuous(breaks = c(0.0, 0.1, 0.2, 0.3, 0.4,
    mean_w_cu_carc, mean_w_cu_healthy)) +
  theme(axis.text.x = element_text(angle = 55, hjust = 1,
    color = "#3a3a3a"),
    axis.text.y = element_text(color = "#3a3a3a", size = 14),
    title = element_text(size = 16),
    legend.text = element_text(size = 14)) +
  labs(title = "Die gewichtete Code Usage pro Gewebe",
    x = "gewichtete Code Usage", y = "Gewebe")

ggsave(filename = "weighted_code_usage_plot.png",
  plot = weighted_code_usage_plot,
  path = "~/home2")


exp_mut_z_scores_table <- exp_mutants %>% select(Primary.site, GENE_NAME,
  SCALED_Z_SCORE) %>% mutate(SOURCE = "exp_mutants")
hpa_z_score_table <- new_hpa_data %>% select(Primary.site, GENE_NAME,
  SCALED_Z_SCORE) %>% mutate(SOURCE = "hpa_data")
z_scores_table <- rbind(exp_mut_z_scores_table, hpa_z_score_table)
