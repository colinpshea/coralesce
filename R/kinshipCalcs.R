#' Calculate pairwise kinship for all possible pairwise comparisons.
#' @description This function calculates pairwise kinship based on all alleles present at all possible combinations of colonies and loci.
#' @param dataset A data frame with `Coral_ID` column as unique identifier for `1:nrow(dataset)` colonies and SNP data at `2:ncol(dataset)` loci; this is a data frame resulting from `convertBasePairstoCodes()`.
#' @returns This function calculates pairwise kinship for all possible comparisons of alleles among individuals. Note that these are simply kinship calculations and this function does *not* exclude invariant loci (`kinshipCalcsNoInvar()` uses `omitInvariantLoci()`and `kinshipCalcs()`). 
#' @importFrom matrixStats rowProds
#' @export
kinshipCalcs <- function(dataset){
  dataset$Aprob1 <- 0  
  dataset$Cprob1 <- 0  
  dataset$Gprob1 <- 0  
  dataset$Tprob1 <- 0  
  dataset$Aprob1 <- ifelse(dataset$allele1=="A", 2/2, dataset$Aprob1)
  dataset$Aprob1 <- ifelse(dataset$allele1=="M"|dataset$allele1=="R"|dataset$allele1=="W", 1/2, dataset$Aprob1)
  dataset$Cprob1 <- ifelse(dataset$allele1=="C", 2/2, dataset$Cprob1)
  dataset$Cprob1 <- ifelse(dataset$allele1=="M"|dataset$allele1=="S"|dataset$allele1=="Y", 1/2, dataset$Cprob1)
  dataset$Gprob1 <- ifelse(dataset$allele1=="G", 2/2, dataset$Gprob1)
  dataset$Gprob1 <- ifelse(dataset$allele1=="R"|dataset$allele1=="S"|dataset$allele1=="K", 1/2, dataset$Gprob1)
  dataset$Tprob1 <- ifelse(dataset$allele1=="T", 2/2, dataset$Tprob1)
  dataset$Tprob1 <- ifelse(dataset$allele1=="W"|dataset$allele1=="Y"|dataset$allele1=="K", 1/2, dataset$Tprob1)
  dataset$Aprob2 <- 0  
  dataset$Cprob2 <- 0  
  dataset$Gprob2 <- 0  
  dataset$Tprob2 <- 0  
  dataset$Aprob2 <- ifelse(dataset$allele2=="A", 2/2, dataset$Aprob2)
  dataset$Aprob2 <- ifelse(dataset$allele2=="M"|dataset$allele2=="R"|dataset$allele2=="W", 1/2, dataset$Aprob2)
  dataset$Cprob2 <- ifelse(dataset$allele2=="C", 2/2, dataset$Cprob2)
  dataset$Cprob2 <- ifelse(dataset$allele2=="M"|dataset$allele2=="S"|dataset$allele2=="Y", 1/2, dataset$Cprob2)
  dataset$Gprob2 <- ifelse(dataset$allele2=="G", 2/2, dataset$Gprob2)
  dataset$Gprob2 <- ifelse(dataset$allele2=="R"|dataset$allele2=="S"|dataset$allele2=="K", 1/2, dataset$Gprob2)
  dataset$Tprob2 <- ifelse(dataset$allele2=="T", 2/2, dataset$Tprob2)
  dataset$Tprob2 <- ifelse(dataset$allele2=="W"|dataset$allele2=="Y"|dataset$allele2=="K", 1/2, dataset$Tprob2)
  dataset$NumAlleles <- 4 - 2*rowSums(is.na(dataset[, c("allele1", "allele2")]))
  dataset$mult_probA <- rowProds(as.matrix(dataset[, c("Aprob1","Aprob2")]))
  dataset$mult_probC <- rowProds(as.matrix(dataset[, c("Cprob1","Cprob2")]))
  dataset$mult_probG <- rowProds(as.matrix(dataset[, c("Gprob1","Gprob2")]))
  dataset$mult_probT <- rowProds(as.matrix(dataset[, c("Tprob1","Tprob2")]))
  
  dataset %>% filter(NumAlleles==4) %>% mutate(combo = interaction(coral1, coral2)) %>% group_by(combo) %>% mutate(totalProbLoci = sum(mult_probA, mult_probC, mult_probG, mult_probT), totalScorable = n()) %>% select(combo, totalProbLoci, totalScorable) %>% mutate(avg_kinship = totalProbLoci/totalScorable) %>% separate(combo, into = c("coral1", "coral2"), remove = T, sep = "[.]") %>% distinct(.)
} 
