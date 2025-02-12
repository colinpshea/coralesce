#' Calculate population- and individual-level kinship (and gene diversity)
#'
#' @description This function calculates population- and individual-level mean kinship and gene diversity (1 - kinship) across all individuals and loci. These calculations exclude invariant loci and are based on all possible among-locus pairwise comparisons at each of `N` loci. Within-individual comparisons (i.e., comparisons of individuals with themselves) are also excluded from these calculations. Optionally, individual-level kinship can be calculated for a subset of individuals using the `targetN` argument.  
#' @param dataset A data frame with `Coral_ID` column as unique identifier (rows) and loci as columns; this is a data frame resulting from the  `convertBasePairstoCodes()` function.
#' @param subset Do you want to subset the data and calculate kinship for the `targetN` least related colonies in the data set? The default value is FALSE, in which case `targetN` is ignored and defaults to NULL. If `subset = TRUE` then you MUST enter a value for `targetN`; otherwise the function will quit and print an error message.
#' @param targetN The desired number of individuals over which population-average kinship and gene diversity are calculated. This is ignored if left as `NULL` or if its value is greater than or equal to `nrow(dataset)` i.e., the number of individuals in a data set. If this value is < `nrow(dataset)`, then kinship is recalculated repeatedly, removing the individual with the highest average kinship, one-by-one, until `targetN` individuals with the lowest average kinship remain.
#' @returns This function returns three objects, `PopAvgMKGD`, `MK_init`, and `MK_final`. `PopAvgMKGD` contains the population-level average kinship and gene diversity; `MK_init` contains average kinship for each individual in dataset (i.e., `now(dataset)` individuals), and `MK_final` contains average kinship for the `targetN` individuals with the lowest average kinship. 
#' 
#' @importFrom matrixStats rowProds
#' @export
kinshipCalcsNoInvar <- function(dataset, subset = FALSE, targetN = NULL){
  if (subset == FALSE){
    if (is.null(targetN)==FALSE) {stop(cat(paste("subset = FALSE but you have entered a value for targetN."), paste("You either mistakenly entered a value for targetN or meant to specify subset = TRUE."), sep = "\n"))
      }
  dat1 <- omitInvariantLoci(dataset = dataset)
  dat2 <- determineAllAlleleMatchesOthers(dataset = dat1)
  dat3 <- kinshipCalcs(dataset = dat2)
  #### Calculate Population-level mean kinship and GD
  PopAvgMKGD <- dat3 %>% pivot_longer(cols = c(coral1, coral2), values_to = "Coral_ID") %>% select(Coral_ID, avg_kinship) %>% arrange(Coral_ID, desc(avg_kinship)) %>% group_by(Coral_ID) %>% summarize(ind_mean_kinship = mean(avg_kinship)) %>% arrange(desc(ind_mean_kinship)) %>% ungroup() %>% summarise(PopAvgMK = mean(ind_mean_kinship), PopAvgGD = 1 - PopAvgMK)
  
  #### Calculate individual-level mean kinship
  MK_init <- dat3 %>% pivot_longer(cols = c(coral1, coral2), values_to = "Coral_ID") %>% select(Coral_ID, avg_kinship) %>% arrange(Coral_ID, desc(avg_kinship)) %>% group_by(Coral_ID) %>% summarize(ind_mean_kinship = mean(avg_kinship)) %>% arrange(desc(ind_mean_kinship))

return(list(PopAvgMKGD = PopAvgMKGD, MK_init = MK_init, MK_final = NULL))
    }  
  if (subset==TRUE){
    if (is.null(targetN)==TRUE) {stop(cat(paste("subset = TRUE but you have not entered a value for targetN."), paste("When subset = TRUE, a targetN value must be entered for this function to work properly."), paste("You either meant to enter subset = FALSE or forgot to enter a value for targetN."), sep = "\n"))
    }
    
  if (targetN <2) {stop(cat("When subset = TRUE, targetN must be ≥ 2", sep = "\n"))
    }

    dat1 <- omitInvariantLoci(dataset = dataset)
    dat2 <- determineAllAlleleMatchesOthers(dataset = dat1)
    dat3 <- kinshipCalcs(dataset = dat2)
    
    #### Calculate Population averaged mean kinship and GD
    PopAvgMKGD <- dat3 %>% pivot_longer(cols = c(coral1, coral2), values_to = "Coral_ID") %>% select(Coral_ID, avg_kinship) %>% arrange(Coral_ID, desc(avg_kinship)) %>% group_by(Coral_ID) %>% summarize(ind_mean_kinship = mean(avg_kinship)) %>% arrange(desc(ind_mean_kinship)) %>% ungroup() %>% summarise(PopAvgMK = mean(ind_mean_kinship), PopAvgGD = 1 - PopAvgMK)
    
    #### Calculate mean kinship for each colony
    MK_init <- dat3 %>% pivot_longer(cols = c(coral1, coral2), values_to = "Coral_ID") %>% select(Coral_ID, avg_kinship) %>% arrange(Coral_ID, desc(avg_kinship)) %>% group_by(Coral_ID) %>% summarize(ind_mean_kinship = mean(avg_kinship)) %>% arrange(desc(ind_mean_kinship))
    
    #### Set "highest_individual to "none" yet because it needs a value for the while loop below but we don't want to omit anyone (yet) and in the event that we want to keep ALL individuals (i.e., targetN = total number of individuals), this will ensure that we can do that; otherwise, we could only ever return total number - 1 individuals becuase one will always be the highest, and in the while look below that individual would be identified and omitted straight away. 
  highest_individual <- "none yet"
  
  #### Then run the above kinship calculations over and over, each time identifying and removing the individual with the highest mean kinship until we're left with targetN individuals. As long as the number of rows in coralAlleleDataRed is > targetN, the process will repeat until target N is reached  
  while (nrow(dat1) > targetN){
    dat1 <- dat1 %>% filter(!(Coral_ID %in% highest_individual))
    dat2 <- determineAllAlleleMatchesOthers(dataset = dat1)
    dat3 <- kinshipCalcs(dataset = dat2) 
    highest_individual <- dat3 %>% pivot_longer(cols = c(coral1, coral2), values_to = "Coral_ID") %>% select(Coral_ID, avg_kinship) %>% arrange(Coral_ID, desc(avg_kinship)) %>% group_by(Coral_ID) %>% summarize(ind_mean_kinship = mean(avg_kinship)) %>% arrange(desc(ind_mean_kinship)) %>% slice(1) %>% pull(Coral_ID)
  }
  #### Calculate mean kinship for each of the least related colonies
  MK_final <- dat3 %>% pivot_longer(cols = c(coral1, coral2), values_to = "Coral_ID") %>% select(Coral_ID, avg_kinship) %>% arrange(Coral_ID, desc(avg_kinship)) %>% group_by(Coral_ID) %>% summarize(ind_mean_kinship = mean(avg_kinship)) %>% arrange(desc(ind_mean_kinship))
  return(list(PopAvgMKGD = PopAvgMKGD, MK_init = MK_init, MK_final = MK_final))
  } 
}
