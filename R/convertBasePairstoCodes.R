#' Convert base pairs to single character IUPAC codes
#' @description This function replaces the base pairs with `IUPAC` codes and converts the data to wide format so that there is a column per locus. If any invalid characters are present in a column, then the entire column is deleted, a warning is issued, and the column names are reported.
#' @param initdata A data frame containing a column `Coral_ID` with unique value per row, as well as at least one column containing base pair data with at least one value matching valid allele data in the `IUPAC` table.
#' @returns This function returns a data frame identical to the input data frame except the paired alleles (e.g., G:) are replace by single allele codes (e.g., G). Note that only "allowable data" are included in the resulting data frame, and any columns not adhering to the data formatting requirements are omitted by `checkforAllowableData()` and listed in a warning issued by `handleError_ProhibitedData()`. 
#' @importFrom dplyr select left_join
#' @importFrom tidyr pivot_longer pivot_wider
#' @export
convertBasePairstoCodes <- function(initdata) {
  # create file raw and select only columns that (1) contain valid allele pairing data as listed in the IUPAC data file and (2) are named Coral_ID.
    raw <- initdata %>%
    select(
      names(initdata)[colSums(apply(initdata, 2, checkforAllowableData)) == nrow(initdata)|names(initdata)=="Coral_ID"]
    )
  #pivots the raw data to only 3 columns
  longraw <- raw %>%
    pivot_longer(
      !Coral_ID,
      names_to = "Locus",
      values_to = "Allelepairs"
    )
  #joins the raw table to the cipher table
  translated <- longraw %>%
    left_join(
      IUPAC,
      by = "Allelepairs"
    )
  #removes the column with the 3 character alleles (colin deleted .data$Allelepairs because apparently it's deprecated 01/16/25)
  clean <- translated %>%
    select(-.data$Allelepairs)
  #pivots the table back out to wide form with each locus as a column
  wideclean <- clean %>%
    pivot_wider(
      names_from = "Locus",
      values_from = "IUPACallele"
    )
  wideclean[wideclean=="?"] <- NA # Changes all the "?" to nulls for the sake of the genet calculations.
  return(wideclean)
}
