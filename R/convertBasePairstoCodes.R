#' Convert to long format and then convert base pairs to IUPAC codes
#' @description This function replaces the base pairs with IUPAC codes and
#' converts the data to wide format so that there is a column per code.
#' # Make sure that the data we're using are ONLY either the
#' Coral_ID OR one of the character strings included in OnlyDataAllowed. If
#' there are no such columns then nothing happens. There's probably a more
#' elegant way to do this but it works fine.
#' @param initdata A dataset containing a column "Coral_ID" with unique value per row, as well as at least one column containing base pair data with at least one value matching base pair data in IUPAC.
#' @importFrom dplyr select left_join
#' @importFrom tidyr pivot_longer pivot_wider
#' @export
convertBasePairstoCodes <- function(initdata) {
  raw <- initdata %>%
    select(
      names(initdata)[colSums(apply(initdata, 2, filterOnlyDataAllowed)) > 0 |
                        names(initdata)=="Coral_ID"]
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
      by="Allelepairs"
    )
  #removes the column with the 3 character alleles (colin deleted .data$Allelepairs because apparently it's deprecated 11625)
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
