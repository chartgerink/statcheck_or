#' Extract statistics and recompute p-values.
#' 
#' This function extracts statistics from strings and returns the extracted values, reported p-values and recomputed p-values. The package relies on the program "pdftotext", see the paragraph "Note" for details on the installation.
#' 
#' @param x A vector of strings.
#' @param stat "t" to extract t-values, "F" to extract F-values, "cor" to extract correlations, "chisq"to extract chi-square values, and "Z" to extract Z-values.
#' @param OneTailedTests Logical. Do we assume that all reported tests are one tailed (TRUE) or two tailed (FALSE, default)?
#' @param alpha Assumed level of significance in the scanned texts. Defaults to .05. 
#' @param pEqualAlphaSig Logical. If TRUE, statcheck counts p <= alpha as significant (default), if FALSE, statcheck counts p < alpha as significant
#' @param OneTailedTxt Logical. If TRUE, statcheck searches the text for "one-sided", "one-tailed", and "directional" to identify the possible use of one-sided tests. If one or more of these strings is found in the text AND the result would have been correct if it was a one-sided test, the result is assumed to be indeed one-sided and is counted as correct.
#' @param AllPValues Logical. If TRUE, the output will consist of a dataframe with all detected p values, also the ones that were not part of the full results in APA format
#'
#' @details Statcheck uses regular expressions to find statistical results in APA format. When a statistical result deviates from APA format, statcheck will not find it. The APA formats that statcheck uses are: t(df) = value, p = value; F(df1,df2) = value, p = value; r(df) = value, p = value; [chi]2 (df, N = value) = value, p = value (N is optional, delta G is also included); Z = value, p = value. All regular expressions take into account that test statistics and p values may be exactly (=) or inexactly (< or >) reported. Different spacing has also been taken into account. 
#'  This function can be used if the text of articles has already been imported in R. To import text from pdf files and automatically send the results to this function use \code{\link{checkPDFdir}} or \code{\link{checkPDF}}. To import text from HTML files use the similar functions \code{\link{checkHTMLdir}} or \code{\link{checkHTML}}. Finally, \code{\link{checkdir}} can be used to import text from both PDF and HTML files in a folder.
#'  Note that the conversion from PDF (and sometimes also HTML) to plain text and extraction of statistics can result in errors. Some statistical values can be missed, especially if the notation is unconventional. It is recommended to manually check some of the results.
#'  PDF files should automatically be converted to plain text files. However, if this does not work, it might help to manually install the program "pdftotext". You can obtain pdftotext from \code{http://www.foolabs.com/xpdf/download.html}. Download and unzip the precompiled binaries. Next, add the folder with the binaries to the PATH variables so that this program can be used from command line.
#'  Also, note that a seemingly inconsistent p value can still be correct when we take into account that the test statistic might have been rounded after calculating the corresponding p value. For instance, a reported t value of 2.35 could correspond to an actual value of 2.345 to 2.354 with a range of p values that can slightly deviate from the recomputed p value. Statcheck will not count cases like this as errors.
#' @seealso \code{\link{checkPDF}}, \code{\link{checkHTMLdir}}, \code{\link{checkHTML}}, \code{\link{checkdir}}
#'
#' @return A data frame containing for each extracted statistic:
#' \item{Source}{Name of the file of which the statistic is extracted}
#' \item{Statistic}{Character indicating the statistic that is extracted}
#' \item{df1}{First degree of freedom}
#' \item{df2}{Second degree of freedom (if applicable)}
#' \item{Test.Comparison}{Reported comparison of the test statistic, when importing from pdf this will often not be converted properly}
#' \item{Value}{Reported value of the statistic}
#' \item{Reported.Comparison}{Reported comparison, when importing from pdf this might not be converted properly}
#' \item{Reported.P.Value}{The reported p-value, or NA if the reported value was NS}
#' \item{Computed}{The recomputed p-value}
#' \item{Raw}{Raw string of the statistical reference that is extracted}
#' \item{Error}{The computed p value is not congruent with the reported p value}
#' \item{DecisionError}{The reported result is significant whereas the recomputed result is not, or vice versa.}
#' \item{OneTail}{Logical. Is it likely that the reported p value resulted from a correction for one-sided testing?}
#' \item{OneTailedInTxt}{Logical. Does the text contain the string "sided", "tailed", and/or "directional"?}
#' \item{CopyPaste}{Logical. Does the exact string of the extracted raw results occur anywhere else in the article?}
#' @export
#'
#' @examples
statcheck <- function(
  x,
  contextlength = 200
){
  # Create empty objects
  raw_result = NULL
  genotype = NULL
  snp = NULL
  dv = NULL
  or_comparison = NULL
  or_result = NULL
  ci_confidence = NULL
  ci_lb = NULL
  ci_ub = NULL
  p_comparison = NULL
  p_result  =  NULL
  p_recalc  =  NULL
  test_stat_recalc <- NULL
  se_recalc = NULL
  
  if (length(x)==0) return(NULL)
  
  if (is.null(names(x))) names(x) <-  1:length(x)
  
  # create counter for loops
  index <- 1
  
  message("Extracting statistics...")
  pb <- txtProgressBar(max=length(x),style=3)
  for (i in 1:length(x)){
    
    txt <- x[i]
    
    # identify sequence of results and extract their text
    resLoc <- gregexpr("(genotype).*?rs[0-9]{1,10}.*?(associate|relate|correlate)(d)?\\s(with|to).*?(odds\\sratio|\\(?OR\\)?).*?[0-9]{2}\\%\\s(confidence\\sinterval|\\(?ci\\)?).*?(p.*?\\s?[0-9]?.[0-9]{1,5})",
                       txt,
                       ignore.case = TRUE)[[1]]
    resContext <- substring(txt,
                            resLoc,
                            resLoc + attr(resLoc, "match.length") + 3)
    
    # locate data in context strings
    # genotype
    locator_genotype <- gregexpr("genotype.?[A-Za-z]{2}", resContext, ignore.case = TRUE)
    
    # SNP
    locator_snp <- gregexpr("\\srs[0-9]{1,10}\\s", resContext, ignore.case = TRUE)
    
    # DV
    locator_ass <- gregexpr("(associat[a-z]{1,}|
                           correlat[a-z]{1,}|
                           relat[a-z]{1,})\\s(with|to)\\s", resContext, ignore.case = TRUE)
    
    locator_dv <- gregexpr("(associat[a-z]{1,}|
                           correlat[a-z]{1,}|
                           relat[a-z]{1,})\\s(with|to)\\s[a-z]{1,}", resContext, ignore.case = TRUE)
    
    # Odds ratio
    locator_or <- gregexpr("(odds ratio.*|or.*?)[<>=]", resContext, ignore.case = TRUE)
    
    # Confidence interval
    locator_ci <- gregexpr("[0-9]{2}\\%.*?(confidence interval.*|ci.*?)[=:;].*?([0-9]{1,}?[.]?[0-9]{1,}.*?[0-9]{1,}?[.]?[0-9]{1,})", resContext, ignore.case = TRUE)
    
    # P-value (if present)
    locator_p <- gregexpr("p.?[<>=]", resContext, ignore.case = TRUE)
    
    for (j in 1:length(resContext)){
      # extract data from each context with use of locators
      # genotype
      genotype_ind <- substr(resContext[[j]], locator_genotype[[j]] + 10 - 1, locator_genotype[[j]] + 10)
      # SNP
      snp_temp <- str_match_all(resContext[[j]],
                                "\\srs[0-9]{1,10}\\s")[[1]][1]
      snp_ind <- gsub(pattern = "\\s", x = snp_temp, replacement = "")
      # DV
      dv_ind <- gsub(pattern = "\\s",
                     x = substr(resContext[[j]], locator_ass[[j]] + attr(locator_ass[[j]], "match.length"), locator_dv[[j]] + attr(locator_dv[[j]], "match.length")),
                     "")
      # OR
      or_temp <- substr(resContext[[j]], locator_or[[j]], locator_or[[j]] + attr(locator_or[[j]], "match.length") + 5)
      or_comp <- str_match_all(or_temp,
                               "[<>=]")[[1]][1]
      or_ind <- as.numeric(str_match_all(or_temp,
                              "[0-9]{1,}?[.][0-9]{1,3}")[[1]][1])
      # CI
      ci_temp <- substr(resContext[[j]], locator_ci[[j]], locator_ci[[j]] + attr(locator_ci[[j]], "match.length") + 5)
      ci_conf <- as.numeric(substr(str_match_all(ci_temp, "[0-9]{1,2}\\%\\s?[cC]")[[1]][,1], 0, 2))
      ci_conc <- str_match_all(ci_temp, "[0-9]{1,}?[.][0-9]{1,3}")[[1]][,1]
      ci_lb_ind <- as.numeric(min(ci_conc))
      ci_ub_ind <- as.numeric(max(ci_conc))
      # P-value
      p_temp <- substr(resContext[[j]], locator_p[[j]], locator_p[[j]] + attr(locator_p[[j]], "match.length") + 5)
      p_comp <- str_match_all(p_temp,
                              "[<>=]")[[1]][1]
      p_ind <- as.numeric(str_match_all(p_temp,
                             "[0-9]{1,}?[.][0-9]{1,3}")[[1]][1])
      
      # recalculate p-value
      # from Altman, D. G., & Bland, J. M. (2011). How to obtain the P value from a confidence interval. BMJ , 343, d2304.
      se_ind <- (ci_ub_ind - ci_lb_ind) / (2 * qnorm((1 - (ci_conf / 100)) / 2, lower.tail = FALSE))
      z_ind <- or_ind / se_ind
      p_recalc_ind <- pnorm(z_ind, lower.tail = FALSE) * 2
      
      
      # write back results into main objects
      raw_result[index] = resContext[[j]]
      genotype[index] = genotype_ind
      snp[index] = snp_ind
      dv[index] = dv_ind
      or_comparison[index] = or_comp
      or_result[index] = or_ind
      ci_confidence[index] = ci_conf
      ci_lb[index] = ci_lb_ind
      ci_ub[index] = ci_ub_ind
      p_comparison[index] = p_comp
      p_result [index] =  p_ind
      p_recalc [index] =  p_recalc_ind
      test_stat_recalc[index] = z_ind
      se_recalc[index] = se_ind
      
      index = index + 1
    }
  }
  
  
  # data frame
  Res <- data.frame(Source = names(x),
                    raw_result = raw_result,
                    genotype = genotype,
                    snp = snp,
                    dv = dv,
                    or_comparison = or_comparison,
                    or_result = or_result,
                    ci_confidence = ci_confidence,
                    ci_lb = ci_lb,
                    ci_ub = ci_ub,
                    p_comparison = p_comparison,
                    p_result  = p_result,
                    p_recalc  = p_recalc,
                    test_stat_recalc = test_stat_recalc,
                    se_recalc = se_recalc)

class(Res) <- c("statcheck","data.frame")

###--------------------------------------------------------------------- 
# Return message when there are no results
if(nrow(Res)>0){
  
  
  return(Res) 
} else {
  Res <- cat("statcheck did not find any results\n")
}
}


###########################

r2t <- function(# Transform r values into t values
  ### Function to transform r values into t values by use of raw r and degrees of freedom.
  r,
  ### Raw correlation value
  df
  ### Degrees of freedom (N-1)
){
  r / (sqrt((1-r^2)/df))
}