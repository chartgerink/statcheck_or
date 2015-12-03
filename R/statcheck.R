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
  # Create empty data frame for main result:
  Res <- data.frame(Source = NULL,
                    context = NULL,
                    genotype = NULL,
                    SNP = NULL,
                    gene = NULL,
                    dependent = NULL,
                    or = NULL,
                    ci_lb = NULL,
                    ci_ub = NULL,
                    pval  =  NULL,
                    se_recalc = NULL,
                    pval_recalc  =  NULL)
  class(Res) <- c("statcheck","data.frame")
  
  if (length(x)==0) return(Res)
  
  if (is.null(names(x))) names(x) <-  1:length(x)
  
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
    locator_dv <- gregexpr("(associat[a-z]{1,2}|
                           correlat[a-z]{1,2}|
                           relat[a-z]{1,2})\\s(with|to)", resContext, ignore.case = TRUE)
    
    # Odds ratio
    locator_or <- gregexpr("(odds ratio.*|or.*?)[<>=]", resContext, ignore.case = TRUE)
    
    # Confidence interval
    locator_ci <- gregexpr("[0-9]{2}\\%.*?(confidence interval.*|ci.*?)[=:;]", resContext, ignore.case = TRUE)
    
    # P-value (if present)
    locator_p <- gregexpr("p.?[<>=]", resContext, ignore.case = TRUE)
    
    for (j in 1:length(resContext)){
      # extract data from each context with use of locators
      # genotype
      genotype_ind <- str_sub(resContext[[j]], locator_genotype[[j]] + 10 - 1, locator_genotype[[j]] + 10)
      print(genotype_ind)
      # SNP
      snp_ind <- str_sub(resContext[[j]], locator_snp[[j]] + 1, locator_snp[[j]] + attr(locator_snp[[j]], "match.length") - 2)
      # DV
      dv_ind <- str_sub(resContext[[j]], locator_dv[[j]] + 5 - 1, locator_snp[[j]] + attr(locator_snp[[j]], "match.length") - 1)
      # OR
      # CI
      # P-value
      
    }
    
    genotype
    
    locator_genotype <- gregexpr("genotype\\s[A-Za-z]{2}", txt, ignore.case = TRUE)
    
    # SNP
    locator_snp <- gregexpr("SNP\\srs[0-9]{1,10}", txt, ignore.case = TRUE)
    
    # DV
    locator_dv <- gregexpr("associat.*", txt, ignore.case = TRUE)
    
    # Odds ratio
    locator_or <- gregexpr("(odds ratio.*|or.*?)[<>=]", txt, ignore.case = TRUE)
    
    # Confidence interval
    locator_ci <- gregexpr("[0-9]{2}\\%.*?(confidence interval.*|ci.*?)[=:;]", txt, ignore.case = TRUE)
    
    # P-value (if present)
    locator_p <- gregexpr("p.?[<>=]", txt, ignore.case = TRUE)
    
    
      # extract data from each context with use of locators
      # genotype
      genotype_ind <- str_sub(txt, locator_genotype[[1]] + 10 - 1, locator_genotype[[1]] + 10)
      print(genotype_ind)
      
      snp_ind <- str_sub(txt, locator_snp[[1]], locator_snp[[1]] + 20)
      print(snp_ind)
      # SNP
      # DV
      # OR
      # CI
      # P-value
      
    }
    
    
    # final data frame
    Res <- data.frame(Source = names(x),
                      context = resContext,
                      genotype = NULL,
                      SNP = NULL,
                      gene = NULL,
                      dependent = NULL,
                      or = NULL,
                      ci_lb = NULL,
                      ci_ub = NULL,
                      pval  =  NULL,
                      se_recalc = NULL,
                      pval_recalc  =  NULL)
    )
    
    class(Res) <- c("statcheck","data.frame")
  }
  
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