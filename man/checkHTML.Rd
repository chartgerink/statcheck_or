\name{checkHTML}
\alias{checkHTML}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extract test statistics from HTML file.
}
\description{
Extracts statistical references from given HTML files.
}
\usage{
checkHTML(files, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{files}{
Vector of strings containing file paths to HTML files to check.
}
  \item{\dots}{
Arguments sent to  \code{\link{statcheck}}.
}
}
\details{
See \code{\link{statcheck}} for more details. Use \code{\link{checkHTMLdir}} to import al HTML files in a given directory at once. 

Note that the conversion to plain text and extraction of statistics can result in errors. Some statistical values can be missed, especially if the notation is unconvetional. It is recommended to manually check some of the results.
}
\value{
A data frame containing for each extracted statistic:
\item{Source}{Name of the file of which the statistic is extracted}
\item{Statistic}{Character indicating the statistic that is extracted}
\item{df1}{First degree of freedom}
\item{df2}{Second degree of freedom (if applicable)}
\item{Value}{Reported value of the statistic}
\item{Reported.Comparison}{Reported comparison, when importing from pdf this will often not be converted properly}
\item{Reported.P.Value}{The reported p-value, or NA if the reported value was NS}
\item{Computed}{The re-computed p-value}
\item{OneTail}{Smallest one-tailed recoputed p-value (if applicable)}
\item{Location}{Location of the statistic in the string.}
\item{Raw}{Raw string of the statistical reference that is extracted}
}
\author{
Sacha Epskamp <mail@sachaepskamp.com> & Michele Nuijten
}
\seealso{
\code{\link{statcheck}}, \code{\link{checkPDF}}, \code{\link{checkPDFdir}}, \code{\link{checkHTMLdir}}, \code{\link{checkdir}}
}