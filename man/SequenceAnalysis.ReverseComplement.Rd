% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ReverseComplement.R
\name{SequenceAnalysis.ReverseComplement}
\alias{SequenceAnalysis.ReverseComplement}
\title{SequenceAnalysis.ReverseComplement}
\usage{
SequenceAnalysis.ReverseComplement(Nucleotide_Sequence = NULL,
  UniprotKB = NULL, CDS = FALSE)
}
\arguments{
\item{Nucleotide_Sequence}{Nucleotide Sequence}

\item{UniprotKB}{UniProt ID of desired sequence}

\item{CDS}{if TRUE Reverse-Complement of CDS Region will be calculated}
}
\value{
Reverse-Complement
}
\description{
Reverse-Complement of desired nucleotide sequence.
}
\details{
Reverse-Complement of desired nucleotide sequence. (atcagt -> Reverse: tgacta -> Reverse-Complement: actgat)
}
\examples{
SequenceAnalysis.ReverseComplement("actagtcacgatcag")
SequenceAnalysis.ReverseComplement(UniprotKB="O15131")
SequenceAnalysis.ReverseComplement(UniprotKB="O15131",CDS=TRUE)
}
\author{
Babak Khorsand
}

