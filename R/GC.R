#' @title SequenceAnalysis.GC
#' @description GC content : percentage of nucleotide g and c
#' @details Calculate GC content for input sequence or by passing UniProt ID it will catch the sequence online from UniProt database and then calculate GC content for it.
#' @author Babak Khorsand
#' @export SequenceAnalysis.GC
#' @param Nucleotide_Sequence Nucleotide Sequence
#' @param UniprotKB UniProt ID of desired sequence
#' @param CDS if TRUE GC content of CDS Region will be calculated
#' @return GC content
#' @examples
#' SequenceAnalysis.GC("actagtcacgatcag")
#' SequenceAnalysis.GC(UniprotKB="O15131")
#' SequenceAnalysis.GC(UniprotKB="O15131",CDS=TRUE)
SequenceAnalysis.GC = function(Nucleotide_Sequence=NULL,UniprotKB=NULL,CDS=FALSE)
{
  if (is.null(Nucleotide_Sequence))
  {
    if (is.null(UniprotKB))
    {
      stop("Nucleotide Sequence or UniprotKB must be set")
    }else
    {
      Nucleotide_Sequence=SequenceAnalysis.GetNucleotideSequence(UniprotKB)
      if (CDS)
      {
        Nucleotide_Sequence=Nucleotide_Sequence[[3]]
        if (Nucleotide_Sequence=="")
          Nucleotide_Sequence=Nucleotide_Sequence[[4]]
      }else
      {
        Nucleotide_Sequence=Nucleotide_Sequence[[4]]
      }
    }
  }
  if (!is.null(Nucleotide_Sequence))
  {
    Nucleotide_Sequence=tolower(Nucleotide_Sequence)
    GC=round((table(strsplit(Nucleotide_Sequence,""))["g"]+table(strsplit(Nucleotide_Sequence,""))["c"])/nchar(Nucleotide_Sequence),3)
    names(GC)="GC"
  }
  return(GC)
}
