#' @title SequenceAnalysis.Reverse
#' @description Reverse of desired nucleotide sequence.
#' @details Reverse of desired nucleotide sequence. (atcagt -> Reverse: tgacta)
#' @author Babak Khorsand
#' @export SequenceAnalysis.Reverse
#' @param Nucleotide_Sequence Nucleotide Sequence
#' @param UniprotKB UniProt ID of desired sequence
#' @param CDS if TRUE Reverse of CDS Region will be calculated
#' @return Reverse
#' @examples
#' SequenceAnalysis.Reverse("actagtcacgatcag")
#' SequenceAnalysis.Reverse(UniprotKB="O15131")
#' SequenceAnalysis.Reverse(UniprotKB="O15131",CDS=TRUE)
SequenceAnalysis.Reverse = function(Nucleotide_Sequence=NULL,UniprotKB=NULL,CDS=FALSE)
{
  Reverse=NULL
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
    Nucleotide_Sequence=toupper(Nucleotide_Sequence)
    Nucleotide_Sequence=unlist(strsplit(Nucleotide_Sequence,""))
    Reverse=Nucleotide_Sequence[length(Nucleotide_Sequence):1]
    Reverse=paste(Reverse,collapse = "")
  }
  return(Reverse)
}
