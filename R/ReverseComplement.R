#' @title SequenceAnalysis.ReverseComplement
#' @description Reverse-Complement of desired nucleotide sequence.
#' @details Reverse-Complement of desired nucleotide sequence. (atcagt -> Reverse: tgacta -> Reverse-Complement: actgat)
#' @author Babak Khorsand
#' @export SequenceAnalysis.ReverseComplement
#' @param Nucleotide_Sequence Nucleotide Sequence
#' @param UniprotKB UniProt ID of desired sequence
#' @param CDS if TRUE Reverse-Complement of CDS Region will be calculated
#' @return Reverse-Complement
#' @examples
#' SequenceAnalysis.ReverseComplement("actagtcacgatcag")
#' SequenceAnalysis.ReverseComplement(UniprotKB="O15131")
#' SequenceAnalysis.ReverseComplement(UniprotKB="O15131",CDS=TRUE)
SequenceAnalysis.ReverseComplement = function(Nucleotide_Sequence=NULL,UniprotKB=NULL,CDS=FALSE)
{
  ReverseComplement=NULL
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
    Nucleotide_Sequence=Nucleotide_Sequence[length(Nucleotide_Sequence):1]
    a=grep("A",Nucleotide_Sequence)
    t=grep("T",Nucleotide_Sequence)
    c=grep("C",Nucleotide_Sequence)
    g=grep("G",Nucleotide_Sequence)
    Nucleotide_Sequence[a]="T"
    Nucleotide_Sequence[t]="A"
    Nucleotide_Sequence[c]="G"
    Nucleotide_Sequence[g]="C"
    ReverseComplement=paste(Nucleotide_Sequence,collapse = "")
  }
  return(ReverseComplement)
}
