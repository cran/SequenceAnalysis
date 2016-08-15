#' @title SequenceAnalysis.Complement
#' @description Complement of desired nucleotide sequence.
#' @details Complement of desired nucleotide sequence. (A for T, C for G and vice versa) (tgacta -> Complement: actgat)
#' @author Babak Khorsand
#' @export SequenceAnalysis.Complement
#' @param Nucleotide_Sequence Nucleotide Sequence
#' @param UniprotKB UniProt ID of desired sequence
#' @param CDS if TRUE Stacking Energy of CDS Region will be calculated
#' @return Complement
#' @examples
#' SequenceAnalysis.Complement("actagtcacgatcag")
#' SequenceAnalysis.Complement(UniprotKB="O15131")
#' SequenceAnalysis.Complement(UniprotKB="O15131",CDS=TRUE)
SequenceAnalysis.Complement = function(Nucleotide_Sequence=NULL,UniprotKB=NULL,CDS=FALSE)
{
  Complement=NULL
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
        Nucleotide_Sequence=Nucleotide_Sequence[[2]]
        if (Nucleotide_Sequence=="N/A")
          Nucleotide_Sequence=Nucleotide_Sequence[[3]]
      }else
      {
        Nucleotide_Sequence=Nucleotide_Sequence[[3]]
      }
    }
  }
  if (Nucleotide_Sequence!="N/A")
  {
    Nucleotide_Sequence=toupper(Nucleotide_Sequence)
    Nucleotide_Sequence=unlist(strsplit(Nucleotide_Sequence,""))
    a=grep("A",Nucleotide_Sequence)
    t=grep("T",Nucleotide_Sequence)
    c=grep("C",Nucleotide_Sequence)
    g=grep("G",Nucleotide_Sequence)
    Nucleotide_Sequence[a]="T"
    Nucleotide_Sequence[t]="A"
    Nucleotide_Sequence[c]="G"
    Nucleotide_Sequence[g]="C"
    Complement=paste(Nucleotide_Sequence,collapse = "")
  }
  if (is.null(Complement))
    Complement="N/A"
  names(Complement)="Complement"
  return(Complement)
}
