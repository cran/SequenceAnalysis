#' @title SequenceAnalysis.StackingEnergy
#' @description Stacking Energy obtained by NN model
#' @details The NN model for nucleic acids assumes that the stability of a given base pair depends on the identity and orientation of neighboring base pairs. Stacking Energy = DeltaG(total) = Sigma (n(i)*DeltaG(i)) + DeltaG(init) + DeltaG(end) + DeltaG(sym), which DeltaG for i, init and end is obtained by Unified NN free energy parameter. Symmetry of self-complementary duplexes is also included by DeltaG(sym) equals to +0.43 kcal/mol if the duplex is self-complementary and zero if it is non-self-complementary.
#' @author Babak Khorsand
#' @export SequenceAnalysis.StackingEnergy
#' @param Nucleotide_Sequence Nucleotide Sequence
#' @param UniprotKB UniProt ID of desired sequence
#' @param CDS if TRUE Stacking Energy of CDS Region will be calculated
#' @return StackingEnergy
#' @examples
#' SequenceAnalysis.StackingEnergy("actagtcacgatcag")
#' SequenceAnalysis.StackingEnergy(UniprotKB="O15131")
#' SequenceAnalysis.StackingEnergy(UniprotKB="O15131",CDS=FALSE)
SequenceAnalysis.StackingEnergy = function(Nucleotide_Sequence=NULL,UniprotKB=NULL,CDS=TRUE)
{
  StackingEnergy="N/A"
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
    Nuc_Length=nchar(Nucleotide_Sequence)
    StackingEnergy=sapply(1:(Nuc_Length-1), function(i) StackingEnergy_Table[substr(Nucleotide_Sequence,i,i+1)])
    ReverseComplement = SequenceAnalysis.ReverseComplement(Nucleotide_Sequence)
    StackingEnergy = ifelse(substr(Nucleotide_Sequence,1,1) %in% c("A","T"),StackingEnergy_Table["init_AT"],StackingEnergy_Table["init_GC"])+
      ifelse(substr(Nucleotide_Sequence,Nuc_Length,Nuc_Length) %in% c("A","T"),StackingEnergy_Table["init_AT"],StackingEnergy_Table["init_GC"])+
      ifelse(Nucleotide_Sequence==ReverseComplement,StackingEnergy_Table["sym"],0)-
      sum(StackingEnergy)
    if (is.null(StackingEnergy))
      StackingEnergy="N/A"
    names(StackingEnergy)="StackingEnergy"
  }
  return(StackingEnergy)
}
