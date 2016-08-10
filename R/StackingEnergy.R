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
#' SequenceAnalysis.StackingEnergy(UniprotKB="O15131",CDS=TRUE)
SequenceAnalysis.StackingEnergy = function(Nucleotide_Sequence=NULL,UniprotKB=NULL,CDS=FALSE)
{
  StackingEnergy=NULL
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
    Nuc_Length=nchar(Nucleotide_Sequence)
    StackingEnergy_Table=c(AA=1,AT=0.88,AC=1.44,AG=1.28,TA=0.58,TT=1,TC=1.3,TG=1.45,CA=1.45,CT=1.28,CC=1.84,CG=2.17,GA=1.3,GT=1.44,GC=2.24,GG=1.84,init_GC=0.98,init_AT=1.03,sym=0.43)
    StackingEnergy=sapply(1:(Nuc_Length-1), function(i) StackingEnergy_Table[substr(Nucleotide_Sequence,i,i+1)])
    ReverseComplement = SequenceAnalysis.ReverseComplement(Nucleotide_Sequence)
    StackingEnergy = ifelse(substr(Nucleotide_Sequence,1,1) %in% c("A","T"),StackingEnergy_Table["init_AT"],StackingEnergy_Table["init_GC"])+
      ifelse(substr(Nucleotide_Sequence,Nuc_Length,Nuc_Length) %in% c("A","T"),StackingEnergy_Table["init_AT"],StackingEnergy_Table["init_GC"])+
      ifelse(Nucleotide_Sequence==ReverseComplement,StackingEnergy_Table["sym"],0)-
      sum(StackingEnergy)
  }
  return(StackingEnergy)
}
