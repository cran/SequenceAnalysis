#' @title SequenceAnalysis.CodonToAminoAcid
#' @description Converting desired codons to related amino acids.
#' @details Converting nucleotide sequence to protein sequence.
#' @author Babak Khorsand
#' @export SequenceAnalysis.CodonToAminoAcid
#' @param Nucleotide_Sequence Nucleotide Sequence
#' @return Related amino acid
#' @examples
#' SequenceAnalysis.CodonToAminoAcid(Nucleotide_Sequence="atggctgcagcggccagtcacgatcagaggtaagttgtc")
SequenceAnalysis.CodonToAminoAcid = function(Nucleotide_Sequence=NULL)
{
  if (is.null(Nucleotide_Sequence))
  {
    stop("Nucleotide Sequence must be set")
  }
  Nucleotide_Sequence=toupper(Nucleotide_Sequence)
  Codon=sapply(0:((nchar(Nucleotide_Sequence)/3)-1), function(x) substr(Nucleotide_Sequence,(x*3)+1,(x*3)+3))
  codon_List=strsplit(codon_List,",")
  Codon=sapply(Codon, function(x) codon_List[grep(x,codon_List)][[1]][1])
  if (Codon[length(Codon)]=="Z")
    Codon=Codon[-length(Codon)]
  Codon=paste(Codon,collapse = "")
  if (is.null(Codon))
    Codon="N/A"
  Codon=toupper(Codon)
  names(Codon)="ProteinSequence"
  return(Codon)
}
