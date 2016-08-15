#' @title SequenceAnalysis.AminoAcidToCodon
#' @description All related codons of desired amino acid.
#' @details All related codons of desired amino acid.
#' @author Babak Khorsand
#' @export SequenceAnalysis.AminoAcidToCodon
#' @param AminoAcid Amino Acid
#' @return Related Codons
#' @examples
#' SequenceAnalysis.AminoAcidToCodon(AminoAcid="A")
SequenceAnalysis.AminoAcidToCodon = function(AminoAcid=NULL)
{
  if (is.null(AminoAcid))
  {
    stop("Enter Amino Acid in one letter code")
  }
  if (nchar(AminoAcid)>1)
  {
    stop("Enter Amino Acid in one letter code")
  }
  AminoAcid=tolower(AminoAcid)
  Codon=codon_List[grep(AminoAcid,codon_List)]
  Codon=unlist(strsplit(Codon,","))
  Codon=Codon[-1]
  if (is.null(Codon))
    Codon="N/A"
  return(Codon)
}
