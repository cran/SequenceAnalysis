#' @title SequenceAnalysis.CodonUsage
#' @description Frequency of occurrence of synonymous codons
#' @details Codon Usage will be calculated for input sequence or by passing UniProt ID it will catch the sequence online from UniProt database and then calculate Codon Usage for it.
#' @author Babak Khorsand
#' @export SequenceAnalysis.CodonUsage
#' @param Nucleotide_Sequence Nucleotide Sequence
#' @param UniprotKB UniProt ID of desired sequence
#' @return CodonUsage
#' @examples
#' SequenceAnalysis.CodonUsage(Nucleotide_Sequence="atggctgcagcggccagtcacgatcagaggtaagttgtc")
#' SequenceAnalysis.CodonUsage(UniprotKB="O15131")
SequenceAnalysis.CodonUsage = function(Nucleotide_Sequence=NULL,UniprotKB=NULL)
{
  CodonUsage_List=NULL
  if (is.null(Nucleotide_Sequence))
  {
    if (is.null(UniprotKB))
    {
      stop("Nucleotide Sequence or UniprotKB must be set")
    }else
    {
      Nucleotide_Sequence=SequenceAnalysis.GetNucleotideSequence(UniprotKB)
      Nucleotide_Sequence=Nucleotide_Sequence[[4]]
    }
  }
  if (!is.null(Nucleotide_Sequence))
  {
    Nucleotide_Sequence=toupper(Nucleotide_Sequence)
    Codon=sapply(0:((nchar(Nucleotide_Sequence)/3)-1), function(x) substr(Nucleotide_Sequence,(x*3)+1,(x*3)+3))
    Codon_Table=table(Codon)
    codon_List=c("A,GCT,GCC,GCA,GCG","R,CGT,CGC,CGA,CGG,AGA,AGG","N,AAT,AAC",
                 "D,GAT,GAC","C,TGT,TGC","Q,CAA,CAG","E,GAA,GAG","G,GGT,GGC,GGA,GGG",
                 "H,CAT,CAC","I,ATT,ATC,ATA","L,TTA,TTG,CTT,CTC,CTA,CTG","K,AAA,AAG",
                 "M,ATG","F,TTT,TTC","P,CCT,CCC,CCA,CCG","S,TCT,TCC,TCA,TCG,AGT,AGC",
                 "T,ACT,ACC,ACA,ACG","W,TGG","Y,TAT,TAC","V,GTT,GTC,GTA,GTG","STOP,TAA,TGA,TAG")
    codon_List=strsplit(codon_List,",")
    CodonUsage_List=sapply(codon_List, function(x) sapply(x[-1], function(y) ifelse(is.na(Codon_Table[y]),0,Codon_Table[y])))
    for (i in 1:21)
    {
      names(CodonUsage_List[[i]])= codon_List[[i]][-1]
    }
    names(CodonUsage_List)=sapply(codon_List, function(x) x[1])
    CodonUsage_List=sapply(CodonUsage_List, function(x) {sum=sum(x); sapply(x, function(y) ifelse(sum==0,0,y/sum))})
    CodonUsage_List=unlist(CodonUsage_List)
    names(CodonUsage_List)=gsub(".*\\.(...)","\\1",names(CodonUsage_List))
  }
  return(CodonUsage_List)
}
