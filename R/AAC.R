#' @title SequenceAnalysis.AAC
#' @description Amino Acid Composition
#' @details Amino acid Composition is calculated by four different methods: a) Twenty-two independent categories are considered, with one amino acid for each category. b) Five categories (Nonpolar Aliphatic, Nonpolar Aromatic, Polar Uncharged, Polar Positively Charged, Polar Negatively Charged) are considered according to the standard chemical structures of amino acids. c) Six categories (Nonpolar Aliphatic, Nonpolar Aromatic, Polar Uncharged, Polar Positively Charged, Polar Negatively Charged, Special cases) are considered which Cysteine, Selenocysteine, Glycine and Proline are placed in Special cases group. d) Eight categories are clustered via k-means algorithm on Physicochemical index of amino acids.
#' @author Babak Khorsand
#' @export SequenceAnalysis.AAC
#' @param Protein_Sequence Protein Sequence or Uniprot ID of a sequence
#' @param Sequence If False it will get sequence from uniprot
#' @param Groups 1 for usual 22 groups amino acid composition, 2 for 5 groups (Nonpolar Aliphatic, nonpolar Aromatic, Polar Uncharged, Polar Positively Charged, Polar Negatively Charged), 3 for 6 groups (nonpolar Aliphatic, nonpolar Aromatic, Polar Uncharged, Polar Positively Charged, Polar Negatively Charged, Special cases) and 4 for 8 group (k-means clustering result) amino acid composition.
#' @return Aminio acid composition
#' @examples
#' SequenceAnalysis.AAC("AKMNAAKMQWYVIGLPCERTDRSCTRQWYVPIG")
#' SequenceAnalysis.AAC("O15131",Sequence=FALSE)
#' SequenceAnalysis.AAC("O15131",Sequence=FALSE,Groups=2)
SequenceAnalysis.AAC=function(Protein_Sequence,Sequence=TRUE,Groups=1)
{
  if (!Sequence)
  {
    Protein_Sequence=SequenceAnalysis.GetProteinSequence(Protein_Sequence)
    if (Protein_Sequence=="N/A")
      stop("Protein Sequence is not available")
  }
  if (length(Protein_Sequence)>0)
  {
    Protein_Sequence=toupper(Protein_Sequence)
    SeqLength=nchar(Protein_Sequence)
    Protein_Seq = strsplit(Protein_Sequence,"")
    Protein_Seq=table(Protein_Seq)
    AAC=c(A=0,R=0,N=0,D=0,C=0,E=0,Q=0,G=0,H=0,I=0,L=0,K=0,M=0,F=0,P=0,S=0,T=0,W=0,Y=0,V=0,U=0,X=0)
    AAC[1]=round(ifelse(is.na(Protein_Seq["A"]),0,Protein_Seq["A"])/SeqLength,5)
    AAC[2]=round(ifelse(is.na(Protein_Seq["R"]),0,Protein_Seq["R"])/SeqLength,5)
    AAC[3]=round(ifelse(is.na(Protein_Seq["N"]),0,Protein_Seq["N"])/SeqLength,5)
    AAC[4]=round(ifelse(is.na(Protein_Seq["D"]),0,Protein_Seq["D"])/SeqLength,5)
    AAC[5]=round(ifelse(is.na(Protein_Seq["C"]),0,Protein_Seq["C"])/SeqLength,5)
    AAC[6]=round(ifelse(is.na(Protein_Seq["E"]),0,Protein_Seq["E"])/SeqLength,5)
    AAC[7]=round(ifelse(is.na(Protein_Seq["Q"]),0,Protein_Seq["Q"])/SeqLength,5)
    AAC[8]=round(ifelse(is.na(Protein_Seq["G"]),0,Protein_Seq["G"])/SeqLength,5)
    AAC[9]=round(ifelse(is.na(Protein_Seq["H"]),0,Protein_Seq["H"])/SeqLength,5)
    AAC[10]=round(ifelse(is.na(Protein_Seq["I"]),0,Protein_Seq["I"])/SeqLength,5)
    AAC[11]=round(ifelse(is.na(Protein_Seq["L"]),0,Protein_Seq["L"])/SeqLength,5)
    AAC[12]=round(ifelse(is.na(Protein_Seq["K"]),0,Protein_Seq["K"])/SeqLength,5)
    AAC[13]=round(ifelse(is.na(Protein_Seq["M"]),0,Protein_Seq["M"])/SeqLength,5)
    AAC[14]=round(ifelse(is.na(Protein_Seq["F"]),0,Protein_Seq["F"])/SeqLength,5)
    AAC[15]=round(ifelse(is.na(Protein_Seq["P"]),0,Protein_Seq["P"])/SeqLength,5)
    AAC[16]=round(ifelse(is.na(Protein_Seq["S"]),0,Protein_Seq["S"])/SeqLength,5)
    AAC[17]=round(ifelse(is.na(Protein_Seq["T"]),0,Protein_Seq["T"])/SeqLength,5)
    AAC[18]=round(ifelse(is.na(Protein_Seq["W"]),0,Protein_Seq["W"])/SeqLength,5)
    AAC[19]=round(ifelse(is.na(Protein_Seq["Y"]),0,Protein_Seq["Y"])/SeqLength,5)
    AAC[20]=round(ifelse(is.na(Protein_Seq["V"]),0,Protein_Seq["V"])/SeqLength,5)
    AAC[21]=round(ifelse(is.na(Protein_Seq["U"]),0,Protein_Seq["U"])/SeqLength,5)
    AAC[22]=round(ifelse(is.na(Protein_Seq["X"]),0,Protein_Seq["X"])/SeqLength,5)
    if(Groups==2)
    {
      AAC2=c(NonPolar_Aliphatic=0,NonPolar_Aromatic=0,Polar_Uncharged=0,Polar_PositivelyCharged=0,Polar_NegativelyCharged=0)
      AAC2[1]=sum(AAC[c(1,20,13,10,11,8)])
      AAC2[2]=sum(AAC[c(18,19,14)])
      AAC2[3]=sum(AAC[c(16,17,15,5,3,7)])
      AAC2[4]=sum(AAC[c(2,9,12)])
      AAC2[5]=sum(AAC[c(4,6)])
      return(AAC2)
    }
    if(Groups==3)
    {
      AAC2=c(NonPolar_Aliphatic=0,NonPolar_Aromatic=0,Polar_Uncharged=0,Polar_PositivelyCharged=0,Polar_NegativelyCharged=0,Special=0)
      AAC2[1]=sum(AAC[c(1,20,13,10,11)])
      AAC2[2]=sum(AAC[c(18,19,14)])
      AAC2[3]=sum(AAC[c(16,17,3,7)])
      AAC2[4]=sum(AAC[c(2,9,12)])
      AAC2[5]=sum(AAC[c(4,6)])
      AAC2[6]=sum(AAC[c(5,8,15,21,22)])
      return(AAC2)
    }
    if(Groups==4)
    {
      AAC2=c(AE=0,ILFMV=0,NDTS=0,G=0,P=0,RKQH=0,YW=0,C=0)
      AAC2[1]=sum(AAC[c(1,6)])
      AAC2[2]=sum(AAC[c(10,11,13,14,20)])
      AAC2[3]=sum(AAC[c(16,17,3,4)])
      AAC2[4]=sum(AAC[c(8)])
      AAC2[5]=sum(AAC[c(15)])
      AAC2[6]=sum(AAC[c(2,12,7,9)])
      AAC2[7]=sum(AAC[c(18,19)])
      AAC2[8]=sum(AAC[c(5)])
      return(AAC2)
    }
    return(AAC)
  }else
  {
    stop("There is n't any sequence")
  }
}
