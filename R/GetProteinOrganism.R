#' @title SequenceAnalysis.GetProteinOrganism
#' @description Get Protein Organism
#' @details Get Protein, Gene and Organism of desired protein from UniProt
#' @author Babak Khorsand
#' @export SequenceAnalysis.GetProteinOrganism
#' @param UniprotKB UniProt ID of a sequence
#' @return Protein, Gene and Organism of desired UniProt ID.
#' @examples
#' SequenceAnalysis.GetProteinOrganism("P03485")
SequenceAnalysis.GetProteinOrganism = function(UniprotKB)
{
  if (is.null(UniprotKB))
    stop("Please Enter the uniprot id of your desired protein")
  else
  {
    url = paste(Protein_Address,UniprotKB,sep = "")
    mydata=tryCatch(readLines(url,n=5,warn=FALSE),error=function(err){print("Connection Error"); return(NULL)})
    if (!is.null(mydata))
    {
      Line=ifelse(length(grep("protein",mydata))>0,grep("protein",mydata),2)
      mydata = gsub(".*<title>(.*)</title>.*","\\1", mydata[Line])
      mydata=unlist(strsplit(mydata," - "))
      Organism=ifelse(length(mydata)>3,3,length(mydata))
      Result=c(Protein = mydata[Organism-1],Gene = ifelse(Organism>2,mydata[Organism-2],"N/A"),Organism=mydata[Organism])
    }
    return(Result)
  }
}
