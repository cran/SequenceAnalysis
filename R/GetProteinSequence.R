#' @title SequenceAnalysis.GetProteinSequence
#' @description Get Protein Sequence From UniProt
#' @details Get Protein Sequence From UniProt by UniProt ID
#' @author Babak Khorsand
#' @import XML
#' @export SequenceAnalysis.GetProteinSequence
#' @param UniprotKB UniProt ID of a sequence
#' @return Protein Sequence
#' @examples
#' SequenceAnalysis.GetProteinSequence("O15131")
SequenceAnalysis.GetProteinSequence = function(UniprotKB)
{
  doc.text="N/A"
  url = paste(Protein_Address,UniprotKB,sep = "")
  doc.html = tryCatch({htmlTreeParse(url, useInternalNodes = TRUE)},error=function(err){return (NULL)})
  if (!is.null(doc.html))
  {
    doc.text = unlist(xpathApply(doc.html, '//pre', xmlValue))
    if (!is.null(doc.text))
    {
      doc.text = gsub("[0-9]","",doc.text)
      doc.text = gsub(" ","",doc.text)
      Protein=doc.text[1]
      if (is.null(Protein))
        Protein="N/A"
      names(Protein)="Protein"
      return(Protein)
    }
  }
}
