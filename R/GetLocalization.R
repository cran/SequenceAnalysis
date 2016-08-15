#' @title SequenceAnalysis.GetLocalization
#' @description Subcellular location of that protein
#' @details Getting Localization From Gene Ontology Cellular Component inside Uniprot.
#' @author Babak Khorsand
#' @export SequenceAnalysis.GetLocalization
#' @param UniprotKB Uniprot ID of a sequence
#' @return Subcellular location of that protein
#' @examples
#' SequenceAnalysis.GetLocalization("O15131")
SequenceAnalysis.GetLocalization = function(UniprotKB)
{
  Localization=NULL
  url = paste(Protein_Address,UniprotKB,sep = "")
  doc.text = tryCatch({readLines(url,warn=FALSE)},error=function(err){return (NULL)})
  if (!is.null(doc.text))
  {
    doc.text = doc.text[grep(Localization_Address,doc.text)]
    if (length(doc.text)>0)
    {
      Localization = character()
      pattern=".*>(.*?)</a>.*"
      Localization = gsub(pattern,"\\1",doc.text)
      while (length(grep("\\.\\.\\.",Localization))==0)
      {
        pattern=paste(".*>.*?</a>",pattern,sep = "")
        Localization = c(Localization,gsub(pattern,"\\1",doc.text))
      }
      Localization=Localization[-length(Localization)]
    }
  }
  if (is.null(Localization))
    Localization="N/A"
  return(Localization)
}
