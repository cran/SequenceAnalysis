#' @title SequenceAnalysis.GetNucleotideSequence
#' @description Get Nucleotide Sequence
#' @details Get Nucleotide Sequence from EBI database by UniProt ID
#' @author Babak Khorsand
#' @import XML
#' @export SequenceAnalysis.GetNucleotideSequence
#' @param UniprotKB UniProt ID of a sequence
#' @return List The first element specify UniprotKB, The second element specify nucleotide sequence of coding region and the third element return whole nucleotide sequence.
#' @examples
#' SequenceAnalysis.GetNucleotideSequence("O15131")
SequenceAnalysis.GetNucleotideSequence = function(UniprotKB)
{
  Nucleotide_CDS=NULL
  Nucleotide_Sequence=NULL
  url = paste(Protein_Address,UniprotKB,sep = "")
  doc.html = tryCatch({htmlTreeParse(url, useInternalNodes = TRUE)},error=function(err){return (NULL)})
  if (!is.null(doc.html))
  {
    doc.text = unlist(xpathApply(doc.html, '//td', xmlValue))
    if (!is.null(doc.text))
    {
      doc.text = doc.text[grep(Nuclotide_Hint,doc.text)[1]+1]
      if (length(doc.text)>0)
      {
        Translation=sub("Translation: ([^\\.]*).*","\\1",doc.text)
        if (length(Translation)>0)
        {
          Translation=unlist(strsplit(sub("Translation: ([^\\.]*).*","\\1",doc.text)," "))
          Translation=Translation[length(Translation)]
          url2=paste(Nucleotide_Address,Translation,"&display=text",sep = "")
          doc.translation=tryCatch({readLines(url2,warn=FALSE)},error=function(err){return (NULL)})
          if (!is.null(doc.translation))
          {
            Nucleotide_CDS=grep("^SQ",doc.translation)
            if (length(Nucleotide_CDS)>0)
            {
              Nucleotide_CDS = doc.translation[(Nucleotide_CDS+1):(length(doc.translation)-1)]
              Nucleotide_CDS = gsub("[[:space:]]","",Nucleotide_CDS)
              Nucleotide_CDS = paste(Nucleotide_CDS,sep = "",collapse = "")
              Nucleotide_CDS = gsub("[0-9]","",Nucleotide_CDS)
            }else
            {
              Nucleotide_CDS="N/A"
            }
          }
        }
        url2=paste(Nucleotide_Address,sub(" .*","", doc.text),"&display=text",sep = "")
        doc.text=tryCatch({readLines(url2,warn=FALSE)},error=function(err){return (NULL)})
        if (!is.null(doc.text))
        {
          Nucleotide_Sequences=grep("^SQ",doc.text)
          if (length(Nucleotide_Sequences)>0)
          {
            Nucleotide_Sequence = doc.text[(Nucleotide_Sequences+1):(length(doc.text)-1)]
            Nucleotide_Sequence = gsub("[[:space:]]","",Nucleotide_Sequence)
            Nucleotide_Sequence = paste(Nucleotide_Sequence,sep = "",collapse = "")
            Nucleotide_Sequence = gsub("[0-9]","",Nucleotide_Sequence)
          }else
          {
            Nucleotide_Sequence="N/A"
          }
        }
      }
    }
  }
  if (is.null(Nucleotide_CDS))
  {
    Nucleotide_CDS="N/A"
  }
  if (is.null(Nucleotide_Sequence))
  {
    Nucleotide_Sequence="N/A"
  }
  returnList = list(UniprotKB=UniprotKB,
                    Nucleotide_CDS=Nucleotide_CDS,
                    Nucleotide_Sequence=Nucleotide_Sequence)
  return(returnList)
}
