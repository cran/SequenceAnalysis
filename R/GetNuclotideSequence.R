#' @title SequenceAnalysis.GetNucleotideSequence
#' @description Get Nucleotide Sequence
#' @details Get Nucleotide Sequence from EBI database by UniProt ID
#' @author Babak Khorsand
#' @import XML
#' @export SequenceAnalysis.GetNucleotideSequence
#' @param UniprotKB UniProt ID of a sequence
#' @return List The first element specify UniprotKB, The second element specify CDS place e.g. 64:1250, the third element specify nucleotide sequence of coding region and the fourth element return whole nucleotide sequence.
#' @examples
#' SequenceAnalysis.GetNucleotideSequence("O15131")
SequenceAnalysis.GetNucleotideSequence = function(UniprotKB)
{
  CDS=NULL
  Nucleotide_CDS=NULL
  Nucleotide_Sequence=NULL
  url = paste("http://www.uniprot.org/uniprot/",UniprotKB,sep = "")
  doc.html = tryCatch({htmlTreeParse(url, useInternalNodes = TRUE)},error=function(err){return (NULL)},warning=function(warn){})
  if (!is.null(doc.html))
  {
    doc.text = unlist(xpathApply(doc.html, '//td', xmlValue))
    if (!is.null(doc.text))
    {
      doc.text = doc.text[grep("EMBL nucleotide sequence",doc.text)[1]+1]
      if (length(doc.text)>0)
      {
        url2=paste("http://www.ebi.ac.uk/ena/data/view/",sub(" .*","", doc.text),"&display=text",sep = "")
        doc.text=tryCatch({readLines(url2)},error=function(err){return (NULL)},warning=function(warn){})
        if (!is.null(doc.text))
        {
          Nucleotide_Sequences=grep("^SQ",doc.text)
          if (length(Nucleotide_Sequences)>0)
          {
            Nucleotide_Sequence = doc.text[(Nucleotide_Sequences+1):(length(doc.text)-1)]
            Nucleotide_Sequence = gsub("[[:space:]]","",Nucleotide_Sequence)
            Nucleotide_Sequence = paste(Nucleotide_Sequence,sep = "",collapse = "")
            Nucleotide_Sequence = gsub("[0-9]","",Nucleotide_Sequence)
          }
          Uniprot=grep("^FT.*/db_xref=.*UniProtKB",doc.text)
          Uniprot=grep(UniprotKB,doc.text[Uniprot])
          if (length(Uniprot)>0)
          {
            CDS=grep("^FT.*CDS.*[0-9]*\\.\\.[0-9]*",doc.text)
            CDS_Line=CDS[Uniprot]
            while (is.na(CDS_Line))
            {
              Uniprot=Uniprot - 1
              CDS_Line=CDS[Uniprot]
            }
            CDS=doc.text[CDS_Line]
            for(i in 1:length(CDS))
            {
              while (substr(CDS[i],nchar(CDS[i]),nchar(CDS[i]))==",")
              {
                CDS[i]=paste(CDS[i],doc.text[CDS_Line[i]+1])
                CDS_Line[i]=CDS_Line[i]+1
              }
            }
            CDS=paste(CDS,collapse = ",")
            CDS=gsub("[[:space:]]","",CDS)
            CDS=gsub("CDS","",CDS)
            CDS=gsub("FT","",CDS)
            if (length(grep(":",CDS))==0)
            {
              CDS=gsub("\\.\\.",":",CDS)
              CDS2=gsub("<","",CDS)
              CDS2=gsub(">","",CDS)
              if (length(grep("^complement",CDS2))>0)
              {
                k=gsub("complement\\((.*)\\)","\\1",CDS2)
                k=strsplit(k,",")
                k=unlist(k)
                k=strsplit(k,":")
                k=unlist(k)
                k=gsub("[^0-9]","",k)
                k=as.integer(k)
                t=integer()
                for (i in seq(1,length(k),2))
                {
                  t=c(t,k[i]:k[i+1])
                }
                t=unique(t)
                CDS_Nuc_Char=strsplit(Nucleotide_Sequence,"")[[1]][t]
                g=rev(CDS_Nuc_Char)
                nuc_a=which(g=="a")
                nuc_t=which(g=="t")
                nuc_c=which(g=="c")
                nuc_g=which(g=="g")
                g[nuc_a]="t"
                g[nuc_t]="a"
                g[nuc_c]="g"
                g[nuc_g]="c"
                g=paste(g,collapse = "")
                Nucleotide_CDS=g
              }else
              {
                k=sub(",complement.*","",CDS2)
                k=strsplit(k,",")
                k=unlist(k)
                k=strsplit(k,":")
                k=unlist(k)
                k=gsub("[^0-9]","",k)
                k=as.integer(k)
                t=integer()
                for (i in seq(1,length(k),2))
                {
                  t=c(t,k[i]:k[i+1])
                }
                t=unique(t)
                Nucleotide_CDS=paste(strsplit(Nucleotide_Sequence,"")[[1]][t],collapse = "")
              }
            }
          }
        }
      }
    }
  }
  if (is.null(CDS))
  {
    CDS="..."
  }else if (length(CDS)==0)
  {
    CDS="..."
  }
  if (is.null(Nucleotide_CDS))
  {
    Nucleotide_CDS="..."
  }else if (length(CDS)==0){
    Nucleotide_CDS="..."
  }
  if (is.null(Nucleotide_Sequence))
  {
    Nucleotide_Sequence="..."
  }else if (length(Nucleotide_Sequence)==0)
  {
    Nucleotide_Sequence="..."
  }
  returnList = list(UniprotKB=UniprotKB,CDS=CDS,
                    Nucleotide_CDS=Nucleotide_CDS,
                    Nucleotide_Sequence=Nucleotide_Sequence)
  return(returnList)
}
