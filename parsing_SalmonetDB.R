library(stringr)
library(RCurl)
library(XML)

getInteraction <- function(prot){
  # prot <- "e1w6v8"
  doc <- htmlParse(paste0("http://salmonet.org/protein/",prot,"/"))
  doc_content <- xpathSApply(doc, "////input/@value")
  doc_content_useful <- str_split(doc_content[[1]],"\n")
  tempList <- str_split(doc_content_useful[[1]],",")
  tempoutput <- do.call(rbind, lapply(tempList, FUN = function(x){
         return(x[c(7,8)])
     }))
  return(tempoutput)
}

protList <- c("e1w6v8", "e1w6w2", "e1w6w7")
mat <- do.call(rbind, sapply(protList, FUN = getInteraction))
