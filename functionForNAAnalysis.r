# Functions for missing data analysis

produceHtmlMissingByBlock=function(dfList,selection)
{
  html="<html>"
  for(df in 1:length(dfList))
  {
    visits=c("V01","V05","V09","V29","V30")
    listBlocks=list()
    for(visit in visits) #V13 pour un patient
    {
      print(visit)
      dataVisit=dfList[[df]][dfList[[df]][,"Visit"]==visit,]
      listBlocks[[visit]]=removeDuplicatedRownames(df=dataVisit,varToSelect=selection[[df]])
      html=paste0(html,"<h3>",visit,"</h3><p><img src='",visit,".png'></p>")
    }
    p=get_patternNA(listBlocks)
    png(file=paste0(names(dfList)[df],".png"),width=800,height=400)
    plot(p)
    dev.off()
    html=paste0(html,"<h3>",names(dfList)[df],"</h3><p><img src='",names(dfList)[df],".png'></p>")
    
  }
  write(html,file="res.html")
  
  
}


produceHtmlMissingByVisit=function()
{
  html="<html>"
  visits=c("V01","V05","V09","V29","V30")
  for(visit in visits) #V13 pour un patient
  {
    print(visit)
    cym_V01=df_cytometrie[df_cytometrie[,"Visit"]==visit,]
    cyk_V01=df_cytokine[df_cytokine[,"Visit"]==visit,];
    cym_V01_f=removeDuplicatedRownames(df=cym_V01,varToSelect=varnum_cytometrie)
    cyk_V01_f=removeDuplicatedRownames(df=cyk_V01,varToSelect=varnum_cytokine)
    p=get_patternNA(list(cytometrie=cym_V01_f,cytokine=cyk_V01_f))
    
    png(file=paste0(visit,".png"),width=800,height=400)
    plot(p)  
    dev.off()
    html=paste0(html,"<h3>",visit,"</h3><p><img src='",visit,".png'></p>")
  }
  write(html,file="res.html")
}



removeDuplicatedRownames=function(df,varToSelect)
{
  indToSup=which(duplicated(df[,"Subject"]))  
  if(length(indToSup)!=0)
  { print("Warning ! duplicated subject names ! ! ")
    print(df[indToSup,"Subject"])
    df2=df[-indToSup,varToSelect ];rownames(df2)=df[-indToSup,"Subject"]
  }
  else
  {
    df2=df[,varToSelect ];rownames(df2)=df[,"Subject"]
  }
  return(df2)
}

