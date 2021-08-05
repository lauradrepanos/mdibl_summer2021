#PURPOSE: generates a PNG summarizing the time-dependent effect of the klf9-knockout upon genes in a user-defined GO group in jcoffman008 
#USAGE (sample GO): Rscript jcoffman008GOSummary.R "GO:0007601" 
      # add tag "zscore" after GO term if you wish to see the Z score of transcript counts of genes in this GO rather than log fold change 
      ##technically any tag after the GO term that isn't "log2fc" will result in the script showing z score 
#edit input_jcoffman008GOsummary.json to change filepaths

options(stringsAsFactors = FALSE)

library(jsonlite)
library(tidyr)
library(ggplot2)
library(reshape2)


args = base::commandArgs(trailingOnly = TRUE)
if(length(args)<1){
  myGO<- "GO:0007601" #default GO group to summarize
}else{
  myGO<- args[1]
}
if(length(args)==2){
  yaxis<- args[2] #will actually result in giving the zscore if this is anything but "log2fc" 
}else{
  yaxis <- "log2fc" #change to true if you wish to see gene expression instead of log2foldchange
}
filepath = read_json("input_jcoffman008GOSummary.json")


allGOterms<- read.delim(file=filepath$ens_GOterms)
genesInGO<- subset(allGOterms,GO.term.accession== myGO)
if(nrow(genesInGO)==0){
  error_message1<- paste("\n***The GO you entered, ", myGO, ", is not recognized by the ensembl reference file. Please make sure you wrote it in the right format (i.e. GO:0007601 or enter another GO term\n",sep="")
  stop(error_message1)
}
myGOterm <- genesInGO[1,]$GO.term.name
GeneIDsinGO<- genesInGO$Gene.stable.ID

if(yaxis=="log2fc"){
  dsq.3hb <- read.table(file=filepath$deseq3hb, sep = "\t",header=TRUE)
  dsq.lightson<- read.table(file=filepath$deseqlightson, sep = "\t",header=TRUE)
  dsq.2.5ha<- read.table(file=filepath$deseq2.5ha, sep = "\t",header=TRUE)

  dsq.3hb<- dsq.3hb[,c("X","log2FoldChange")]
  dsq.lightson<- dsq.lightson[,c("X","log2FoldChange")]
  dsq.2.5ha<- dsq.2.5ha[,c("X","log2FoldChange")]
  names(dsq.3hb)[names(dsq.3hb) == "log2FoldChange"] <- "3h before lights on"
  names(dsq.lightson)[names(dsq.lightson) == "log2FoldChange"] <- "lights on"
  names(dsq.2.5ha)[names(dsq.2.5ha) == "log2FoldChange"] <- "2.5 hours after lights on"
  dsq.merge1 <-merge(dsq.3hb,dsq.lightson)
  dsq.merge2<- merge(dsq.merge1,dsq.2.5ha)
  long_dsq<- melt(dsq.merge2)
  names(long_dsq)[names(long_dsq)=="variable"]<- "time"
  names(long_dsq)[names(long_dsq)=="value"]<- "log2FoldChange"
  long_dsq<- separate(data = long_dsq, col = X, into = c("GeneID", "GeneName"), sep = "\\_")
  
  lfcInGO<- long_dsq[long_dsq$GeneID %in% GeneIDsinGO,]
  if(nrow(lfcInGO)<5){
    error_message2<- paste('\n***The GO you entered, ', myGO, ', has too few genes sufficiently expressed in jcoffman008 to be represented in the DESEQ analysis output and summarized for statistical trends. Try adding the "zscore" tag to see Z-normalized transcript counts of genes in this GO instead***\n\n',sep='')
    stop(error_message2)
  }
  dir.create("GO_summaries",showWarnings=FALSE)
  png(filename=paste("GO_summaries/",myGOterm," jcoffman008.png",sep=""))

  ggplot(data=lfcInGO,aes(x=time,y=log2FoldChange,fill=time))+ scale_fill_manual(values=c("steelblue3","tomato1","mediumseagreen"))+
    geom_violin()+ geom_boxplot(width=0.25,alpha=0.2)+
    geom_hline(yintercept=0,linetype="dashed")+
    theme_bw()+
    labs( x="",y="log2FoldChange WT vs. klf9 -/-")+
    theme(legend.position="none", aspect.ratio = .75,
          axis.ticks.x=element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    ggtitle(paste("Gene expression in ",myGO,": ",myGOterm,sep=""))
  
}else{
  zscore <- read.delim(file=filepath$Z_normalized)
  zscore<- separate(data = zscore, col = X, into = c("GeneID", "GeneName"), sep = "\\_",extra="merge")
  zscoreInGO_wide<- zscore[zscore$GeneID %in% GeneIDsinGO,]
  zscoreInGO <- melt(zscoreInGO_wide)
  names(zscoreInGO)[names(zscoreInGO)=="variable"]<- "sample"
  names(zscoreInGO)[names(zscoreInGO)=="value"]<- "Zscore"
  zscoreInGO<- separate(data = zscoreInGO, col = sample, into = c("genetics", "time", "replicate"), sep = "\\_")
  zscoreInGO$time[zscoreInGO$time=="m30"]<--3
  zscoreInGO$time[zscoreInGO$time=="z00"]<-0
  zscoreInGO$time[zscoreInGO$time=="p25"]<-2.5
  zscoreInGO$genetics[zscoreInGO$genetics=="ko"]<-"klf9-/-"
  zscoreInGO$genetics<- factor(zscoreInGO$genetics,levels=c("wt","klf9-/-"))
  
  if(nrow(zscoreInGO)==0){
    error_message3<- paste('\n***The GO you entered, ', myGO, ', contains no genes sufficiently expressed in jcoffman008 to calculate normalized counts ***\n\n',sep='')
    stop(error_message3)
  }

  dir.create("GO_summaries",showWarnings=FALSE)
  png(filename=paste("GO_summaries/",myGOterm," znorm jcoffman008.png",sep=""))

  if(nrow(zscoreInGO_wide)<35){
    ggplot(data=zscoreInGO,aes(x=GeneName,y=Zscore))+ 
      geom_point(data= zscoreInGO[c(1,2,4,5,6)],color="grey85")+
      geom_point(aes(color=genetics))+ scale_color_manual(values = c("black", "red")) +
      facet_grid(.~ time+ genetics, scales='free',labeller=label_both)+
      theme_bw()+
      labs(x="gene", y="Z-normalized count")+
      theme(strip.background = element_blank(),
            legend.position="none",
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=7),axis.ticks.x=element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())+
      ggtitle(paste("Expression of genes in ",myGO,": ",myGOterm,sep=""))
  }else{
    ggplot(data=zscoreInGO,aes(x=GeneName,y=Zscore))+ 
      geom_point(data= zscoreInGO[c(1,2,4,5,6)],color="grey85")+
      geom_point(aes(color=genetics))+ scale_color_manual(values = c("black", "red")) +
      facet_grid(.~ time+ genetics, scales='free',labeller=label_both)+
      theme_bw()+
      labs(x="gene", y="Z-normalized count")+
      theme(strip.background = element_blank(),
            legend.position="none",
            axis.text.x = element_blank(),axis.ticks.x=element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())+
      ggtitle(paste("Expression of genes in ",myGO,": ",myGOterm,sep=""))
  }
}

dev.off() 



          