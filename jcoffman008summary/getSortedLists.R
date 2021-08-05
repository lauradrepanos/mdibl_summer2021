#Getting lists for jcoffman008 downstream analysis. Alternative to sorting in Excel, allows easy updates if change of mind in cutoffs/ procedure + tracks filtering for reproducibility  
#USAGE: Rscript getSortedLists.R
# edit SE and baseMean cutoffs and filepaths in input_getSortedLists.json 
# output in sorted_lists directory 

options(stringsAsFactors = FALSE)
library(jsonlite)
library(tidyr) #needed for separate()
library(ggplot2)

#create subdirectories for output 
dir.create("sorted_lists",showWarnings=FALSE)
dir.create("sorted_lists/histograms",showWarnings=FALSE)

json = read_json("input_getSortedLists.json")

##establish parameters 
lfcSE_cutoff_pct <- json$lfcSE_cutoff_pct #genes in the top __ % of the distribution of lfcSE are excluded from downstream analysis
baseMean_cutoff_pct <- json$baseMean_cutoff_pct #genes in the bottom __% of the distribution of baseMean values are excluded from downstream analysis


print("loading files...")
#DESEQ that uses PC1 and time as covariant on all KO vs WT samples (negative fold change = downregulated in KOs)
DESeq2_out <- read.delim(json$filepath$deseq)
#DESEQ on just samples collected at 3hours before lights on KO vs WT with PC1 (from PCA on just 3hb samples) as covariant
DESeq2_out_3hb <- read.delim(json$filepath$deseq3hb)
DESeq2_out_lightson <- read.delim(json$filepath$deseqlightson)
DESeq2_out_2.5ha <- read.delim(json$filepath$deseq2.5ha)
#loading scores from PCA on all samples 
PCA <- read.delim(json$filepath$loading_scores)

##filter files: make low count, p-value cutoffs and separate gene name and ID (since downstream analysis tools do not take GeneID_GeneName format)
print("filtering out low counts and p values...")
DESeq2_out<- separate(data = DESeq2_out, col = X, into = c("GeneID", "GeneName"), sep = "\\_",extra="merge",fill="right")
baseMean_cutoff<- quantile(DESeq2_out$baseMean,probs=baseMean_cutoff_pct)
lfcSE_cutoff <- quantile(DESeq2_out$lfcSE,probs=(1-lfcSE_cutoff_pct))
#plots that demonstrate distribution and cutoff of lfcSE and baseMean
png("sorted_lists/histograms/DESEQ_PC1_time_se.png")
ggplot()+
  geom_histogram(aes(log(DESeq2_out$lfcSE)), binwidth = 0.1, col ="black", fill = "white")+
  geom_vline(xintercept =  log(lfcSE_cutoff), linetype = "dashed", col = 'blue')+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  labs(title = "Distribution of lfcSE of genes in DESeq of all KO vs. WT",
       caption = "The blue line represents the lfcSE threshold above which genes are excluded from downstream analysis.")+
  xlab("Standard error of log fold change (log scale)")+
  ylab("gene count")
dev.off()
png("sorted_lists/histograms/DESEQ_PC1_time_baseMean.png")
ggplot()+
  geom_histogram(aes(log(DESeq2_out$baseMean)), binwidth = 0.1, col ="black", fill = "white")+
  geom_vline(xintercept =  log(baseMean_cutoff), linetype = "dashed", col = 'blue')+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  labs(title = "Distribution of baseMean of genes in all samples",
       caption = "The blue line represents the baseMean threshold below which genes are excluded from downstream analysis.")+
  xlab("baseMean transcript count (log scale)")+
  ylab("gene count")
dev.off()
DESeq2_out <- subset(DESeq2_out,baseMean>baseMean_cutoff)
DESeq2_out <- subset(DESeq2_out,lfcSE<lfcSE_cutoff)


DESeq2_out_3hb<- separate(data = DESeq2_out_3hb, col = X, into = c("GeneID", "GeneName"), sep = "\\_",extra="merge",fill="right")
baseMean_cutoff_3hb<- quantile(DESeq2_out_3hb$baseMean,probs=baseMean_cutoff_pct)
lfcSE_cutoff_3hb <- quantile(DESeq2_out_3hb$lfcSE,probs=(1-lfcSE_cutoff_pct))
png("sorted_lists/histograms/DESEQ_3hb_se.png")
ggplot()+
  geom_histogram(aes(log(DESeq2_out_3hb$lfcSE)), binwidth = 0.1, col ="black", fill = "white")+
  geom_vline(xintercept =  log(lfcSE_cutoff_3hb), linetype = "dashed", col = 'blue')+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  labs(title = "Distribution of lfcSE of genes in DESeq of KO vs. WT 3h before lights on",
       caption = "The blue line represents the lfcSE threshold above which genes are excluded from downstream analysis.")+
  xlab("Standard error of log fold change (log scale)")+
  ylab("gene count")
dev.off()
png("sorted_lists/histograms/DESEQ_3hb_baseMean.png")
ggplot()+
  geom_histogram(aes(log(DESeq2_out_3hb$baseMean)), binwidth = 0.1, col ="black", fill = "white")+
  geom_vline(xintercept =  log(baseMean_cutoff_3hb), linetype = "dashed", col = 'blue')+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  labs(title = "Distribution of baseMean of genes 3h before lights on",
       caption = "The blue line represents the baseMean threshold below which genes are excluded from downstream analysis.")+
  xlab("baseMean transcript count (log scale)")+
  ylab("gene count")
dev.off()
DESeq2_out_3hb <- subset(DESeq2_out_3hb,baseMean>baseMean_cutoff_3hb)
DESeq2_out_3hb <- subset(DESeq2_out_3hb,lfcSE<lfcSE_cutoff_3hb)

DESeq2_out_lightson<- separate(data = DESeq2_out_lightson, col = X, into = c("GeneID", "GeneName"), sep = "\\_",extra="merge",fill="right")
baseMean_cutoff_lightson<- quantile(DESeq2_out_lightson$baseMean,probs=baseMean_cutoff_pct)
lfcSE_cutoff_lightson <- quantile(DESeq2_out_lightson$lfcSE,probs=(1-lfcSE_cutoff_pct))
png("sorted_lists/histograms/DESEQ_lightson_se.png")
ggplot()+
  geom_histogram(aes(log(DESeq2_out_lightson$lfcSE)), binwidth = 0.1, col ="black", fill = "white")+
  geom_vline(xintercept =  log(lfcSE_cutoff_lightson), linetype = "dashed", col = 'blue')+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  labs(title = "Distribution of lfcSE of genes in DESeq of KO vs. WT at lights on",
       caption = "The blue line represents the lfcSE threshold above which genes are excluded from downstream analysis.")+
  xlab("Standard error of log fold change (log scale)")+
  ylab("gene count")
dev.off()
png("sorted_lists/histograms/DESEQ_lightson_baseMean.png")
ggplot()+
  geom_histogram(aes(log(DESeq2_out_lightson$baseMean)), binwidth = 0.1, col ="black", fill = "white")+
  geom_vline(xintercept =  log(baseMean_cutoff_lightson), linetype = "dashed", col = 'blue')+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  labs(title = "Distribution of baseMean of genes at lights on",
       caption = "The blue line represents the baseMean threshold below which genes are excluded from downstream analysis.")+
  xlab("baseMean transcript count (log scale)")+
  ylab("gene count")
dev.off()
DESeq2_out_lightson <- subset(DESeq2_out_lightson,baseMean>baseMean_cutoff_lightson)
DESeq2_out_lightson <- subset(DESeq2_out_lightson,lfcSE<lfcSE_cutoff_lightson)

DESeq2_out_2.5ha<- separate(data = DESeq2_out_2.5ha, col = X, into = c("GeneID", "GeneName"), sep = "\\_",extra="merge",fill="right")
baseMean_cutoff_2.5ha<- quantile(DESeq2_out_2.5ha$baseMean,probs=baseMean_cutoff_pct)
lfcSE_cutoff_2.5ha <- quantile(DESeq2_out_2.5ha$lfcSE,probs=(1-lfcSE_cutoff_pct))
png("sorted_lists/histograms/DESEQ_2.5ha_se.png")
ggplot()+
  geom_histogram(aes(log(DESeq2_out_2.5ha$lfcSE)), binwidth = 0.1, col ="black", fill = "white")+
  geom_vline(xintercept =  log(lfcSE_cutoff_2.5ha), linetype = "dashed", col = 'blue')+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  labs(title = "Distribution of lfcSE of genes in DESeq of KO vs. WT 2.5h after lights on",
       caption = "The blue line represents the lfcSE threshold above which genes are excluded from downstream analysis.")+
  xlab("Standard error of log fold change (log scale)")+
  ylab("gene count")
dev.off()
png("sorted_lists/histograms/DESEQ_2.5ha_baseMean.png")
ggplot()+
  geom_histogram(aes(log(DESeq2_out_2.5ha$baseMean)), binwidth = 0.1, col ="black", fill = "white")+
  geom_vline(xintercept =  log(baseMean_cutoff_2.5ha), linetype = "dashed", col = 'blue')+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  labs(title = "Distribution of baseMean of genes 2.5 hours after lights on",
       caption = "The blue line represents the baseMean threshold below which genes are excluded from downstream analysis.")+
  xlab("baseMean transcript count (log scale)")+
  ylab("gene count")
dev.off()
DESeq2_out_2.5ha <- subset(DESeq2_out_2.5ha,baseMean>baseMean_cutoff_2.5ha )
DESeq2_out_2.5ha <- subset(DESeq2_out_2.5ha,lfcSE<lfcSE_cutoff_2.5ha)

PCA<- separate(data = PCA, col = X, into = c("GeneID", "GeneName"), sep = "\\_",extra="merge",fill="right")
#do not need to filter out low counts since rpca already does this 


print("making sorted lists...")

DESeq2_out_upreg <- DESeq2_out[order(-DESeq2_out$log2FoldChange),]
DESeq2_out_downreg <- DESeq2_out[order(DESeq2_out$log2FoldChange),]
#sidenote: uncommenting line below allows you to see the stats of one gene: 
#print(DESeq2_out[which(DESeq2_out$GeneName=="hmgcra"),])

DESeq2_out_3hb_upreg <- DESeq2_out_3hb[order(-DESeq2_out_3hb$log2FoldChange),]
DESeq2_out_3hb_downreg <- DESeq2_out_3hb[order(DESeq2_out_3hb$log2FoldChange),]

DESeq2_out_lightson_upreg <- DESeq2_out_lightson[order(-DESeq2_out_lightson$log2FoldChange),]
DESeq2_out_lightson_downreg <- DESeq2_out_lightson[order(DESeq2_out_lightson$log2FoldChange),]

DESeq2_out_2.5ha_upreg <- DESeq2_out_2.5ha[order(-DESeq2_out_2.5ha$log2FoldChange),]
DESeq2_out_2.5ha_downreg <- DESeq2_out_2.5ha[order(DESeq2_out_2.5ha$log2FoldChange),]

PC3_ascending <- PCA[order(PCA$PC3),]
PC3_descending <- PCA[order(-PCA$PC3),]

print("generating output files and placing in sorted_lists directory...")


write.csv(DESeq2_out_upreg,"sorted_lists/jcoffman008_DESEQ_PC1_time_upreg_in_KO.csv",row.names=FALSE)
write.csv(DESeq2_out_downreg,"sorted_lists/jcoffman008_DESEQ_PC1_time_downreg_in_KO.csv",row.names=FALSE)


write.csv(DESeq2_out_3hb_upreg,"sorted_lists/jcoffman008_DESEQ_3hb_upreg_in_KO.csv",row.names=FALSE)
write.csv(DESeq2_out_3hb_downreg,"sorted_lists/jcoffman008_DESEQ_3hb_downreg_in_KO.csv",row.names=FALSE)


write.csv(DESeq2_out_lightson_upreg,"sorted_lists/jcoffman008_DESEQ_lightson_upreg_in_KO.csv",row.names=FALSE)
write.csv(DESeq2_out_lightson_downreg,"sorted_lists/jcoffman008_DESEQ_lightson_downreg_in_KO.csv",row.names=FALSE)


write.csv(DESeq2_out_2.5ha_upreg,"sorted_lists/jcoffman008_DESEQ_2.5ha_upreg_in_KO.csv",row.names=FALSE)
write.csv(DESeq2_out_2.5ha_downreg,"sorted_lists/jcoffman008_DESEQ_2.5ha_downreg_in_KO.csv",row.names=FALSE)

write.csv(PC3_ascending,"sorted_lists/jcoffman008_sortbyPC3_downreg_in_KO.csv",row.names=FALSE)
write.csv(PC3_descending,"sorted_lists/jcoffman008_sortbyPC3_upreg_in_KO.csv",row.names=FALSE)

