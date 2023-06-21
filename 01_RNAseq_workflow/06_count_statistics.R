args   <- commandArgs(TRUE)
setwd(args[1])
#setwd("/Users/admin/Analysis/17_20191202_Mm_RNAseq/counts/")


library(RColorBrewer)
library(ggplot2)
palette(unique(c(brewer.pal(9,"Set1"),brewer.pal(8,"Set2"))))

library(DESeq2)

count.table<-read.delim("count_table.txt", skip=1)
names1<-sapply(strsplit(colnames(count.table[,7:ncol(count.table)]),
                        split = ".", fixed = T), function(x) 
                            paste(x[(grep("bams", x)+1)]))
colnames(count.table)[7:ncol(count.table)]<-names1

######################### get gene Annotation ###################################
library(biomaRt)
listMarts()
ensembl <- useMart("ensembl")
ensembl = useDataset("mmusculus_gene_ensembl",mart=ensembl)

musGenes<-as.character(count.table$Geneid)

searchAttributes(mart = ensembl, pattern = "ensembl")


biotype.raw<-getBM(attributes = c("ensembl_gene_id","gene_biotype","mgi_symbol"),
      filters="ensembl_gene_id",values=musGenes ,mart=ensembl)
biotype<-biotype.raw[(!duplicated(biotype.raw$ensembl_gene_id)),]

write.table(biotype,file = "~/Annotation/MusMusculus/protein_coding_and_lincRNA_INFO.txt",
            sep="\t", quote=F, row.names = F)

#biotype<-read.delim("~/Annotation/MusMusculus/protein_coding_and_lincRNA_INFO.txt")

lincRNA.table<-subset(biotype, gene_biotype=="lncRNA")
lincRNA<-lincRNA.table$ensembl_gene_id

ProtCode.table<-subset(biotype, gene_biotype=="protein_coding")
ProtCode<-as.character(ProtCode.table$ensembl_gene_id)

##################################################################################
counts<-count.table[,7:ncol(count.table)]

############## histograms ######################
dir.create("plots")
setwd("plots/")

max.y<-max(apply(counts, 2, function(x) sum(log2(x+1)==0)))

pdf("count_histograms_withZeros.pdf", width=5, height =5 )
    for( i in 1:ncol(counts)){
        hist(log2(counts[,i]+1), breaks = 100, col="lightsalmon",
        ylim=c(0,max.y), main=colnames(counts)[i], xlab="log2(counts+1)")
    }
dev.off()


################



############## plot ############

max.x<-c()
max.y<-c()
for(i in 1:ncol(counts)){
  max.x<-c(max.x,max(density(log2(counts[,i]))$x))
  max.y<-c(max.y,max(density(log2(counts[,i]))$y))
}

pdf("count_density_raw.pdf", width=7, height = 5)
  par(mar=c(5.1, 4.1, 4.1, 12.1))
  plot(NULL,NULL,xlab="log2(counts)", ylab="Density", xlim=c(-2,max(max.x)), ylim=c(0,max(max.y)),
       main="protein coding and lincRNA")
  
  for(i in 1:ncol(counts)){
    lines(density(log2(counts[,i])), col=i, lwd=2.5)
  }
  legend(max(max.x)+1,max(max.y),legend =colnames(counts), xpd = T,  col=1:ncol(counts),
         bty="n", lwd = 2, cex = 0.75)
  #legend("topright", legend = colnames(counts),  col=1:ncol(counts),
   #      bty="n", lwd = 2, cex = 0.75 )

dev.off()



################## count distribution ##########################

#plot count distribution
totalCounts<-apply(counts,2,sum)
#get max for plot
max.x<-c()
max.y<-c()
for(i in 1:ncol(counts)){
  max.x<-c(max.x,max(density(log2(counts[,i]/(totalCounts[i]/1e6)))$x))
  max.y<-c(max.y,max(density(log2(counts[,i]/(totalCounts[i]/1e6)))$y))
}



pdf("count_density_CPM.pdf", width=7, height = 5)
  par(mar=c(5.1, 4.1, 4.1, 12.1))
  plot(NULL,NULL,xlab="log2(counts per million)", ylab="Density", xlim=c(-7,max(max.x)), ylim=c(0,max(max.y)),
       main="protein coding and lincRNA")
  
  for(i in 1:ncol(counts)){
    lines(density(log2(counts[,i]/(totalCounts[i]/1e6))), col=i, lwd=2.5)
  }
  
   legend(max(max.x)+1,max(max.y),legend =colnames(counts), xpd = T,  col=1:ncol(counts),
           bty="n", lwd = 2, cex = 0.75)
  abline(v=0, lty=2)
dev.off()


#################### detection ##############################

# plot feature detection
normCounts<-t(t(counts)/(totalCounts/1e6))
rownames(normCounts)<-as.character(count.table$Geneid)

#######################################################################################

#in at least 3 samples >1 CPM
table(apply(log2(normCounts),1,function(x) sum(x>0)>6))
exp.table<-(apply(log2(normCounts),1,function(x) sum(x>0)>6))

table((names(exp.table) %in% as.character(lincRNA)))


#######################################################################################

# at least 1 CPM (log2(1 CPM)=0)
lincRNA_det<-table(log2(rowMeans(normCounts[lincRNA,])) > 0)
lincRNA_rate<-lincRNA_det["TRUE"]/sum(lincRNA_det)



ProtCode_det<-table(log2(rowMeans(normCounts[ProtCode,])) > 0)
ProtCode_rate<-ProtCode_det["TRUE"]/sum(ProtCode_det)

df.BARplot<-data.frame(BioTypes=c("protein coding","lincRNA"),
                       CoveredFeatures=c(ProtCode_rate,lincRNA_rate), 
                       detected=c(ProtCode_det["TRUE"], lincRNA_det["TRUE"]))


g<-ggplot(data = df.BARplot, aes(x=BioTypes, y=CoveredFeatures))+
  geom_bar(stat="identity", fill="steelblue")+theme_minimal()+
  geom_text(aes(label=detected), color="black", size=6, vjust=-0.75)+
  ylab("% Covered Features (min. 1 CPM)") + ylim(0,1)+
  theme(axis.text=element_text(size=14), axis.title = element_text(size = 18),
        axis.title.y = element_text(margin=unit(c(0, 7, 0, 0), "mm")),
        axis.title.x = element_text(margin=unit(c(7, 0, 0, 0), "mm")),
        legend.key.size =unit(10,"mm"), legend.text = element_text(size=14),
        legend.title = element_text(size=18))

pdf("covered_Features.pdf", width=5, height = 5)
  print(g)
dev.off()

sink("covered_Features.txt")
  print(df.BARplot)
sink()

###################################
df.stats<-read.delim("../count_table.txt.summary")
rownames(df.stats)<-df.stats$Status


df.assignments<-data.frame(sample=names1, assigned_reads=as.numeric(df.stats["Assigned",-1]))

g<-ggplot(data = df.assignments, aes(x=sample, y=(assigned_reads/1e6)))+
  geom_bar(stat="identity",color="black", fill="slategrey", position="dodge")+theme_bw()+
  ylab("assigned reads (Mio)")+xlab("")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

pdf("stats_assigned_reads.pdf", width=(2+nrow(df.assignments)/4), height=6)
  print(g)
dev.off()

################### stacked plot #################

df.red<-df.stats[df.stats[,2]!=0,-1]
colnames(df.red)<-names1
percentages<-t(df.red)/apply(df.red,2,sum)

library(reshape2)

df.plot<-melt(percentages)
colnames(df.plot)<-c("sample", "assignment", "proportion")
df.plot$assignment<-factor(df.plot$assignment, levels = 
                             rev(levels(df.plot$assignment)))

pdf("stats_assignment.pdf",height=(nrow(df.assignments)/4), width=8)
  ggplot(df.plot, aes(x=sample,y=proportion, fill=assignment))+
  geom_bar(stat = "identity", position = "stack")+ylim(0,1)+
  xlab("")+theme_bw()+coord_flip()
dev.off()
