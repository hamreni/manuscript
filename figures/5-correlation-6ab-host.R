library(ggplot2)
library(DESeq2)
library(tidyverse)
library(rtracklayer)
library(FinCal)
library(cubature)
library(miscTools)
library(psych)
library(vsn)
library(ggrepel)
library(RColorBrewer)

snodata <- read.csv("data/sno_genes.csv", sep=";")
normdata <- read.csv(file = "data/all_norm_total_counts.csv", sep=";")
interest <- normdata %>% filter(row.names(normdata) %in% append(snodata$id, snodata$parent_id))
corr <- corr.test(t(interest), adjust = "fdr", method = 'pearson', ci=F)

corr$p[lower.tri(corr$p,diag=TRUE)]=NA
pval <- as.data.frame(as.table(corr$p), stringsAsFactors=F)
pval <- na.omit(pval)
padj <- as.data.frame(as.table(corr$p.adj), stringsAsFactors=F)
padj <- na.omit(padj)
corr$r[lower.tri(corr$r,diag=TRUE)]=NA
Correlation <- as.data.frame(as.table(corr$r), stringsAsFactors=F)
Correlation <- na.omit(Correlation)

cortable <- (cbind( Correlation, padj, pval))
cortable <- cortable[,c(1,2,3,5,7)]

colnames(cortable) <- c("gene1","gene2","cor","p.adj", "p.val")
cortablefilt <-data.frame(cortable,stringsAsFactors=F) #[abs(cortable[,3])> 0 & cortable[,4] < 0.05 ,]
padj <- cortablefilt[,4]
padj[padj==0] <- as.numeric(unlist(format(.Machine)))[1]
cortablefilt <- cbind(cortablefilt, -log10(padj))

parentcorr <- snodata
parentcorr$ids <- paste(parentcorr$id, ";", parentcorr$parent_id)
cortablefilt$ids <- paste(cortablefilt$gene2, ";", cortablefilt$gene1)
parentcorr$corr <- cortablefilt$cor[match(parentcorr$ids, cortablefilt$ids)]
parentcorr$pad <- cortablefilt$`-log10(padj)`[match(parentcorr$ids, cortablefilt$ids)]
cortablefilt$ids <- paste(cortablefilt$gene1, ";", cortablefilt$gene2)
parentcorr$corr2 <- cortablefilt$cor[match(parentcorr$ids, cortablefilt$ids)]
parentcorr$pad2 <- cortablefilt$`-log10(padj)`[match(parentcorr$ids, cortablefilt$ids)]
parentcorr$pearson <- coalesce(parentcorr$corr2, parentcorr$corr)
parentcorr$padj <- coalesce(parentcorr$pad2, parentcorr$pad)
parentcorr[c("corr","corr2","pad","pad2")] <- NULL
parentcorr <- na.omit(parentcorr)

p <- data.frame()
p <- subset(parentcorr,padj > 1.30103)
p <- p[order(p$pearson),] 
min <- head(p,2)
max <- tail(p, 2)
p <- rbind(min,max)

ggplot(parentcorr, aes(x=pearson, y=padj, colour = padj > (-log10(0.001))))+ geom_point() + 
  scale_colour_manual(name = 'p-value < 0.05', values = setNames(c('grey','#33A02C'),c(F, T)))+
  theme_minimal()+
  xlab("Correlation of abundance between snoRNAs and their parent gene (Pearson's r)")+
  ylab("-log10(FDR-adjusted p-value)")

###################
parentcorr$group <- round(parentcorr$pearson, digits = 0)

parent_type <- c( rep('non-coding', 7), rep('protein-coding', 7)) 
pearson_cat <- rep(c('<-0.5', '-0.5;-0.25', '-0.25;0', '0;0.25', '0.25;0.5', '0.5;0.75', '0.75;1'), 2)

value<- c(nrow(parentcorr[parentcorr$has_parent == 'non-coding' & parentcorr$pearson < -0.5,]),
          nrow(parentcorr[parentcorr$has_parent == 'non-coding' & -0.5 <= parentcorr$pearson & parentcorr$pearson < -0.25,]),
          nrow(parentcorr[parentcorr$has_parent == 'non-coding' & -0.25 <= parentcorr$pearson & parentcorr$pearson < 0,]),
          nrow(parentcorr[parentcorr$has_parent == 'non-coding' & 0 <= parentcorr$pearson & parentcorr$pearson < 0.25,]),
          nrow(parentcorr[parentcorr$has_parent == 'non-coding' & 0.25 <= parentcorr$pearson & parentcorr$pearson < 0.5,]),
          nrow(parentcorr[parentcorr$has_parent == 'non-coding' & 0.5 <= parentcorr$pearson & parentcorr$pearson < 0.75,]),
          nrow(parentcorr[parentcorr$has_parent == 'non-coding' & 0.75 <= parentcorr$pearson,]),
          nrow(parentcorr[parentcorr$has_parent == 'protein-coding' & parentcorr$pearson < -0.5,]),
          nrow(parentcorr[parentcorr$has_parent == 'protein-coding' & -0.5 <= parentcorr$pearson & parentcorr$pearson < -0.25,]),
          nrow(parentcorr[parentcorr$has_parent == 'protein-coding' & -0.25 <= parentcorr$pearson & parentcorr$pearson < 0,]),
          nrow(parentcorr[parentcorr$has_parent == 'protein-coding' & 0 <= parentcorr$pearson & parentcorr$pearson < 0.25,]),
          nrow(parentcorr[parentcorr$has_parent == 'protein-coding' & 0.25 <= parentcorr$pearson & parentcorr$pearson < 0.5,]),
          nrow(parentcorr[parentcorr$has_parent == 'protein-coding' & 0.5 <= parentcorr$pearson & parentcorr$pearson < 0.75,]),
          nrow(parentcorr[parentcorr$has_parent == 'protein-coding' & parentcorr$pearson > 0.75,]))
bardata <-data.frame(parent_type,pearson_cat,value)
bardata$perc <- c(0,0,0, 1/40, 12/117, 17/54, 15/24, 1, 1, 1, 39/40, 105/117, 37/54, 9/24)*100 
bardata$pearson_cat <- as.character(bardata$pearson_cat )
bardata$pearson_cat <- factor(bardata$pearson_cat, levels=unique(bardata$pearson_cat))
ggplot(bardata, aes(x=pearson_cat, y=perc, fill=parent_type)) + 
  geom_bar(stat="identity", width=0.5, position='fill')+
  #xlim(as.character(bardata$pearson_cat))+
  theme_minimal()+
  xlab("Correlation of abundance between snoRNAs and their parent gene (Pearson's r)")+
  ylab("Proportion of snoRNAs")#+ 
  #scale_fill_brewer(palette="Paired")
bardata
###################################################GO

non <- unique(parentcorr$parent_id[which(parentcorr$has_parent == 'non-coding')])
#add function based on ZFLNC database
GO <- c("other", "proteasome assembly", "ribosome biogenesis and translation", "poorly characterized",  "ribosome biogenesis and translation", "poorly uncharacterized",
        "other", "RNA binding, processing, splicing", "RNA binding, processing, splicing", "ribosome biogenesis and translation", "poorly uncharacterized", 
        "RNA binding, processing, splicing", "ribosome biogenesis and translation", "poorly characterized")
#add function based on biomart GO names
prot <- unique(parentcorr$parent_id[which(parentcorr$has_parent == 'protein-coding')])
ribo_prot <- c("ENSDARG00000007320",
              "ENSDARG00000010516",
              "ENSDARG00000015862",
              "ENSDARG00000019181",
              "ENSDARG00000019230",
              "ENSDARG00000036298",
              "ENSDARG00000036316",
              "ENSDARG00000036875",
              "ENSDARG00000041182",
              "ENSDARG00000044093",
              "ENSDARG00000053058",
              "ENSDARG00000054818",
              "ENSDARG00000103007",
              "ENSDARG00000116649")

rna_smg <- c("ENSDARG00000002710",
             "ENSDARG00000003599",
             "ENSDARG00000006200",
             "ENSDARG00000008292",
             "ENSDARG00000010137",
             "ENSDARG00000014244",
             "ENSDARG00000016443",
             "ENSDARG00000016484",
             "ENSDARG00000019230",
             "ENSDARG00000021140",
             "ENSDARG00000022430",
             "ENSDARG00000029248",
             "ENSDARG00000032175",
             "ENSDARG00000040184",
             "ENSDARG00000041182",
             "ENSDARG00000043873",
             "ENSDARG00000057167",
             "ENSDARG00000063050",
             "ENSDARG00000075113",
             "ENSDARG00000092115",
             "ENSDARG00000099256",
             "ENSDARG00000099430",
             "ENSDARG00000103163",
             "ENSDARG00000104609"
)

ribogen <- c( "ENSDARG00000006316",
              "ENSDARG00000006691",
              "ENSDARG00000012820",
              "ENSDARG00000025073",
              "ENSDARG00000025581",
              "ENSDARG00000034291",
              "ENSDARG00000042065",
              "ENSDARG00000042094",
              "ENSDARG00000042566",
              "ENSDARG00000043304",
              "ENSDARG00000046157",
              "ENSDARG00000051783",
              "ENSDARG00000053457",
              "ENSDARG00000055868",
              "ENSDARG00000055996",
              "ENSDARG00000099380",
              "ENSDARG00000100371",
              "ENSDARG00000100392",
              "ENSDARG00000104353"
)

other <- c("ENSDARG00000005058",
           "ENSDARG00000005513",
           "ENSDARG00000005867",
           "ENSDARG00000006963",
           "ENSDARG00000012040",
           "ENSDARG00000013522",
           "ENSDARG00000014817",
           "ENSDARG00000016477",
           "ENSDARG00000017891",
           "ENSDARG00000018698",
           "ENSDARG00000020007",
           "ENSDARG00000020711",
           "ENSDARG00000025566",
           "ENSDARG00000026842",
           "ENSDARG00000030700",
           "ENSDARG00000033916",
           "ENSDARG00000034457",
           "ENSDARG00000035751",
           "ENSDARG00000036500",
           "ENSDARG00000036864",
           "ENSDARG00000039007",
           "ENSDARG00000040245",
           "ENSDARG00000040666",
           "ENSDARG00000044521",
           "ENSDARG00000056318",
           "ENSDARG00000062058",
           "ENSDARG00000063295",
           "ENSDARG00000068992",
           "ENSDARG00000078069",
           "ENSDARG00000078473",
           "ENSDARG00000079148",
           "ENSDARG00000086927",
           "ENSDARG00000090286",
           "ENSDARG00000091187",
           "ENSDARG00000099184",
           "ENSDARG00000100623",
           "ENSDARG00000101652",
           "ENSDARG00000102624",
           "ENSDARG00000103553",
           "ENSDARG00000104889",
           "ENSDARG00000106560",
           "ENSDARG00000109513",
           "ENSDARG00000116881"
)

GO2 <- c()

for (pr in prot){
  if (pr %in% ribogen){
    GO2 <- c(GO2, "ribosome biogenesis and translation")
  }
  else if (pr %in% ribo_prot){
    GO2 <- c(GO2, "ribosomal protein")
  }
  else if (pr %in% rna_smg){
    GO2 <- c(GO2, "RNA binding, processing, splicing")
  }
  else if (pr %in% other){
    GO2 <- c(GO2, "other")
  }
  else {
    GO2 <- c(GO2, "poorly characterized")
    #print(pr)
  }
}

parentGO <- c()
parentGO$id <- c(non,prot) 
parentGO$GO <- c(GO,GO2)
parentGO <- as.data.frame(parentGO)

table1$val2 <- table2$val2[match(table1$pid, table2$pid)]
parentcorr$GO <- parentGO$GO[match(parentcorr$parent_id, parentGO$id)]


###########################################################

pear <- c(rep("negative", 5), rep("positive", 5))
biological_function <- rep(c("ribosome biogenesis and translation",
           "ribosomal protein",
           "RNA binding, processing, splicing",
           "other",
           "poorly characterized"),2)
perc <- c(round(nrow(parentcorr[which(parentcorr$pearson < -0.25 & parentcorr$GO == "ribosome biogenesis and translation"),])/nrow(parentcorr[which(parentcorr$pearson < -0.25),])*100, 1),
          round(nrow(parentcorr[which(parentcorr$pearson < -0.25 & parentcorr$GO == "ribosomal protein"),])/nrow(parentcorr[which(parentcorr$pearson < -0.25),])*100, 1),
          round(nrow(parentcorr[which(parentcorr$pearson < -0.25 & parentcorr$GO == "RNA binding, processing, splicing"),])/nrow(parentcorr[which(parentcorr$pearson < -0.25),])*100, 1),
          round(nrow(parentcorr[which(parentcorr$pearson < -0.25 & parentcorr$GO == "other"),])/nrow(parentcorr[which(parentcorr$pearson < -0.25),])*100, 1),
          round(nrow(parentcorr[which(parentcorr$pearson < -0.25 & parentcorr$GO == "poorly characterized"),])/nrow(parentcorr[which(parentcorr$pearson < -0.25),])*100, 1),
          round(nrow(parentcorr[which(parentcorr$pearson > 0.25 & parentcorr$GO == "ribosome biogenesis and translation"),])/nrow(parentcorr[which(parentcorr$pearson > 0.25),])*100, 1),
          round(nrow(parentcorr[which(parentcorr$pearson > 0.25 & parentcorr$GO == "ribosomal protein"),])/nrow(parentcorr[which(parentcorr$pearson > 0.25),])*100, 1),
          round(nrow(parentcorr[which(parentcorr$pearson > 0.25 & parentcorr$GO == "RNA binding, processing, splicing"),])/nrow(parentcorr[which(parentcorr$pearson > 0.25),])*100, 1),
          round(nrow(parentcorr[which(parentcorr$pearson > 0.25 & parentcorr$GO == "other"),])/nrow(parentcorr[which(parentcorr$pearson > 0.25),])*100, 1),
          round(nrow(parentcorr[which(parentcorr$pearson > 0.25 & parentcorr$GO == "poorly characterized"),])/nrow(parentcorr[which(parentcorr$pearson > 0.25),])*100, 1))

perc <- c(round(nrow(parentcorr[which(parentcorr$pearson < -0 & parentcorr$GO == "ribosome biogenesis and translation"),])/nrow(parentcorr[which(parentcorr$pearson < -0),])*100, 1),
          round(nrow(parentcorr[which(parentcorr$pearson < -0 & parentcorr$GO == "ribosomal protein"),])/nrow(parentcorr[which(parentcorr$pearson < -0),])*100, 1),
          round(nrow(parentcorr[which(parentcorr$pearson < -0 & parentcorr$GO == "RNA binding, processing, splicing"),])/nrow(parentcorr[which(parentcorr$pearson < -0),])*100, 1),
          round(nrow(parentcorr[which(parentcorr$pearson < -0 & parentcorr$GO == "other"),])/nrow(parentcorr[which(parentcorr$pearson < -0),])*100, 1),
          round(nrow(parentcorr[which(parentcorr$pearson < -0 & parentcorr$GO == "poorly characterized"),])/nrow(parentcorr[which(parentcorr$pearson < -0),])*100, 1),
          round(nrow(parentcorr[which(parentcorr$pearson > 0 & parentcorr$GO == "ribosome biogenesis and translation"),])/nrow(parentcorr[which(parentcorr$pearson > 0),])*100, 1),
          round(nrow(parentcorr[which(parentcorr$pearson > 0 & parentcorr$GO == "ribosomal protein"),])/nrow(parentcorr[which(parentcorr$pearson > 0),])*100, 1),
          round(nrow(parentcorr[which(parentcorr$pearson > 0 & parentcorr$GO == "RNA binding, processing, splicing"),])/nrow(parentcorr[which(parentcorr$pearson > 0),])*100, 1),
          round(nrow(parentcorr[which(parentcorr$pearson > 0 & parentcorr$GO == "other"),])/nrow(parentcorr[which(parentcorr$pearson > 0),])*100, 1),
          round(nrow(parentcorr[which(parentcorr$pearson > 0 & parentcorr$GO == "poorly characterized"),])/nrow(parentcorr[which(parentcorr$pearson > 0),])*100, 1))

barplotGO <- data.frame(pear,biological_function, perc)
barplotGO
ggplot(barplotGO, aes(x=pear, y=perc, fill=biological_function)) + 
  geom_bar(position = "fill",stat="identity", width=0.4 )+
  theme_minimal()+
  xlab("Correlation of abundance between snoRNAs and their parent gene (Pearson's r)")+
  ylab("Proportion of snoRNAs")+
  scale_y_continuous(labels = scales::percent)+ 
  scale_fill_brewer(palette="Paired")

dat <- data.frame(
  "pos" = c(nrow(parentcorr[which(parentcorr$pearson > 0 & parentcorr$GO == "ribosome biogenesis and translation"),]),
            nrow(parentcorr[which(parentcorr$pearson > 0 & parentcorr$GO == "ribosomal protein"),]),
            nrow(parentcorr[which(parentcorr$pearson > 0 & parentcorr$GO == "RNA binding, processing, splicing"),]),
            nrow(parentcorr[which(parentcorr$pearson > 0 & parentcorr$GO == "other"),]),
            nrow(parentcorr[which(parentcorr$pearson > 0 & parentcorr$GO == "poorly characterized"),])),
  "neg" = c(nrow(parentcorr[which(parentcorr$pearson < -0 & parentcorr$GO == "ribosome biogenesis and translation"),]),
            nrow(parentcorr[which(parentcorr$pearson < -0 & parentcorr$GO == "ribosomal protein"),]),
            nrow(parentcorr[which(parentcorr$pearson < -0 & parentcorr$GO == "RNA binding, processing, splicing"),]),
            nrow(parentcorr[which(parentcorr$pearson < -0 & parentcorr$GO == "other"),]),
            nrow(parentcorr[which(parentcorr$pearson < -0 & parentcorr$GO == "poorly characterized"),])),
  row.names = c("ribosome biogenesis and translation",
                "ribosomal protein",
                "RNA binding, processing, splicing",
                "other",
                "poorly characterized"),
  stringsAsFactors = FALSE
)

dat
row_wise_fisher_test(dat, detailed = FALSE)
