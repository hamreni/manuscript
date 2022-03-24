library(ggplot2)
library(DESeq2)
library(tidyverse)
library(dplyr)
library(ggrepel)
library(RColorBrewer)
library(rtracklayer)

coldata <- read.csv(file = "data/OURcoldata.csv", sep=";")
countdata <- read.csv(file = "data/OURexpression.csv", sep=";", row.names = 1)

gtf <- import("data/expanded_GRCz11-104-annot-sorted.gtf") %>% as_tibble %>% distinct(gene_id, gene_name, gene_biotype)

idx <- match(row.names(countdata), gtf$gene_id )
countdata$type <-  gtf$gene_biotype [ idx ]
countdata <- countdata[1:ncol(countdata)]

#set selected characteristic (based on coldata)
select <- which(coldata$tissue %in% "brain")
countdata <- countdata[,select]

total <- countdata %>% 
  group_by(type) %>% 
  summarise_all(sum)
total <- total %>% mutate(total = rowSums(.[2:ncol(countdata)]))

pie <- 1
pie$type <- total$type
pie$values <- total$total
pie <- as.data.frame(pie[2:3])
pie<- pie[order(pie$values),] 
pie <- subset(pie, type!= "rRNA")
pie <- pie%>%group_by(type)%>%mutate(Percentage=as.numeric(paste0(round(values/sum(pie$values)*100, 2))))
sno <- subset(pie, type== "snoRNA")
pie <- subset(pie, type!= "snoRNA")

threshold <- 10
up <- pie[which(pie$Percentage > threshold),]
up <- rbind(up,sno)
low <- subset(pie, Percentage < threshold)
low$type <- "other"
low <- aggregate(Percentage~type, low, FUN=sum)
low$values <- 0
dat <- rbind(up,low)

dat<- dat[order(dat$type),] 
dat$type <- as.factor(dat$type)
dat <- dat %>%
  arrange(desc(type)) %>%
  mutate(lab.ypos = cumsum(Percentage) - 0.5*Percentage)

#set label position
dat$pos <- c(0.64, 1.4, 46, 95, 98)

brewer.pal(8, "Pastel2")
mycols <- c("lincRNA"="#A6CEE3", 
            "miRNA"= "#1F78B4", 
            "ncRNA"="#B2DF8A", 
            "snoRNA"="#33A02C", 
            "other"="#FB9A99", 
            "pre_miRNA"="#E31A1C", 
            "protein_coding"="#FDBF6F", 
            "tRNA"="#FF7F00",
            "lncRNA" = "#CAB2D6",
            "snRNA" = "#6A3D9A",
            "processed_transcript" = "#FFFF99",
            "antisense" = "#B3E2CD",
            "piRNA" = "#FDCDAC",
            "pseudogene" ="#CBD5E8",
            "unprocessed_pseudogene" = "#F4CAE4")

ggplot(dat, aes(x = "", y = Percentage, fill = type)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0)+
  geom_label_repel(data = dat,
                   aes(y = pos, label = paste0(Percentage, "%")),
                   size = 4.5, nudge_x = 1, show.legend = F) +
  scale_fill_manual(values = mycols) +
  theme_void()

