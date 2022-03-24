library(ggplot2)
library(DESeq2)
library(tidyverse)
library(rtracklayer)
library(FinCal)
library(cubature)
library(miscTools)


coldata <- read.csv(file = "data/OURcoldata.csv", sep=";")
countdata <- read.csv(file = "data/OURexpression.csv", sep=";", row.names = 1)


gtf <- import("data/expanded_GRCz11-104-annot-sorted.gtf") %>% as_tibble %>% distinct(gene_id, gene_name, gene_biotype)
snodata <- read.csv("data/sno_genes.csv", sep=";")
idx <- match(row.names(countdata), gtf$gene_id )
countdata$type <-  gtf$gene_biotype [ idx ]
countdata <- countdata[which(countdata$type=='snoRNA'),]
countdata <- countdata[1:length(countdata)-1]

dds <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~tissue)
keep <- rowSums(counts(dds)) > 0
dds <- dds[keep,]
dds <- DESeq(dds)

c <- (counts(dds, normalized=T))
normdata <- data.frame(counts(dds, normalized=T))
normdata$avg <- (rowMeans(normdata[1:length(countdata)]))
normdata$sd <- rowSds(c)
normdata$cv <- coefficient.variation(sd = normdata$sd, avg= normdata$avg)*100

d <- density(normdata$cv)
plot(d,  main="CV distribution of all expressed snoRNAs during development", 
     xlab = "Coefficient of variance", frame = FALSE, col = "darkred")


derivatives <- (ddnorm(d$y))
x0 <- which.min(derivatives)
spl <- smooth.spline(d$y ~ d$x)
lines(spl, col=2)
newx <- x0
pred0 <- predict(spl, x=newx, deriv=0)
pred1 <- predict(spl, x=newx, deriv=1)
yint <- pred0$y - (pred1$y*newx)
xint <- -yint/pred1$y
rug(jitter(normdata$cv))
points(pred0, col=2, pch=19)
lines(d$x, yint + pred1$y*d$x, col=3)

res <- results(dds)
idx <- match(row.names(res), gtf$gene_id )
res$biotype <- gtf$gene_biotype [ idx ]
res <- res[which(res$biotype=='snoRNA'),]
idx2 <- match(row.names(res), snodata$id )
res$parented <- snodata$has_parent [ idx2 ]
idx <- match(row.names(res), row.names(normdata) )
res$cv <-  normdata$cv [ idx ]
 result <- as.data.frame(res)
 
results %>% 
  count(parented) %>% 
  mutate(perc = n / nrow(results)) -> res2

ggplot(res2, aes(fill=parented, y=perc, x="Expressed snoRNas")) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values=c("darkred", "#E69F00", "#999999")) + 
  scale_y_continuous(labels = scales::percent) + 
  theme_minimal()

