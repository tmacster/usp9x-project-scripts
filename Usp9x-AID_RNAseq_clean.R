library(gdata)
library(gplots)
library(ggrepel)
library(DESeq2)
library(GGally)
library(pheatmap)
library(edgeR)
library(tidyverse)
library(viridis)
library(reshape2)

setwd("/Users/tmac/Documents/Usp9x-project/RNAseq/")

##if doing CNN/spike-in norm: for other types of graphs etc, use code from the ERCC sections.

######################################################################################################
############################            voom version of ERCC norm                     ################  
#######################################################################################################

#NB expressed is data removing any zero values. # is.expressed = apply(seqdata, 1, function(row) all(row !=0 ))

####################   READ IN RAW COUNTS  ####################
seqdata <- read.table("seq_readcounts.txt", header=T, row.names=1, quote="") #use if previously exported data table.
seqdata <- seqdata[,6:23] #remove first 5 columns of stuff
#ncol(seqdata) #n = 18 
seqdata <- seqdata[,-8] # remove the sample with super low read count (aid2-8h_high)
colnames(seqdata) <- c("neg1", "neg2", "neg3","low1_8h","low2_8h","low3_8h","high1_8h","high3_8h","flag1","flag2","flag3",
                     "low1_48h","low2_48h","low3_48h","high1_48h","high2_48h","high3_48h")

#head(seqdata)

##############################################################
#get a table from subset of seqdata, only for ercc
ercc <- seqdata[grep("^ERCC", rownames(seqdata)),] # confirm that there are 92 spike-ins

# filter data to remove non-expressed genes/ERCCs 
is.expressed <- rowSums(cpm(seqdata) > 0) >= 3 #only keep genes where >3 samples have a cpm >0. needs edgeR
is.expressed <- apply(seqdata, 1, function(row) all(row !=0 )) #remove any rows where there's a zero value 
raw_expr2 <- seqdata[is.expressed,] # n=15499 
raw_expr <- raw_expr2[,c(1:3,9:11,4:8,12:17)] # reorder to move flag samples next to neg

# get a table from subset of seqdata, for ercc. WROTE TO FILE
ercc_expr <- raw_expr[grep("^ERCC", rownames(raw_expr)),] 
#nrow(ercc_expr) # 59 ERCCs expressed in at least 3 samples
#write.table(ercc_expr, "ERCC_expressed_new.txt", sep="\t", row.names=T, quote=F)
genes_expr <- raw_expr[!rownames(raw_expr) %in% rownames(ercc_expr),] #remove the ERCC data from the table
#write.table(raw_expr, "RAW_expressed_new.txt",sep="\t", quote=F, row.names=T)

#examine ERCC linearity between samples
#pairs(~neg1+neg2+neg3+low1_8h+low2_8h+low3_8h+high1_8h+high3_8h+flag1+flag2+flag3+low1_48h+low2_48h+
 #       low3_48h+high1_48h+low2_48h+low3_48h+high1_48h+high2_48h+high3_48h, data=ercc)
#pairs(ercc[,1:9]) #8h samples

################################
## READ IN for ercc analysis ###
ercc_expr <- read.table("ERCC_expressed_new.txt", header=T, row.names=1, quote="")
genes_expr <- read.table("RAW_expressed_new.txt", header=T, row.names=1, quote="")
#head(ercc_expr)

####

SampleName <- c("neg1", "neg2", "neg3", "flag1","flag2","flag3","low1_8h", "low2_8h", "low3_8h", "high1_8h", "high3_8h", 
                "low1_48h","low2_48h","low3_48h","high1_48h","high2_48h","high3_48h")
Condition <- c("ctrl","ctrl","ctrl","ctrl","ctrl","ctrl","low_8h","low_8h","low_8h","high_8h","high_8h","low_48h","low_48h",
               "low_48h","high_48h","high_48h","high_48h")
info <- cbind(SampleName,Condition)

condition <- factor(info[, "Condition"]) #Compare data by sample type
design <- model.matrix(~0+condition) #creates the matrix that tells you which condition it is. Try it - 1 = yes, 0 = no.
colnames(design) <- levels(condition)

#normalize the data using ERCC SPIKE-INs - see https://support.bioconductor.org/p/74870/
N <- colSums(genes_expr)
nf <- calcNormFactors(ercc_expr, lib.size=N)
voom.norm <- voom(genes_expr, design, lib.size = N * nf, plot=F)

#head(voom.norm$E) #- this is expression counts. Norm count + 0.5 log2. Check it looks ok on plot.

voom_genes <- voom.norm$E

#write.table(voom.norm$E, "voom/VOOM_NORM-expr.genes_log2exprs_flagnoaux_final.txt", sep="\t", quote=F, row.names=T)   #export ERCC norm count data 

#READ IN DATA ########################################################################################
voom_genes <-read.table("voom/VOOM_NORM-expr.genes_log2exprs_flagnoaux_final.txt", header=T, row.names=1, quote="")        
######################################################################################################


############################### Density plots of raw vs ERCC-norm data ##############################
# first call the environment of expressed genes
raw_av <- data.frame(neg=rowMeans(genes_expr[,1:3]), flag=rowMeans(genes_expr[,4:6]), low_8h=rowMeans(genes_expr[,7:9]), high_8h=rowMeans(genes_expr[,10:11]), 
                     low_48h=rowMeans(genes_expr[,12:14]), high_48h=rowMeans(genes_expr[,15:17]))
head(raw_av)
tail(raw_av)

#################  density plot with voom genes.   ##############
voom_av <- data.frame(neg=rowMeans(voom_genes[,1:3]), flag=rowMeans(voom_genes[,4:6]), low_8h=rowMeans(voom_genes[,7:9]), high_8h=rowMeans(voom_genes[,10:11]),
                      low_48h=rowMeans(voom_genes[,12:14]), high_48h=rowMeans(voom_genes[,15:17]))

head(voom_av)

dat_8h <- data.frame(dens = c(voom_av$neg, voom_av$low_8h, voom_av$high_8h), key=rep(c("no auxin", "8h low", "8h high"), each = nrow(raw_av)))
ggplot(dat_8h, aes(x = dens, fill = key)) + 
  geom_density(alpha = 0.3)
dat_48h <- data.frame(dens = c(voom_av$flag, voom_av$low_48h, voom_av$high_48h), key=rep(c("flag", "48h low", "48h high"), each = nrow(raw_av)))
ggplot(dat_48h, aes(x = dens, fill = key)) + 
  geom_density(alpha = 0.3)

################# VOOM PCA - check clustering of data ##########################

pcaV <- prcomp(t(voom_genes), scale.=T, center=T)
scores <- as.data.frame(pcaV$x)

eigs <- pcaV$sdev^2
percentVar <- round(100 * eigs/sum(eigs), 1)


colors = factor(c("grey85","grey85","grey85", "grey35","grey35","grey35", "lightblue","lightblue","lightblue","salmon1","salmon1",
                "#134153","#134153","#134153","#902e2e","#902e2e","#902e2e")) 

ggplot(data=scores, aes(x=PC1, y=PC2, label=rownames(scores))) +
  geom_point(color=colors, size=3) +
  geom_hline(yintercept = 0, color = "black", size=0.5) +
  geom_vline(xintercept = 0, color = "black", size=0.5) +
  #geom_text_repel(colour = "black", alpha = 0.8, size = 4) +
  theme(axis.text.x=element_text(size=12, color = "black"), axis.text.y=element_text(size=12, color = "black")) +
  theme(axis.title.x=element_text(vjust=-0.35, size=16), axis.title.y=element_text(vjust=0.35, size=16)) +
  theme(panel.background=element_rect(fill="transparent", color="black", size=1, linetype="solid"), plot.background = element_rect(fill = "transparent"), panel.grid.major=element_blank(), panel.grid.minor = element_blank()) +
  labs(x = paste0("PC1 (",percentVar[1],"% of variance)"), y = paste0("PC2 (",percentVar[2],"% of variance)"), title = "PCA") +
  ylim(c(-100,110)) + xlim(c(-150,215)) #+
  #ggsave("RNAseq_voom_PCA.pdf", device = "pdf", dpi = 300, width = 4, height = 4, units = "in")


#################################################################################
################################## VOOM CORRELATION HEATMAP #####################
#################################################################################

#use voom_genes (see read in above)
genes_corr <- cor(voom_genes, method = "spearman")  
genes_corr
par(oma = c(3,4,3,4))
heatmap.2(genes_corr, scale="none", col = viridis, key=TRUE, symkey=T, density.info="none", trace="none")
pheatmap(genes_corr, color=viridis(20), border_color = NA)

#################################################################################
####### TOPTABLE ANALYSIS _ DE TESTING on voom-norm ##################
#################################################################################

fit <- lmFit(voom.norm, design) #fits a linear model for each gene. Takes the voom.norm object.

# ---------------- 48h_low vs ctrl (flag/no auxin) --------------------- #
contrasts48h_low <- makeContrasts(low_48h - ctrl, levels = design) #compares up to three things, here is 48h low vs controls
contr48h_low.fit <- eBayes(contrasts.fit(fit, contrasts48h_low)) #fits to the Bayesian model, comparing the chosen groups
diff48h_low <- topTable(contr48h_low.fit, coef=NULL, number=15499) #performs topTable analysis of genes for the defined genes. number is nrow(genes_expr)

diff48h_low_export <- diff48h_low[order(-diff48h_low$logFC),] #sort by descending logFC
diff48h_low_export.padj05 <- subset(diff48h_low_export, adj.P.Val < 0.05) #changed to adj P <0.05 #6008 genes
diff48h_low_export.padj05.logFC0.7 <- subset(diff48h_low_export.padj05, logFC > 0.7 | logFC < -0.7) #3410
#write.table(diff48h_low_export.padj05.logFC0.7, "voom/TopTable-48h_LOWvflagnoaux-padj05,logFC0.7.txt", sep="\t", quote=F, row.names=T)
#write.table(diff48h_low_export, "voom/TopTable-48h_LOWvflagnoaux-ALL.txt", sep="\t", quote=F, row.names=T)
#write.table(diff48h_low_export, "../RNAseq_voom/TopTable-48h_LOWvflagnoaux-ALL_15k.txt", sep="\t", quote=F, row.names=T)

nrow(subset(diff48h_low_export.padj05, logFC < -0.7)) #979
nrow(subset(diff48h_low_export.padj05, logFC > 0.7)) #2434

diff48h_low_export_GSEA <- diff48h_low_export
row.names(diff48h_low_export_GSEA) <- toupper(row.names(diff48h_low_export_GSEA))
head(as.matrix(diff48h_low_export_GSEA)[,1])
#write.table((as.matrix(diff48h_low_export_GSEA)[,1]), "voom/GSEA/voomGSEA_48h_lowvflagnoauxin.RNK",sep="\t", quote=F, row.names=T, col.names = F)

# ---------------- 48h_high vs ctrl (flag/no auxin) --------------------- #
contrasts48h_high <- makeContrasts(high_48h - ctrl, levels = design)
contr48h_high.fit <- eBayes(contrasts.fit(fit, contrasts48h_high))
diff48h_high <- topTable(contr48h_high.fit, coef=NULL, number=15499)

diff48h_high_export <- diff48h_high[order(-diff48h_high$logFC),] 
diff48h_high_export.padj05 <- subset(diff48h_high_export, adj.P.Val < 0.05) #11226 genes
diff48h_high_export.padj05.logFC0.7 <- subset(diff48h_high_export.padj05, logFC > 0.7 | logFC < -0.7) #3979
#write.table(diff48h_high_export.padj05.logFC0.7, "voom/TopTable-48h_HIGHvflagnoaux-padj05,logFC0.7.txt", sep="\t", quote=F, row.names=T)
#write.table(diff48h_high_export, "voom/TopTable-48h_HIGHvflagnoaux-ALL.txt", sep="\t", quote=F, row.names=T)
#write.table(diff48h_high_export, "../RNAseq_voom/TopTable-48h_HIGHvflagnoaux-ALL_15k.txt", sep="\t", quote=F, row.names=T)

nrow(subset(diff48h_high_export.padj05, logFC < -0.7)) #3941
nrow(subset(diff48h_high_export.padj05, logFC > 0.7)) #38

diff48h_high_export_GSEA <- diff48h_high_export
row.names(diff48h_high_export_GSEA) <- toupper(row.names(diff48h_high_export_GSEA))
head(as.matrix(diff48h_high_export_GSEA)[,1])
#write.table((as.matrix(diff48h_high_export_GSEA)[,1]), "voom/GSEA/voomGSEA_48h_highvflagnoauxin.RNK",sep="\t", quote=F, row.names=T, col.names = F)

################################# 8h_low vs ctrl #################################
contrasts8h_low <- makeContrasts(low_8h - ctrl, levels = design)
contr8h_low.fit <- eBayes(contrasts.fit(fit, contrasts8h_low))
diff8h_low <- topTable(contr8h_low.fit, coef=NULL, number=15499)

diff8h_low_export <- diff8h_low[order(-diff8h_low$logFC),] 
diff8h_low_export.padj05 <- subset(diff8h_low_export, adj.P.Val < 0.05) #4002 genes
diff8h_low_export.padj05.logFC0.7 <- subset(diff8h_low_export.padj05, logFC > 0.7 | logFC < -0.7) #103 at log(2), 2155 at log(0.7)
#write.table(diff8h_low.padj05.logFC0.7, "RNAseq_voom/TopTable-8h_LOWvflagnoaux-padj05,logFC0.7.txt", sep="\t", quote=F, row.names=T)
#write.table(diff8h_low_export.padj05, "TopTable-8h_LOW-padj05.txt", sep="\t", quote=F, row.names=T)
#write.table(diff8h_low_export, "RNAseq_voom/TopTable-8h_LOWvflagnoaux-ALL_15k.txt", sep="\t", quote=F, row.names=T)

nrow(subset(diff8h_low_export.padj05, logFC < -0.7)) #845
nrow(subset(diff8h_low_export.padj05, logFC > 0.7)) #1312

diff8h_low_export_GSEA <- diff8h_low_export
row.names(diff8h_low_export_GSEA) <- toupper(row.names(diff8h_low_export_GSEA))
head(as.matrix(diff8h_low_export_GSEA)[,1])
#write.table((as.matrix(diff8h_low_export_GSEA)[,1]), "voom/GSEA/voomGSEA_8h_lowvflagnoauxin_15k.RNK",sep="\t", quote=F, row.names=T, col.names = F)

## compare the two 8h samples --> what gets most hypertranscribed?
contrasts_8h <- makeContrasts(low_8h - high_8h, levels = design)
contr_8h.fit <- eBayes(contrasts.fit(fit, contrasts_8h))
diff_8h <- topTable(contr_8h.fit, coef=NULL, number=15499)

diff_8h_sorted <- diff_8h[order(-diff_8h$logFC),] #sort by descending logFC #12974 genes

################################# 8h_high vs ctrl #################################
contrasts8h_high <- makeContrasts(high_8h - ctrl, levels = design) 
contr8h_high.fit <- eBayes(contrasts.fit(fit, contrasts8h_high)) 
diff8h_high <- topTable(contr8h_high.fit, coef=NULL, number=15499) 

diff8h_high_export <- diff8h_high[order(-diff8h_high$logFC),] #sort by descending logFC 
diff8h_high_export.padj05 = subset(diff8h_high_export, adj.P.Val < 0.05) #changed to adj P <0.05 #423 genes
diff8h_high_export.padj05.logFC0.7 = subset(diff8h_high_export.padj05, logFC > 0.7 | logFC < -0.7) #347
#write.table(diff8h_high_export.padj05.logFC0.7, "RNAseq_voom/TopTable-8h_HIGHvflagnoaux-padj05,logFC0.7.txt", sep="\t", quote=F, row.names=T)
#write.table(diff8h_high_export.padj05, "TopTable-8h_HIGH-padj05.txt", sep="\t", quote=F, row.names=T)
#write.table(diff8h_high_export, "RNAseq_voom/TopTable-8h_HIGHvflagnoaux-ALL_15k.txt", sep="\t", quote=F, row.names=T)

nrow(subset(diff8h_high_export.padj05, logFC < -0.7)) #277
nrow(subset(diff8h_high_export.padj05, logFC > 0.7)) #70

diff8h_high_export_GSEA <- diff8h_high_export
row.names(diff8h_high_export_GSEA) <- toupper(row.names(diff8h_high_export_GSEA))
head(as.matrix(diff8h_high_export_GSEA)[,1])
#write.table((as.matrix(diff8h_high_export_GSEA)[,1]), "voom/GSEA/voomGSEA_8h_highvflagnoauxin_15k.RNK",sep="\t", quote=F, row.names=T, col.names = F)

# ---------------- flag vs no-auxin --------------------- #
# only 63 genes are significantly de between the flag cells & no-auxin control cells
SampleName <- c("neg1", "neg2", "neg3", "flag1","flag2","flag3","low1_8h", "low2_8h", "low3_8h", "high1_8h", "high3_8h", 
                "low1_48h","low2_48h","low3_48h","high1_48h","high2_48h","high3_48h")
Condition <- c("ctrl","ctrl","ctrl","flag","flag","flag","low_8h","low_8h","low_8h","high_8h","high_8h","low_48h","low_48h",
               "low_48h","high_48h","high_48h","high_48h")
info <- cbind(SampleName,Condition)
condition <- factor(info[, "Condition"]) #Compare data by sample type
design <- model.matrix(~0+condition) #creates the matrix that tells you which condition it is. Try it - 1 = yes, 0 = no.
colnames(design) <- levels(condition)

fit_flag <- lmFit(voom.norm, design) #fits a linear model for each gene. Takes the voom.norm object.

contrasts_flag <- makeContrasts(flag - ctrl, levels = design) 
contrasts_flag.fit <- eBayes(contrasts.fit(fit_flag, contrasts_flag)) #fits to the Bayesian model, comparing the chosen groups
flag_diff <- topTable(contrasts_flag.fit, coef=NULL, number=15499)
flag_diff_tbl <- as_tibble(flag_diff, rownames = "GeneID") %>% rename(padj = adj.P.Val) %>% filter(padj < 0.05 & logFC > 0.7)

##################################################################################
################################# HEATMAPS  ######################################
##################################################################################

#heatmaps of average data
my_palette = viridis()
my_palette2 <- colorRampPalette(c("dodgerblue3", "khaki2", "darkred"))(n=299)

heatmap.2(as.matrix(voom_genes), dendrogram="column", col=my_palette2, scale="row", key=T, density.info="none", 
          trace="none", cexCol=1, Colv=T, labRow=F, distfun=pcorr, main="voom row means expression")

#### collect ribosomal gene data
ribo3 <-grep("Rpl|Rps", rownames(voom_genes)) #get row numbers which contain Rpl or Rps in the title.
voom_ribo <- voom_genes[ribo3,]

voom_ribo_av <- filter(voom_av, str_detect(GeneID, "Rps|Rpl")) %>% filter(str_detect(GeneID, "Rps6k", negate = T)) %>% 
  transmute(GeneID, ctrl = log2(ctrl + 0.5), low_48h = log2(low_48h + 0.5), high_48h = log2(high_48h + 0.5)) %>% na.omit()


ribo_av <- grep("Rpl|Rps", rownames(voom_av)) #get row numbers which contain Rpl or Rps in the title.
voom_ribo_av <- voom_av[ribo_av,]
r6k_av <- voom_ribo_av[grep("Rps6k", rownames(voom_ribo_av)),]
voom_ribo_av2 <- voom_ribo_av[!rownames(voom_ribo_av) %in% rownames(r6k_av),] #remove weird Rps6k
#nrow(voom_ribo_av2) #90


# heatmaps of ribosomal genes
pcorr <- function(x) as.dist(1-cor(t(x),method="pearson"))

# all samples
heatmap.2(as.matrix(voom_ribo[,c(1:6,12:17)]), dendrogram="column", col=viridis, scale="row", key=T, density.info="none", 
          trace="none", cexCol=1, Colv=T, labRow=rownames(voom_ribo), distfun=pcorr, main = "Rps/Rpl voom")

# averages (only plots some of ribo genes....)
#heatmap.2(as.matrix(voom_ribo_av2[,c(-2,-5,-6)]), dendrogram = "column", col=viridis, scale="row", key=T, density.info="none", 
 #         trace="none", cexCol=1, Colv=T, labRow=rownames(voom_ribo_av2), distfun=pcorr, main = "Rps/Rpl row means voom")

#heatmap.2(as.matrix(voom_ribo_av[,-1]), dendrogram = "column", col=viridis, scale="row", key=T, density.info="none", 
 #         trace="none", cexCol=1, Colv=T, labRow=voom_ribo_av$GeneID, distfun=pcorr, main = "Rps/Rpl row means voom")

#pheatmap that essentially plots z-score
pheatmap(voom_ribo_av[-1], color = viridis(20), labels_row = voom_ribo_av$GeneID, scale = "row", clustering_distance_rows = "correlation", border_color = NA, 
         fontsize_row = 5, cluster_cols = TRUE, treeheight_row = 0, treeheight_col = 0,
         filename = "../RNAseq_voom/printed_plots/ribo_heatmap_48h_voom.tiff", width = 2, height = 6)

# now do the same for rlog norm data
ribo3de <- grep("Rpl|Rps", rownames(exprsR))
de_ribo <- exprsR[ribo3de,]
#heatmap.2(as.matrix(de_ribo), dendrogram="column", col=viridis, scale="row", key=T, density.info="none", 
 #         trace="none", cexCol=1, Colv=T, labRow=T, distfun=pcorr, main = "Rps/Rpl rlog de")



##################################################################################
# ----- heatmaps to compare DEGs between samples
##################################################################################
# read in toptable outupts
#diff48h_low_export <-read_tsv("voom/TopTable-48h_LOWvflagnoaux-ALL.txt", col_names = c("GeneID", "logFC", "AveExpr", "t", "pval", "padj", "B")) %>% slice(-1) %>% select(GeneID, logFC)        
#diff48h_high_export <-read_tsv("voom/TopTable-48h_HIGHvflagnoaux-ALL.txt", col_names = c("GeneID", "logFC", "AveExpr", "t", "pval", "padj", "B")) %>% slice(-1) %>% select(GeneID, logFC)        
#diff8h_high_export <-read_tsv("voom/TopTable-8h_HIGHvflagnoaux-ALL.txt", col_names = c("GeneID", "logFC", "AveExpr", "t", "pval", "padj", "B")) %>% slice(-1) %>% select(GeneID, logFC)        
#diff8h_low_export <-read_tsv("voom/TopTable-8h_LOWvflagnoaux-ALL.txt", col_names = c("GeneID", "logFC", "AveExpr", "t", "pval", "padj", "B")) %>% slice(-1) %>% select(GeneID, logFC)        

# select the DEG at 8h low
diff8h_low_export.padj05.logFC0.7 <- read_tsv("voom/TopTable-8h_LOWvflagnoaux-padj05,logFC0.7.txt", col_names = c("GeneID", "logFC", "AveExpr", "t", "pval", "padj", "B")) %>% slice(-1)
up <- filter(diff8h_low_export.padj05.logFC0.7, logFC > 0)
down <- filter(diff8h_low_export.padj05.logFC0.7, logFC < 0)

voom_genes_8hlowDE_up <- subset(voom_genes, rownames(voom_genes) %in% up$GeneID)
voom_genes_8hlowDE_down <- subset(voom_genes, rownames(voom_genes) %in% down$GeneID)

# now plot
pcorr <- function(x) as.dist(1-cor(t(x),method="pearson"))
my_palette = viridis(30)
heatmap.2(as.matrix(voom_genes_8hlowDE_up), dendrogram="column", col=viridis, scale="row", key=T, density.info="none", 
          trace="none", cexCol=1, Colv=T, labRow=T, distfun=pcorr, main = "genes UP 8h low")

heatmap.2(as.matrix(voom_genes_8hlowDE_down), dendrogram="column", col=viridis, scale="row", key=T, density.info="none", 
          trace="none", cexCol=1, Colv=T, labRow=T, distfun=pcorr, main = "genes DOWN 8h low")

#### create voom_av object
# pool together controls in a voom_av object
voom_av <- data.frame(ctrl=rowMeans(voom_genes[,1:6]), low_8h=rowMeans(voom_genes[,7:9]), high_8h=rowMeans(voom_genes[,10:11]),
                      low_48h=rowMeans(voom_genes[,12:14]), high_48h=rowMeans(voom_genes[,15:17])) %>% as_tibble(rownames = "GeneID")

#plot averages 
voom_genes_av_8hlowDE_up <- subset(voom_av, rownames(voom_genes) %in% up$GeneID)
voom_genes_av_8hlowDE_down <- subset(voom_av, rownames(voom_genes) %in% down$GeneID)

heatmap.2(as.matrix(voom_genes_av_8hlowDE_up), dendrogram="column", col=viridis, scale="row", key=T, density.info="none", 
          trace="none", cexCol=1, Colv=T, labRow=T, distfun=pcorr, main = "genes UP 8h low")
heatmap.2(as.matrix(voom_genes_av_8hlowDE_down), dendrogram="column", col=viridis, scale="row", key=T, density.info="none", 
          trace="none", cexCol=1, Colv=T, labRow=T, distfun=pcorr, main = "genes DOWN 8h low")


##################################################################################
# comparing 8h low to 48h low
de_48h <- as_tibble(diff48h_low, rownames = "GeneID")  %>% dplyr::rename(log2FC = logFC, pvalue = P.Value, padj = adj.P.Val) %>% filter(padj < 0.05)
de_8h <- as_tibble(diff8h_low, rownames = "GeneID") %>% dplyr::rename(log2FC = logFC, pvalue = P.Value, padj = adj.P.Val) %>% filter(padj < 0.05)

de_usp9x_low <- inner_join(de_48h_2, de_8h_2, by = "GeneID", suffix = c("_48h", "_8h")) %>% arrange(-log2FC_48h)

# heatmap of genes de at 48h
my_palette2 <- colorRampPalette(c("dodgerblue3", "white", "khaki2", "khaki3"))(n = 50)
pheatmap(de_usp9x_low[,-1], cluster_cols = FALSE, cluster_rows = FALSE, color = my_palette2) # see https://stackoverflow.com/questions/31677923/set-0-point-for-pheatmap-in-r for breaks

de_usp9x_low %>% melt %>% ggplot() + 
  geom_tile(aes(variable, value, fill = value)) +
  scale_fill_gradient2(low = "cadetblue3", mid = "white", high = "goldenrod2", midpoint = 0) +
  theme(panel.background = element_rect(fill = "transparent", color = "black", size = 1), plot.background = element_rect(fill = "transparent"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme_bw()

ggplot(toptable_low8_long, aes(x = variable, y = reorder(GeneID, value), fill = value)) + 
  geom_tile(width = 0.4) + 
  scale_fill_gradient2(low = "cadetblue3", mid = "white", high = "goldenrod2", midpoint = 0) +
  geom_text(aes(label = GeneID, color = category), nudge_x = 0.3, size = 2.5, hjust = 0, fontface = ifelse(abs(toptable_low8_long$value) > 0.4, "bold", "plain")) +
  scale_color_manual(values = c("forestgreen", "blue")) + 
  theme_void() #theme(axis.text.y = element_text(color = ifelse(toptable_low_long$GeneID %in% primed, "blue", "forestgreen")))

# scatter plot to compare
de_usp9x_low_padj <- inner_join(de_48h, de_8h, by = "GeneID", suffix = c("_48h", "_8h")) %>% select(GeneID, matches("log2FC")) %>% arrange(-log2FC_48h) %>%
  mutate(cluster = case_when(log2FC_8h > 0 & log2FC_48h > 0 ~ "up_in_both", log2FC_8h < 0 & log2FC_48h < 0 ~ "down_in_both", log2FC_8h > 0 & log2FC_48h < 0 ~ "up_in_8h", log2FC_8h < 0 & log2FC_48h > 0 ~ "up_in_48h"))

de_usp9x_low_padj %>% ggplot(aes(log2FC_8h, log2FC_48h, color = cluster)) + 
  geom_point(size = 0.7) + 
  #scale_color_viridis(discrete = TRUE) +
  scale_color_manual(values = c("#2E0854", "grey", "grey", "goldenrod2")) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  #geom_hline(yintercept = 0.7, linetype = "dotted") + geom_vline(xintercept = 0.7, linetype = "dotted") + geom_hline(yintercept = -0.7, linetype = "dotted") + geom_vline(xintercept = -0.7, linetype = "dotted") +
  #xlim(c(-3.5,3.5)) + ylim(c(-4,6)) +
  theme(axis.text.x = element_text(size = 12, color = "black"), axis.text.y = element_text(size = 12, color = "black")) +
  theme(axis.title.x = element_text(vjust = -0.35, size = 12), axis.title.y = element_text(vjust = 0.35, size = 12)) +
  theme(panel.background = element_rect(fill = "transparent", color = "black", size = 1), plot.background = element_rect(fill = "transparent"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #+
  #ggsave("correlation_8h_48h_low_toptable.pdf", device = "pdf", dpi = 300, width = 4.5, height = 3, units = "in")

de_usp9x_low_padj %>% group_by(cluster) %>% tally
#1 down_in_both   852
#2 up_in_48h       90
#3 up_in_8h        46
#4 up_in_both    1809

de_usp9x_low_padj %>% group_by(cluster) %>% select(-GeneID) %>% summarise_each(mean)
#  cluster      log2FC_48h log2FC_8h
#1 down_in_both     -1.01     -0.969
#2 up_in_48h         0.773    -0.726
#3 up_in_8h         -0.841     0.782
#4 up_in_both        1.35      0.885

cor(de_usp9x_low_padj[,c(2:3)], method = "pearson") # 0.83 between samples

filter(de_usp9x_low_padj, cluster == "up_in_8h") %>% View()

########

##################################################################################
#- hypergeometric tests for Venn diagrams
##################################################################################

a <- 5 #overlap
b <- 66-a #DE, not in ChIP target data
c <- 596-a #ChIP target, not DE
d <- 15499-a-b-c #not target and not DE (15499 as universe of expressed genes)

suz12_venn <- matrix(data = c(a, b, c, d), nrow = 2)
fisher.test(suz12_venn) #p-value < 2.2e-16


ezh2_venn <- matrix(data = c(a, b, c, d), nrow = 2)
fisher.test(ezh2_venn)


#############
# plot naive vs pluripotency gene expr
naive <- c("Fgf4",	"Dppa3", "Tfcp2l1", "Nr0b1", "Dppa5", "Zfp42",	"Tbx3",	"Nr5a2", "Prdm14", "Klf2", "Klf5",	"Klf4",	"Nanog",	"Tet2",	"Tead4", "Spp1", "Stat3", "Tcf3",
           "Prdm16", "Sp100", "Dazl", "Trpm1", "Crxos", "Bmp4")
#formative <- c("Pou5f1", "Sox2", "Sall4")
primed <- c("Otx2", "Lef1", "Dnmt3b", "Dnmt3a", "Foxd3", "Sox3", "Fgf5", "Pou3f1", "Sox1", "Nes", "Pax6", "Meis2", "Zfp281", "Zic3", "Lefty2", "Hoxb1", "Zic1", "Zic3", 
            "Dlx3", "Otx2", "Foxa2", "Sox4", "Lef1", "Lin28a") #Mixl1
#primed <- c("Sox6",	"Hoxb3",	"Pou3f2",	"Rfx4",	"Dll1",	"Meis1",	"Lmo2",	"Tet3",	"Sox1",	"Zic1",	"Meis2",	"Znf521",	"Tet1",	"Dusp6",	"Zic3",
#           "Dnmt3a",	"Dnmt3b",	"Zic2", "Zic5", "Cebpb",	"Otx2", "Fgf5", "Foxp1", "Sall2", "Zyx", "Lef1", "Zfp281")

toptable_low8 <- low_8h %>% filter(GeneID %in% naive | GeneID %in% primed) %>% select(GeneID, log2FC) %>% add_column(category = ifelse(.$GeneID %in% primed, "primed", "naive"))
toptable_high8 <- high_8h %>% filter(GeneID %in% naive | GeneID %in% primed) %>% select(GeneID, log2FC)
toptable_low48 <- low_48h %>% filter(GeneID %in% naive | GeneID %in% primed) %>% select(GeneID, log2FC)


toptable_8h <- inner_join(toptable_low8, toptable_high8, by = "GeneID", suffix = c("_low", "_high")) %>% arrange(desc(log2FC_high)) 

# re-plot as a ggplot with text customization
# make GeneID a factor based on ranking in 8h sample & annotate by color: [[#TEXT LABELS DON'T WORK]]
toptable_8h_long <- toptable_8h %>% 
  add_column(color = case_when(.$category == "primed" ~ "blue", TRUE ~ "forestgreen")) %>%
  mutate(GeneID = fct_reorder(.$GeneID, log2FC_high), color = fct_reorder(.$color, log2FC_low)) %>% 
  melt(value.name = "log2FC") 


ggplot(toptable_8h_long, aes(x = variable, y = GeneID, fill = log2FC)) + 
  geom_text(aes(label = GeneID, color = category), hjust = 1, size = 2.5) +
  scale_color_manual(values = c("forestgreen", "blue")) +
  theme_classic() + 
  theme(axis.text.y = element_blank(), 
        axis.text.x = element_text(color = "black", angle = 45, vjust = 0.5), 
        axis.line.x = element_blank(), axis.line.y = element_blank(), 
        panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #+
#ggsave("pluripotency_genes_Usp9x_8h_annotations_NEW.pdf", device = "pdf", dpi = 300, width = 3, height = 6, units = "in")

# add in 48h time points
pluripotency_expr <- list(toptable_low8, toptable_high8, toptable_low48) %>% 
  reduce(inner_join, by = "GeneID", suffix = c("_low8", "_high8")) %>%
  rename(log2FC_low48 = log2FC) %>% 
  arrange(desc(log2FC_high8)) %>%
  mutate(GeneID = fct_reorder(.$GeneID, log2FC_low8)) 

pheatmap(pluripotency_expr[,c(2,4,5)], labels_row = pluripotency_expr$GeneID, color = viridis(20), cluster_rows = FALSE, fontsize = 4, annotation_colors = ann_colors, 
         treeheight_col = 25, filename = "~/Documents/Usp9x-project/RNAseq_voom/printed_plots/pluripotency_genes_Usp9x_8h_heatmap_NEW.pdf", width = 1.5)


