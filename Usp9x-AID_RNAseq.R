# code to analyze CNN RNA-seq with ERCC spike-ins 
# Trisha Macrae, PhD with help from M. Percharde (https://github.com/mpercharde/RNAseq)

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

#setwd("/dir")


####################   READ IN RAW COUNTS  ####################
# read counts available from GEO (GSE146800)

seqdata <- read_tsv("EScell_readcounts.txt", skip = 1, 
                    col_names = c("GeneID", "chr", "start", "end", "length", 
                                  "neg1", "neg2", "neg3","low1_8h","low2_8h","low3_8h","high1_8h", "high2_8h", "flag1","flag2","flag3",
                                  "low1_48h","low2_48h","low3_48h","high1_48h","high2_48h","high3_48h")) %>%
  select(-chr:-length) %>%
  mutate_at(vars(neg1:high3_48h), as.numeric) %>%
  filter(GeneID != "1-Mar" & GeneID != "2-Mar") %>% #gene name error
  column_to_rownames(var = "GeneID")

##############################################################
ercc <- seqdata[grep("^ERCC", rownames(seqdata)),] # confirm that there are 92 spike-ins

# filter data to remove non-expressed genes/ERCCs 
is.expressed <- rowSums(cpm(seqdata) > 0) >= 3 #only keep genes where >3 samples have a cpm >0. needs edgeR
is.expressed <- apply(seqdata, 1, function(row) all(row !=0 )) #remove any rows where there's a zero value 
raw_expr2 <- seqdata[is.expressed,] # n=15499 
raw_expr <- raw_expr2[,c(1:3,9:11,4:8,12:17)] # reorder to move flag samples next to neg

# get ercc expression
ercc_expr <- raw_expr[grep("^ERCC", rownames(raw_expr)),] 
#nrow(ercc_expr) # 59 ERCCs expressed in at least 3 samples
genes_expr <- raw_expr[!rownames(raw_expr) %in% rownames(ercc_expr),] #remove the ERCC data from the table

#write.table(ercc_expr, "ERCC_expressed_new.txt", sep="\t", row.names=T, quote=F)
#write.table(raw_expr, "RAW_expressed_new.txt",sep="\t", quote=F, row.names=T)

#examine ERCC linearity between samples
#pairs(~neg1+neg2+neg3+low1_8h+low2_8h+low3_8h+high1_8h+high3_8h+flag1+flag2+flag3+low1_48h+low2_48h+
 #       low3_48h+high1_48h+low2_48h+low3_48h+high1_48h+high2_48h+high3_48h, data=ercc)
#pairs(ercc[,1:9]) #8h samples

###################################
## read in raw data for analysis ##
#ercc_expr <- read.table("ERCC_expressed_new.txt", header=T, row.names=1, quote="")
#genes_expr <- read.table("RAW_expressed_new.txt", header=T, row.names=1, quote="")
#head(ercc_expr)

####
SampleName <- c("neg1", "neg2", "neg3", "flag1","flag2","flag3","low1_8h", "low2_8h", "low3_8h", "high1_8h", "high3_8h", 
                "low1_48h","low2_48h","low3_48h","high1_48h","high2_48h","high3_48h")
Condition <- c("ctrl","ctrl","ctrl","ctrl","ctrl","ctrl","low_8h","low_8h","low_8h","high_8h","high_8h","low_48h","low_48h",
               "low_48h","high_48h","high_48h","high_48h")
info <- cbind(SampleName, Condition)

condition <- factor(info[, "Condition"]) #Compare data by sample type
design <- model.matrix(~0 + condition) #creates the matrix that tells you which condition it is. Try it - 1 = yes, 0 = no.
colnames(design) <- levels(condition)

#normalize the data using ERCC spike-ins - see https://support.bioconductor.org/p/74870/
N <- colSums(genes_expr)
nf <- calcNormFactors(ercc_expr, lib.size=N)
voom.norm <- voom(genes_expr, design, lib.size = N * nf, plot=F)

voom_genes <- voom.norm$E #- expression counts as Norm count + 0.5 log2, can check that it looks ok

#write.table(voom.norm$E, "voom-norm_log2exprs.txt", sep="\t", quote=F, row.names=T)   #export ERCC norm count data 


############################### 
# quality checks
############################### 

# ----- Density plots of raw vs ERCC-norm data
# first call the environment of expressed genes
raw_av <- data.frame(neg=rowMeans(genes_expr[,1:3]), flag=rowMeans(genes_expr[,4:6]), low_8h=rowMeans(genes_expr[,7:9]), high_8h=rowMeans(genes_expr[,10:11]), 
                     low_48h=rowMeans(genes_expr[,12:14]), high_48h=rowMeans(genes_expr[,15:17]))
#head(raw_av)

voom_av <- data.frame(neg=rowMeans(voom_genes[,1:3]), flag=rowMeans(voom_genes[,4:6]), low_8h=rowMeans(voom_genes[,7:9]), high_8h=rowMeans(voom_genes[,10:11]),
                      low_48h=rowMeans(voom_genes[,12:14]), high_48h=rowMeans(voom_genes[,15:17]))

# density plots
#dat_8h <- data.frame(dens = c(voom_av$neg, voom_av$low_8h, voom_av$high_8h), key=rep(c("no auxin", "8h low", "8h high"), each = nrow(raw_av)))
#ggplot(dat_8h, aes(x = dens, fill = key)) + 
 # geom_density(alpha = 0.3)
#dat_48h <- data.frame(dens = c(voom_av$flag, voom_av$low_48h, voom_av$high_48h), key=rep(c("flag", "48h low", "48h high"), each = nrow(raw_av)))
#ggplot(dat_48h, aes(x = dens, fill = key)) + 
 # geom_density(alpha = 0.3)

###
# ---- PCA

pca <- prcomp(t(voom_genes), scale.=T, center=T)
scores <- as.data.frame(pca$x)
eigs <- pca$sdev^2
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

###
# ---- pairwise sample correlation 

genes_corr <- cor(voom_genes, method = "spearman") 
par(oma = c(3,4,3,4))
pheatmap(genes_corr, color=viridis(20), border_color = NA)

########################################################
####### TOPTABLE ANALYSIS: DE TESTING ##################

fit <- lmFit(voom.norm, design) # fits a linear model for each gene

# ---------------- 48h_low vs ctrl (flag/no auxin) --------------------- #
contrasts48h_low <- makeContrasts(low_48h - ctrl, levels = design) #compares up to three things, here is 48h low vs controls
contr48h_low.fit <- eBayes(contrasts.fit(fit, contrasts48h_low)) #fits to the Bayesian model, comparing the chosen groups
diff48h_low <- topTable(contr48h_low.fit, coef=NULL, number=15499) #performs topTable analysis of genes for the defined genes. number is nrow(genes_expr)

diff48h_low_export <- diff48h_low[order(-diff48h_low$logFC),] #sort by descending logFC
diff48h_low_export.padj05 <- subset(diff48h_low_export, adj.P.Val < 0.05) #changed to adj P <0.05 #6008 genes
diff48h_low_export.padj05.logFC0.7 <- subset(diff48h_low_export.padj05, logFC > 0.7 | logFC < -0.7) #3410
#write_tsv(diff48h_low_export, "xx.txt")

nrow(subset(diff48h_low_export.padj05, logFC < -0.7)) #979
nrow(subset(diff48h_low_export.padj05, logFC > 0.7)) #2434


# ---------------- 48h_high vs ctrl (flag/no auxin) --------------------- #
contrasts48h_high <- makeContrasts(high_48h - ctrl, levels = design)
contr48h_high.fit <- eBayes(contrasts.fit(fit, contrasts48h_high))
diff48h_high <- topTable(contr48h_high.fit, coef=NULL, number=15499)

diff48h_high_export <- diff48h_high[order(-diff48h_high$logFC),] 
diff48h_high_export.padj05 <- subset(diff48h_high_export, adj.P.Val < 0.05) #11226 genes
diff48h_high_export.padj05.logFC0.7 <- subset(diff48h_high_export.padj05, logFC > 0.7 | logFC < -0.7) #3979

nrow(subset(diff48h_high_export.padj05, logFC < -0.7)) #3941
nrow(subset(diff48h_high_export.padj05, logFC > 0.7)) #38


# ----------------8h_low vs ctrl--------------------- #
contrasts8h_low <- makeContrasts(low_8h - ctrl, levels = design)
contr8h_low.fit <- eBayes(contrasts.fit(fit, contrasts8h_low))
diff8h_low <- topTable(contr8h_low.fit, coef=NULL, number=15499)

diff8h_low_export <- diff8h_low[order(-diff8h_low$logFC),] 
diff8h_low_export.padj05 <- subset(diff8h_low_export, adj.P.Val < 0.05) #4002 genes
diff8h_low_export.padj05.logFC0.7 <- subset(diff8h_low_export.padj05, logFC > 0.7 | logFC < -0.7) #103 at log(2), 2155 at log(0.7)
#write_tsv(diff8h_low_export, "xx.txt")

nrow(subset(diff8h_low_export.padj05, logFC < -0.7)) #845
nrow(subset(diff8h_low_export.padj05, logFC > 0.7)) #1312


# ----------------8h_high vs ctrl--------------------- #
contrasts8h_high <- makeContrasts(high_8h - ctrl, levels = design) 
contr8h_high.fit <- eBayes(contrasts.fit(fit, contrasts8h_high)) 
diff8h_high <- topTable(contr8h_high.fit, coef=NULL, number=15499) 

diff8h_high_export <- diff8h_high[order(-diff8h_high$logFC),] #sort by descending logFC 
diff8h_high_export.padj05 = subset(diff8h_high_export, adj.P.Val < 0.05) #changed to adj P <0.05 #423 genes
diff8h_high_export.padj05.logFC0.7 = subset(diff8h_high_export.padj05, logFC > 0.7 | logFC < -0.7) #347
#write_tsv(diff8h_high_export, "xx.txt")

nrow(subset(diff8h_high_export.padj05, logFC < -0.7)) #277
nrow(subset(diff8h_high_export.padj05, logFC > 0.7)) #70

#####
# create the ranked gene file for GSEA: 
diff48h_low_export_GSEA <- diff48h_low_export
row.names(diff48h_low_export_GSEA) <- toupper(row.names(diff48h_low_export_GSEA))
#head(as.matrix(diff48h_low_export_GSEA)[,1])
#write.table((as.matrix(diff48h_low_export_GSEA)[,1]), "voom/GSEA/voomGSEA_48h_lowvflagnoauxin.RNK",sep="\t", quote=F, row.names=T, col.names = F)

################################# 
# GENE HEATMAPS of ercc-norm data

pcorr <- function(x) as.dist(1 - cor(t(x), method = "pearson"))

# heatmaps of average data (all genes)
heatmap.2(as.matrix(voom_genes), dendrogram="column", col=viridis, scale="row", key=T, density.info="none", 
          trace="none", cexCol=1, Colv=T, labRow=F, distfun=pcorr, main="voom row means expression")

# heatmap of ribosomal gene data (Rpl/Rps) - Supp Fig 2b
voom_av_tbl <- as_tibble(voom_av, rownames = "GeneID") %>%
  mutate(ctrl = (neg + flag)/2) %>%
  select(-neg, -flag)

voom_ribo_av <- filter(voom_av_tbl, str_detect(GeneID, "Rps|Rpl")) %>% 
  filter(str_detect(GeneID, "Rps6k", negate = T)) %>% 
  transmute(GeneID, ctrl = log2(ctrl + 0.5), low_48h = log2(low_48h + 0.5), high_48h = log2(high_48h + 0.5)) %>% 
  na.omit()

# pheatmap that essentially plots z-score
pheatmap(voom_ribo_av[-1], color = viridis(20), labels_row = voom_ribo_av$GeneID, scale = "row", 
         clustering_distance_rows = "correlation", border_color = NA, fontsize_row = 5, cluster_cols = TRUE, 
         filename = "ribo_heatmap_48h.tiff",
         treeheight_row = 0, treeheight_col = 0, width = 2, height = 6)

# individual replicates
#heatmap.2(as.matrix(voom_ribo[,c(1:6,12:17)]), dendrogram="column", col=viridis, scale="row", key=T, density.info="none", 
#         trace="none", cexCol=1, Colv=T, labRow=rownames(voom_ribo), distfun=pcorr, main = "Rps/Rpl voom")

#########################################################
#----------------- hypergeometric tests for Venn diagrams
#########################################################

a <- 5 #overlap
b <- 66-a #DE, not in ChIP target data
c <- 596-a #ChIP target, not DE
d <- 15499-a-b-c #not target and not DE (15499 as universe of expressed genes)

suz12_venn <- matrix(data = c(a, b, c, d), nrow = 2)
fisher.test(suz12_venn) #p-value < 2.2e-16

ezh2_venn <- matrix(data = c(a, b, c, d), nrow = 2)
fisher.test(ezh2_venn)

#########################################################
#---------- naive vs pluripotency gene expr (Supp Fig. 1f)

naive <- c("Fgf4",	"Dppa3", "Tfcp2l1", "Nr0b1", "Dppa5", "Zfp42",	"Tbx3",	"Nr5a2", "Prdm14", "Klf2", "Klf5",	"Klf4",	"Nanog",	"Tet2",	
           "Tead4", "Spp1", "Stat3", "Tcf3", "Prdm16", "Sp100", "Dazl", "Trpm1", "Crxos", "Bmp4")
primed <- c("Otx2", "Lef1", "Dnmt3b", "Dnmt3a", "Foxd3", "Sox3", "Fgf5", "Pou3f1", "Sox1", "Nes", "Pax6", "Meis2", "Zfp281", "Zic3", 
            "Lefty2", "Hoxb1", "Zic1", "Zic3", "Dlx3", "Otx2", "Foxa2", "Sox4", "Lef1", "Lin28a")

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

# for text...
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
  purrr::reduce(inner_join, by = "GeneID", suffix = c("_low8", "_high8")) %>%
  dplyr::rename(log2FC_low48 = log2FC) %>% 
  arrange(desc(log2FC_high8)) %>%
  mutate(GeneID = fct_reorder(.$GeneID, log2FC_low8)) 

pheatmap(pluripotency_expr[,c(2,4,5)], labels_row = pluripotency_expr$GeneID, color = viridis(20), cluster_rows = FALSE, fontsize = 4, annotation_colors = ann_colors, 
         treeheight_col = 25, filename = "pluripotency_genes_Usp9x_8h_heatmap.pdf", width = 1.5)

#################################################
## ------ boxplots/violin plots of expression

# expression of all genes at 48h (Fig. 1g)
toptable_long <- toptable %>% 
  pivot_longer(cols = starts_with("log2FC"), names_to = "stage", names_prefix = "log2FC_", values_to = "log2FC") %>%
  mutate(factor = case_when(stage == "high8" ~ 1, stage == "high48" ~ 2, stage == "low8" ~ 3, TRUE ~ 4))

toptable_long %>% filter(grepl("48", stage)) %>% 
  ggplot(aes(fct_reorder(stage, factor), log2FC, fill = stage)) + 
  geom_hline(yintercept = 0) +
  geom_violin(color = "black", scale = "width") + 
  #stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1), geom = "pointrange", color = "black") +
  geom_boxplot(fill = "white", width = 0.2, outlier.shape = NA, color = "black") + 
  scale_fill_manual(values = c("#902E2E", "#436c7f")) +
  #ylim(c(-2.5, 5)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), axis.text = element_text(color = "black", size = "14")) + #, axis.text.x = element_text(angle = 90))
  theme(panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #+ 
#ggsave("RNAseq_log2fc_48h.pdf", device = "pdf", dpi = 300, width = 3, height = 3, units = "in")

wilcox.test(toptable$log2FC_high48, toptable$log2FC_low48)

# Suz12 targets based on Enrichr (Fig. 1i)
suz12 <- read_tsv("Suz12_targets_Pasini.txt", col_names = "GeneID") %>%
  unique()

suz12_toptable <- toptable %>% filter(GeneID %in% suz12$GeneID)

suz12_toptable_long <- suz12_toptable %>% 
  pivot_longer(cols = starts_with("log2FC"), names_to = "stage", names_prefix = "log2FC_", values_to = "log2FC") %>%
  mutate(factor = case_when(stage == "high8" ~ 1, stage == "high48" ~ 2, stage == "low8" ~ 3, TRUE ~ 4))

# get a random gene subset
set.seed(0)

random_long <- toptable %>% 
  sample_n(nrow(suz12_toptable), replace = FALSE) %>%
  pivot_longer(cols = starts_with("log2FC"), names_to = "stage", names_prefix = "log2FC_", values_to = "log2FC") %>%
  mutate(factor = case_when(stage == "high8" ~ 1, stage == "high48" ~ 2, stage == "low8" ~ 3, TRUE ~ 4))

random_suz12 <- bind_rows("suz12" = suz12_toptable_long, "random" = random_long, .id = "subset") 

random_suz12 %>% 
  ggplot(aes(fct_reorder(stage, factor), log2FC, fill = subset)) + 
  geom_hline(yintercept = 0) +
  stat_boxplot(geom = "errorbar", position = position_dodge(width = 0.7), width = 0.3) +
  geom_boxplot(color = "black", position = position_dodge(width = 0.7), width = 0.6, outlier.shape = NA) + #geom_boxplot(fill = "white", width = 0.1, outlier.shape = NA, color = "black") + 
  scale_fill_manual(values = c("grey", "#E793B7")) +
  ylim(c(-2.5, 3.5)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), axis.text = element_text(color = "black", size = "12"), axis.text.x = element_text(angle = 90), 
        panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #+ 
#ggsave("RNAseq_suz12target_expr.pdf", device = "pdf", dpi = 300, width = 4, height = 3, units = "in")

# calculate median for each
random_suz12 %>% group_by(subset, stage) %>%
  summarize(median = median(log2FC), sd = sd(log2FC))

# stats: ANOVA with pair-wise comparison
random_suz12_stats <- random_suz12 %>% unite(col = "group", c(subset, stage)) %>%
  select(-GeneID, -factor)res.aov <- aov(log2FC ~ stage + subset, data = random_suz12)

summary(res.aov) # p <2e-16
pairwise.wilcox.test(random_suz12_stats$log2FC, random_suz12_stats$group, p.adj = "fdr")


# individual tests...
random_suz12_wide <- random_suz12 %>% select(-factor) %>%
  pivot_wider(names_from = c(subset, stage), values_from = log2FC)
wilcox.test(random_suz12_wide$suz12_low8, random_suz12_wide$random_low8) #  p-value < 2.2e-16

