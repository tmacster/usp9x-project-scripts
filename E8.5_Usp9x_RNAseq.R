# code to analyze RNA-seq from Usp9x ko embryos
# Trisha Macrae, PhD

library(tidyverse)
library(gplots)
library(DESeq2)
library(pheatmap)
library(edgeR)
library(RColorBrewer)
library(ggrepel)
library(viridis)
library(reshape2)
library(readxl)
library(writexl)

#setwd("current_dir/")

###### read in raw counts ###### 
seqdata <- read_tsv("E8.5_Usp9x_ko_rawcounts.txt", comment = "#", 
                    col_names = c("GeneID", "chr", "start", "end", "strand", "length", "ctrl1", "ctrl2", "ctrl3", "ctrl4", "ctrl5", "ctrl6", "mut1", "mut2", "mut3", "mut4", "mut5", "mut6")) %>%
  dplyr::slice(-1) %>%
  select(-chr:-length)
head(seqdata)

countdata <- seqdata %>% 
  mutate_at(vars(ctrl1:mut6), as.numeric) %>%
  column_to_rownames("GeneID")

############################################################################################################
#    Normalization and heatmap using Deseq - RAW DATA NORMALIZED TO READ DEPTH. #raw_expr is the table     #
############################################################################################################

samplename <- c("ctrl1", "ctrl2", "ctrl3", "ctrl4", "ctrl5", "ctrl6", "mut1", "mut2", "mut3", "mut4", "mut5", "mut6")
condition <- c(rep(c("ctrl", "mut"), each = 6))
litter <- c(rep(c("A", "B", "A", "B"), each = 3))
info <- as.data.frame(cbind(samplename, condition))
info

matrix <- DESeqDataSetFromMatrix(countData = countdata,  colData = info, design = ~ condition)
keep <- rowSums(counts(matrix)) >= 10 # filter to remove low-count genes (must have at least 10 total reads across samples, http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#pre-filtering)
matrix <- matrix[keep,]

matrixR <- rlog(matrix)
exprsR <- assay(matrixR) # use this for looking at data, clustering, etc; RAW TABLE exprs is for DE testing
head(exprsR)
#write.table(exprsR, "Usp9x_ko_8.5-rlogNormCounts.txt",sep="\t", quote=F, row.names=T)

exprsR_tbl <- as_tibble(exprsR, rownames = "GeneID")

# heatmap and clustering for all genes
pcorr <- function(x) as.dist(1-cor(t(x), method="pearson"))                    
col.cell <- c("grey", "#00429D")[info$condition] # color vector for celltype variable

heatmap.2(as.matrix(exprsR), dendrogram="column", col=inferno(20), scale="row", key=T, density.info="none", ColSideColors=col.cell,
          trace="none", cexCol=1, Colv=T, distfun=pcorr)

# heatmap for pairwise comparisons
counts_cor <- cor(as.matrix(exprsR), method = "pearson") # mut5 clusters among litter 1 controls
pheatmap(counts_cor, color = magma(20))

# heatmap of variable genes
var_genes <- apply(exprsR, 1, var)
select_var <- names(sort(var_genes, decreasing=TRUE))[1:2000] # Get the gene names for the top most variable genes
highly_variable_lcpm <- exprsR[select_var,]
col.cell <- c("grey","#00429D")[info$condition] # color vector for celltype variable

heatmap.2(highly_variable_lcpm, col=inferno(20), trace="none", main="Top 2000 most variable genes", ColSideColors=col.cell,
          density.info="none", distfun = pcorr, scale="row")

#pdf(file="2000_var_genes_blue_heatmap.pdf", width = 6, height = 5) # Save the heatmap
heatmap.2(highly_variable_lcpm, col=viridis(20), trace="none", main="2000 most variable genes", ColSideColors=col.cell, labRow=FALSE, dendrogram = "column", 
          Colv = TRUE, density.info="none", distfun = pcorr, scale="row")
#dev.off()

# PCA plot
pca <- prcomp(t(exprsR), scale.=T, center=T)
scores <- as.data.frame(pca$x)

eigs <- pca$sdev^2
percentVar <- round(100 * eigs/sum(eigs), 1)
names <- c("ctrl1", "ctrl2", "ctrl3", "ctrl4", "ctrl5", "ctrl6", "mut1", "mut2", "mut3", "mut4", "mut5", "mut6")
colors <- factor(c(rep("grey", times = 6), rep("#00429D", times = 6)))

ggplot(data = scores, aes(x = PC1, y = PC2, label = names)) +
  geom_point(color = colors, size = 3, alpha = 0.8) +
  geom_hline(yintercept = 0, color = "black") + 
  geom_vline(xintercept = 0, color= "black") + 
  #geom_text(color = colors, size = 4, hjust = -0.1, vjust = 0) + 
  #geom_text_repel(color = colors, nudge_x = -1, direction = "both", point.padding = 0.1, fontsize = 16) + 
  theme(axis.text.x = element_text(size = 12, color = "black"), axis.text.y = element_text(size = 12, color = "black")) +
  theme(axis.title.x = element_text(vjust = -0.35, size = 14), axis.title.y = element_text(vjust = 0.35, size = 14)) +
  theme(panel.background = element_rect(fill = "transparent", color = "black", size = 1), plot.background = element_rect(fill = "transparent"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  labs(x = paste0("PC1: ",percentVar[1],"% variance"), y = paste0("PC2: ",percentVar[2],"% variance"), title = "PCA plot") +
  coord_fixed() +
  theme(plot.title = element_text(size = 20, lineheight = 0.8, vjust = 2)) #+
  #ggsave("exprsR_PCA.pdf", device = "pdf", dpi = 300)


########## detect Usp9x in each sample 
cpm <- cpm(countdata)
cpm_tbl <- as_tibble(cpm) %>% add_column(GeneID = rownames(countdata)) 

usp9x <- filter(cpm_tbl, GeneID == "Usp9x")
barplot(as.matrix(usp9x[,1:12])) # confirms ko

##################################
####### TOPTABLE ANALYSIS ########
##################################

# for all mut vs ctrls
DEmatrix <- DESeq(matrix) #RAW data for DE analysis 
toptable <- results(DEmatrix, tidy = T) %>% dplyr::rename(GeneID = row) 

toptable2 <- toptable %>% arrange(-log2FoldChange)
#write_tsv(toptable2, "toptable/toptable_Usp9x_ko-vs-ctrl.txt", col_names = TRUE, quote_escape = FALSE)

toptable_sort <- filter(toptable2, padj < 0.1)
#write_tsv(toptable_sort, "toptable/toptable_Usp9x_ko-vs-ctrl_padj0.1.txt", col_names = TRUE, quote_escape = FALSE)

#---- read in toptable data (for all mut vs all controls)
toptable_sort <- read_tsv("toptable/toptable_Usp9x_ko-vs-ctrl_padj0.1.txt", col_names = TRUE)


# create object for GSEA
#toptable_gsea <- toptable_sort %>% select(GeneID, log2FoldChange) %>% mutate(GeneID = toupper(GeneID))
#write_tsv(toptable_gsea, "GSEA.rnk", col_names = T, quote_escape = F)

######################### 
# repeat DE by litter
######################### 

# -- litter 1
countdata_1 <- countdata[,c(1:3,7:9)]

samplename_1 <- c("ctrl1", "ctrl2", "ctrl3", "mut1", "mut2", "mut3")
condition_1 <- c(rep(c("ctrl", "mut"), each = 3))
info_1 <- as.data.frame(cbind(samplename_1, condition_1))

matrix1 <- DESeqDataSetFromMatrix(countData = countdata_1, colData = info_1, design = ~ condition_1) 
keep1 <- rowSums(counts(matrix1)) >= 10 
matrix1 <- matrix1[keep1,]
matrixR1 <- rlog(matrix1)  
exprsR1 <- assay(matrixR1)

exprsR1_tbl <- as_tibble(exprsR1, rownames = "GeneID")

# PCA plot
pca1 <- prcomp(t(exprsR1), center=T) #scale.=T, 
scores1 <- as.data.frame(pca1$x)
eigs <- pca1$sdev^2
percentVar <- round(100 * eigs/sum(eigs), 1)

colors <- factor(c("grey","grey","grey", "#00429D","#00429D","#00429D")) 
ggplot(data = scores1, aes(x = PC1, y = PC2, label = samplename_1)) +
  geom_point(color = colors, size = 2) +
  geom_hline(yintercept = 0, color = "black") + geom_vline(xintercept = 0, color= "black") + 
  geom_text_repel(color = colors, size = 4) + 
  theme(axis.text.x = element_text(size = 12, color = "black"), axis.text.y = element_text(size = 12, color = "black")) +
  theme(axis.title.x = element_text(vjust = -0.35, size = 14), axis.title.y = element_text(vjust = 0.35, size = 14)) +
  labs(x = paste0("PC1 (",percentVar[1],"% of variance)"), y = paste0("PC2 (",percentVar[2],"% of variance"), title = "PCA plot") +
  theme(panel.background = element_rect(fill = "transparent", color = "black", size = 1), plot.background = element_rect(fill = "transparent"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle("PC1 / PC2 on litter 1") + 
  theme(plot.title = element_text(size = 14, lineheight = 0.8, vjust = 2)) #+
  #ggsave("exprsR_litter1_PCA.pdf", device = "pdf", dpi = 300, height = 3, width = 3, units = "in")


## DE testing
DEmatrix1 <- DESeq(matrix1) 
toptable_1 <- results(DEmatrix1, tidy = T)

toptable_1_sort <- toptable_1 %>% dplyr::rename(GeneID = row) %>% arrange(-log2FoldChange)
#toptable_1_sort %>% write_tsv("analysis/toptable/toptable_Usp9x_ko-vs-ctrl_litter1.txt", col_names = TRUE, quote_escape = FALSE)
#toptable_1_sort %>% filter(padj < 0.1) %>% write_tsv("toptable/toptable_Usp9x_ko-vs-ctrl_litter1_padj0.1.txt", col_names = TRUE, quote_escape = FALSE)

#toptable_1_sort <- read_tsv("toptable/toptable_Usp9x_ko-vs-ctrl_litter1_padj0.1.txt")

########################## 
#---- litter 2
countdata_2 <- countdata[,-c(1:3,7:9)]

samplename_2 <- c("ctrl4", "ctrl5", "ctrl6", "mut4", "mut5", "mut6")
condition_2 <- c(rep(c("ctrl", "mut"), each = 3))
info_2 <- as.data.frame(cbind(samplename_2, condition_2))

matrix2 <- DESeqDataSetFromMatrix(countData = countdata_2,  colData = info_2, design = ~ condition_2) 
keep2 <- rowSums(counts(matrix2)) >= 10 
matrix2 <- matrix2[keep2,]
matrixR2 <- rlog(matrix2)  
exprsR2 <- assay(matrixR2)

exprsR2_tbl <- as_tibble(exprsR2, rownames = "GeneID")

###### 
# PCA plot
pca2 <- prcomp(t(exprsR2), scale.=T, center=T)
scores2 <- as.data.frame(pca2$x)
eigs <- pca2$sdev^2
percentVar <- round(100 * eigs/sum(eigs), 1)

colors <- factor(c("grey","grey","grey", "#00429D","#00429D","#00429D")) 
ggplot(data = scores2, aes(x = PC1, y = PC2, label = samplename_2)) +
  geom_point(color = colors, size = 2) +
  geom_hline(yintercept = 0, color = "black") + 
  geom_vline(xintercept = 0, color= "black") + 
  geom_text_repel(color = colors, size = 4) + 
  theme(axis.text.x = element_text(size = 12, color = "black"), axis.text.y = element_text(size = 12, color = "black")) +
  theme(axis.title.x = element_text(vjust = -0.35, size = 14), axis.title.y = element_text(vjust = 0.35, size = 14)) +
  labs(x = paste0("PC1 (",percentVar[1],"% of variance)"), y = paste0("PC2 (",percentVar[2],"% of variance"), title = "PCA plot") +
  theme(panel.background = element_rect(fill = "transparent", color = "black", size = 1), plot.background = element_rect(fill = "transparent"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle("PC1 / PC2 on litter 1") + 
  theme(plot.title = element_text(size = 14, lineheight = 0.8, vjust = 2)) #+
  #ggsave("exprsR_litter2_PCA_blue.pdf", device = "pdf", dpi = 300, height = 3, width = 3, units = "in")

# ----- DE testing
DEmatrix2 <- DESeq(matrix2) 
toptable_2 <- results(DEmatrix2, tidy = T) 

toptable_2_sort <- toptable_2 %>% dplyr::rename(GeneID = row) %>% arrange(-log2FoldChange)
#toptable_2_sort %>% filter(padj < 0.1) %>% write_tsv("toptable/toptable_Usp9x_ko-vs-ctrl_litter2_padj0.1.txt", col_names = TRUE, quote_escape = FALSE)

#toptable_2_sort <- read_tsv("toptable/toptable_Usp9x_ko-vs-ctrl_litter2_padj0.1.txt")

##################
# export shared up/down gene sets to a file
toptable_combined <- inner_join(toptable_1_sort, toptable_2_sort, by = "GeneID", suffix = c("_1", "_2")) %>% select(-matches("base|stat|pvalue|lfc"))
up_in_both <- filter_at(toptable_combined, vars(contains("log2")), all_vars(. > 0)) %>% filter_at(vars(contains("padj")), all_vars(. < 0.1))
down_in_both <- filter_at(toptable_combined, vars(contains("log2")), all_vars(. < 0)) %>% filter_at(vars(contains("padj")), all_vars(. < 0.1))

de_in_both_list <- list("up_in_both" = up_in_both, "down_in_both" = down_in_both)
#write_xlsx(de_in_both_list, "de_both_litters_padj0.1.xlsx")

####
# define gene subsets that are up/down in both litters (after individual toptable analyses)

litter1 <- read_tsv("toptable_Usp9x_ko-vs-ctrl_litter1.txt") %>% mutate(direction = case_when(padj < 0.1 & log2FoldChange > 0 ~ "up", padj < 0.1 & log2FoldChange < 0 ~ "down", TRUE ~ "unchanged")) %>% 
  select(-baseMean, -lfcSE:-padj) %>% filter(direction != "unchanged")
litter2 <- read_tsv("toptable_Usp9x_ko-vs-ctrl_litter2.txt") %>% mutate(direction = case_when(padj < 0.1 & log2FoldChange > 0 ~ "up", padj < 0.1 & log2FoldChange < 0 ~ "down", TRUE ~ "unchanged")) %>% 
  select(-baseMean, -lfcSE:-padj) %>% filter(direction != "unchanged")
de_in_both_full <- inner_join(litter1, litter2, by = "GeneID", suffix = c("_1", "_2")) %>% filter(direction_1 == direction_2) %>% select(-matches("Litter")) 

# define sets for later
up_in_both <- filter(de_in_both_full, direction_1 == "up" & direction_2 == "up") %>% select(-matches("dir")) %>% dplyr::rename(litter1 = log2FoldChange_1, litter2 = log2FoldChange_2) %>% melt(variable = "litter", value.name = "log2FoldChange")
down_in_both <- filter(de_in_both_full, direction_1 == "down" & direction_2 == "down") %>% select(-matches("dir")) %>% dplyr::rename(litter1 = log2FoldChange_1, litter2 = log2FoldChange_2) %>% melt(variable = "litter", value.name = "log2FoldChange")

de_in_both <- de_in_both_full %>% select(GeneID, direction = direction_1)


####################################
# plot genes up in E8.5 mutants over 48h low ES cell data

# import ES celltoptables
low_8h <- read_tsv("TopTable-8h_low.txt", col_names = c("GeneID", "log2FC", "AveExpr", "t", "pvalue", "padj", "B"), col_types = "cnnnnnn", skip = 1)

high_8h <- read_tsv("TopTable-8h_high.txt", col_names = c("GeneID", "log2FC", "AveExpr", "t", "pvalue", "padj", "B"), col_types = "cnnnnnn", skip = 1)

low_48h <- read_tsv("TopTable-48h_low.txt", col_names = c("GeneID", "log2FC", "AveExpr", "t", "pvalue", "padj", "B"), col_types = "cnnnnnn", skip = 1)

high_48h <- read_tsv("TopTable-48h_high.txt", col_names = c("GeneID", "log2FC", "AveExpr", "t", "pvalue", "padj", "B"), col_types = "cnnnnnn", skip = 1)

early <- inner_join(low_8h, high_8h, by = "GeneID", suffix = c("_low", "_high"))
late <- inner_join(low_48h, high_48h, by = "GeneID", suffix = c("_low", "_high"))
toptable <- inner_join(early, late, by = "GeneID", suffix = c("8", "48")) 

# re-import 
embryo <- read_excel("de_both_litters_padj0.1.xlsx")

embryo_genes_toptable <- toptable %>% filter(GeneID %in% de_in_both$GeneID)

set.seed(0)
embryo_genes_random <- toptable %>% sample_n(nrow(embryo_genes_toptable), replace = FALSE)

embryo_genes <- bind_rows("embryo" = embryo_genes_toptable, "random" = embryo_genes_random, .id = "subset") %>%
  pivot_longer(starts_with("log2"), names_to = "usp9x", values_to = "log2FC")

embryo_genes %>% filter(grepl("low48", usp9x)) %>% #as violin
  ggplot(aes(subset, log2FC, fill = subset)) +
  geom_hline(yintercept = 0) +
  geom_violin() + geom_boxplot(width = 0.15, fill = "white") + 
  scale_fill_manual(values = c("goldenrod2", "grey")) +
  ylim(c(-2.5, 3.5)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), axis.text = element_text(color = "black", size = "12"), axis.text.x = element_text(angle = 90), 
        panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #+ 
#ggsave("RNAseq_e8.5_expr_violin.pdf", device = "pdf", dpi = 300, width = 3, height = 3, units = "in")

# stats
wilcox.test(embryo_genes_toptable$log2FC_low48, embryo_genes_random$log2FC_low48)


########################
# compare to Beccari et al. (Fig. 2e, Supp Fig. 3i)

beccari_7.5 <- read_excel("Beccari_embryo_DEseq2.xlsx", sheet = 1, col_names = c("GeneID", "Ens", "baseMean", "log2FC", "lfcSE", "stat", "pvalue", "padj"), skip = 1) %>%
  select(GeneID, baseMean, log2FC, pvalue, padj) 

beccari_8.5 <- read_excel("Beccari_embryo_DEseq2.xlsx", sheet = 2, col_names = c("GeneID", "Ens", "baseMean", "log2FC", "lfcSE", "stat", "pvalue", "padj"), skip = 1) %>% 
  select(GeneID, baseMean, log2FC, pvalue, padj) 

beccari <- inner_join(beccari_7.5, beccari_8.5, by = "GeneID", suffix = c("_7.5", "_8.5")) %>% 
  select(-matches("base|padj|pvalue"))

# plot only up genes
up_in_both_beccari <- inner_join(up_in_both, beccari) %>% transmute(GeneID, log2FC_7.5 = -log2FC_7.5, log2FC_8.5 = -log2FC_8.5)
up_in_both_beccari %>% melt() %>%
  ggplot(aes(variable, value)) +
  geom_hline(yintercept = 0, color = "black") +
  stat_boxplot(geom = "errorbar", width = 0.3) + geom_boxplot(fill = "goldenrod2", color = "black") + 
  ylim(c(-6, 4)) +
  theme_classic() + 
  theme(panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(axis.text.x=element_text(size=12, color="black", angle = 90, hjust = 0), axis.text.y=element_text(size=12, color="black"))

# plot only down genes
down_in_both_beccari <- inner_join(down_in_both, beccari) %>% transmute(GeneID, log2FC_7.5 = -log2FC_7.5, log2FC_8.5 = -log2FC_8.5)
down_in_both_beccari %>% melt() %>%
  ggplot(aes(variable, value)) +
  geom_hline(yintercept = 0, color = "black") +
  stat_boxplot(geom = "errorbar", width = 0.3) + geom_boxplot(fill = "#2E0854", color = "black") + 
  #ylim(c(-6, 4)) +
  theme_classic() + 
  theme(axis.text.x=element_text(size=12, color="black", angle = 90, hjust = 0), axis.text.y=element_text(size=12, color="black"))

#stats
de_in_both_beccari <- inner_join(de_in_both, beccari) %>% transmute(GeneID, direction, log2FC_7.5 = -log2FC_7.5, log2FC_8.5 = -log2FC_8.5)
de_in_both_beccari_up <- de_in_both_beccari %>% filter(direction == "up")
de_in_both_beccari_down <- de_in_both_beccari %>% filter(direction == "down")

t.test(de_in_both_beccari_down$log2FC_7.5, de_in_both_beccari_down$log2FC_8.5)
t.test(de_in_both_beccari_up$log2FC_7.5, de_in_both_beccari_up$log2FC_8.5)


#########################
# K27 coverage over DE genes (Fig. 1f)

# import multiBamSummary counts
scale <- c(32955958, 28259537, 30329242, 24293121, 25977939, 16167853, 17860721, 22826385, 15104926, 17147758, 16991088)/1000000

up_in_both_K27 <- read_tsv("up_in_both_K27_10kb_upst.txt", skip = 1, 
                        col_names = c("chr", "start", "end", "E6.5_input", "E6.5_1", "E6.5_2", "E6.5_3", "E7.5_input", "E7.5_1", "E7.5_2", "E8.5_input", "E8.5_1", "E8.5_2", "E8.5_3")) %>%
  select(-chr:-end) %>%
  do(sweep(., 2, scale, "/")) %>% 
  transmute(E6.5_1 = E6.5_1/E6.5_input, E6.5_2 = E6.5_2/E6.5_input, E6.5_3 = E6.5_3/E6.5_input,
            E7.5_1 = E7.5_1/E7.5_input, E7.5_2 = E7.5_2/E7.5_input,
            E8.5_1 = E8.5_1/E8.5_input, E8.5_2 = E8.5_2/E8.5_input, E8.5_3 = E8.5_3/E8.5_input) %>% 
  select(-contains("input")) 

down_in_both_K27 <- read_tsv("down_in_both_K27_10kb_upst.txt", skip = 1, 
                           col_names = c("chr", "start", "end", "E6.5_input", "E6.5_1", "E6.5_2", "E6.5_3", "E7.5_input", "E7.5_1", "E7.5_2", "E8.5_input", "E8.5_1", "E8.5_2", "E8.5_3")) %>%
  select(-chr:-end) %>%
  do(sweep(., 2, scale, "/")) %>% 
  transmute(E6.5_1 = E6.5_1/E6.5_input, E6.5_2 = E6.5_2/E6.5_input, E6.5_3 = E6.5_3/E6.5_input,
            E7.5_1 = E7.5_1/E7.5_input, E7.5_2 = E7.5_2/E7.5_input,
            E8.5_1 = E8.5_1/E8.5_input, E8.5_2 = E8.5_2/E8.5_input, E8.5_3 = E8.5_3/E8.5_input) %>% 
  select(-contains("input")) 

de_in_both_K27 <- bind_rows("up" = up_in_both_K27, "down" = down_in_both_K27, .id = "direction") %>% 
  mutate(E6.5 = log2((E6.5_1 + E6.5_2 + E6.5_3)/3), E7.5 = log2((E7.5_1 + E7.5_2)/2), E8.5 = log2((E8.5_1 + E8.5_2 + E8.5_3)/3)) %>% 
  select(-contains("_"))

de_in_both_K27 %>% melt %>% ggplot(aes(variable, value, fill = direction)) + 
  geom_hline(yintercept = 0, color = "black") +
  stat_boxplot(lwd = 0.8, geom = "errorbar", width = 0.3) + geom_boxplot(width = 0.5, size = 0.8, color = "black") + #geom_jitter(width = 0.1, alpha = 0.7) + 
  scale_fill_manual(values = c("#2E0854", "goldenrod2")) +
  ylab("log2(K27/input)") + ggtitle("genes de in ko") + 
  facet_grid(~direction) + theme(strip.background = element_rect(color="white", fill="transparent", size=1, linetype="solid")) +
  theme(axis.text.x=element_text(size=12, color="black", angle = 90, hjust = 0), axis.text.y=element_text(size=14, color="black")) +
  theme(panel.background = element_rect(fill = "transparent", color = "black", size = 1), plot.background = element_rect(fill = "transparent"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #+
  #ggsave("de_in_mutants_H3K27me3.pdf", device = "pdf", dpi = 300, width = 4.5, height = 3.5, units = "in")

# to plot up/down separately...
de_in_both_K27 %>% filter(direction == "up") %>% melt() %>% 
  ggplot(aes(variable, value)) +
  geom_hline(yintercept = 0, color = "black") +
  stat_boxplot(geom = "errorbar", width = 0.3) + geom_boxplot(fill = "goldenrod2", color = "black", width = 0.7) + #geom_jitter(width = 0.2, alpha = 0.5) + 
  ylab("log2(K27/input)") + ggtitle("genes up in ko") + 
  theme_classic() +
  theme(panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(axis.text.x=element_text(size=12, color="black", angle = 90, hjust = 0), axis.text.y=element_text(size=14, color="black"))

# stats
up_in_both_K27_long <- de_in_both_K27 %>% filter(direction == "up") %>% melt()

K27_up_stats <- de_in_both_K27 %>% filter(direction == "up")
wilcox.test(K27_up_stats$E6.5, K27_up_stats$E8.5) 
wilcox.test(K27_up_stats$E6.5, K27_up_stats$E7.5)

K27_down_stats <- de_in_both_K27 %>% filter(direction == "down")
wilcox.test(K27_down_stats$E6.5, K27_down_stats$E8.5)
wilcox.test(K27_down_stats$E6.5, K27_down_stats$E7.5)


#########################
# ChEA binding analysis (Fig. 2d, Supp Fig. 3h)
#########################

# read in Enrichr ChEA output - take top 15 factors
chea_up_genes <- read_tsv("chea_up_in_both_padj0.1.txt", col_names = c("Term", "Overlap", "pvalue", "padj", "old_p", "old_padj", "Odds_ratio", "Combined_Score", "Genes"), skip = 1) %>%
  filter(pvalue < 0.01) %>% distinct(Overlap, Genes, .keep_all = TRUE) %>% arrange(padj) %>% dplyr::slice(1:15) 

chea_up_genes2 <- chea_up_genes %>% 
  separate(Overlap, into = c("count", "total")) %>%
  select(-contains("old"), -Odds_ratio) %>% #separate(Term, into = c("ChIP_X", "details"), "_", extra = "merge") %>% #"Expt", "Cells", "Species"
  mutate_at(vars(pvalue, count, Combined_Score), as.numeric) %>% 
  #slice(2:51) %>% 
  group_by(Term) %>% #slice(which.max(count)) %>% # filter(count == max(count) & padj == min(padj)) %>%
  add_column(Direction = "up")

glimpse(chea_up_genes2)

chea_down_genes <- read_tsv("chea_down_in_both_padj0.1.txt", col_names = c("Term", "Overlap", "pvalue", "padj", "old_p", "old_padj", "Odds_ratio", "Combined_Score", "Genes"), skip = 1) %>%
  filter(pvalue < 0.01) %>% distinct(Overlap, Genes, .keep_all = TRUE) %>% arrange(padj) %>% dplyr::slice(1:15)

chea_down_genes2 <- chea_down_genes %>% 
  separate(Overlap, into = c("count", "total")) %>%
  select(-contains("old"), -Odds_ratio) %>% #separate(Term, into = c("ChIP_X", "details"), "_", extra = "merge") %>% #"Expt", "Cells", "Species"
  mutate_at(vars(count:Combined_Score), as.numeric) %>% 
  #group_by(Term) %>% slice(which.max(count)) %>% # filter(count == max(count) & padj == min(padj)) %>%
  add_column(Direction = "down")

# plot individually - top 15 terms
chea_up_genes2 %>% ggplot(aes(x = fct_reorder(tolower(Term), -pvalue), y = -log10(pvalue), size = count)) + 
  geom_point(color = "goldenrod2") + coord_flip() +  
  theme(
    panel.background = element_rect(fill = "transparent", color = "black", size = 0.5), plot.background = element_rect(fill = "transparent"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 12, color = "black"), axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_text(vjust = -0.35, size = 14), axis.title.y = element_text(vjust = 0.35, size = 14),
    legend.background = element_blank(), legend.position = "bottom", legend.title=element_blank()) #+
  #ggsave("chea_up_padj0.1.pdf", device = "pdf", dpi = 300, width = 7, height = 4, units = "in")

chea_down_genes2 %>% ggplot(aes(x = fct_reorder(tolower(Term), -pvalue), y = -log10(pvalue), size = count)) + 
  geom_point(color = "#2E0854") + coord_flip() + 
  theme(
    panel.background = element_rect(fill = "transparent", color = "black", size = 0.5), plot.background = element_rect(fill = "transparent"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 12, color = "black"), axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_text(vjust = -0.35, size = 14), axis.title.y = element_text(vjust = 0.35, size = 14),
    legend.background = element_blank(), legend.position = "bottom", legend.title=element_blank()) #+
  #ggsave("chea_down_padj0.1.pdf", device = "pdf", dpi = 300, width = 7, height = 4, units = "in")

#######
# fin #
#######
