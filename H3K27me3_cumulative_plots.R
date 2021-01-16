library(tidyverse)
library(viridis)
library(viridisLite)
library(pheatmap)
library(sva)
library(ggrepel)
library(reshape2)

#setwd("~/Documents/")
# REVISIONS: code now includs no-auxin samples ("neg")

usp9x_filt <- read_tsv("Usp9x_raw_counts_filt.txt")
spike <- c(1, 0.64, 0.50, 0.53, 0.50, 0.77, 0.84)
usp9x_norm <- usp9x_filt %>% do(sweep(usp9x_filt[,-c(1:3)], 2, spike, "*")) 

usp9x_norm <- usp9x_norm %>% transmute(neg_a = neg_a/inp_a, neg_b = neg_b/inp_a, low_8a = low_8a/inp_a, low_8b = low_8b/inp_a, 
                                       high_8a = high_8a/inp_a, high_8b = high_8b/inp_a) %>% 
  select(-starts_with("inp")) %>% 
  add_column(chr = usp9x_filt$chr, start = usp9x_filt$start, end = usp9x_filt$end, .before = 1)


# now van Mierlo data
vm_filt <- read_tsv("vm_raw_counts_filt.txt")
vm_spike <- c(1, 0.74, 0.70, 1, 0.43, 0.45)
vm_norm <- vm_filt %>% do(sweep(vm_filt[,-c(1:3)], 2, vm_spike, "*")) 
vm_norm <- vm_norm %>% transmute(rep1_2i = rep1_2i/input_2i, rep2_2i = rep2_2i/input_2i, rep1_fbs = rep1_fbs/input_fbs, rep2_fbs = rep2_fbs/input_fbs) %>% 
  select(-starts_with("inp")) %>% 
  mutate(chr = vm_filt$chr, start = vm_filt$start, end = vm_filt$end)

embryo_filt <- read_tsv("embryo_raw_counts_filt.txt")
# only through E7.5
#  based on *deduplicated* bam file size...
scale <- c(32955958, 28259537, 30329242, 24293121, 25977939, 16167853, 17860721, #e6.7-e7.5
  6471401, 8884981, 18428579, 17028823, 15348802, 22937054, 22183869, 11512529, 14678290, 7835380, 20373904, 18715173, 8297461, 18735604, 20239066, #2c - ICM
  8328612, 18551283, 19835285)/1000000 #morula

embryo_norm <- embryo_filt %>% do(sweep(embryo_filt[,-c(1:3)], 2, scale, "/")) %>% 
  transmute(rep1_2c = rep1_2c/inp_2c, rep2_2c = rep2_2c/inp_2c, rep3_2c = rep3_2c/inp_2c, rep1_4c = rep1_4c/inp_4c, rep2_4c = rep2_4c/inp_4c, rep3_4c = rep3_4c/inp_4c, 
            rep1_8c = rep1_8c/inp_8c, rep2_8c = rep2_8c/inp_8c, rep3_8c = rep3_8c/inp_8c, ICM_rep1 = ICM_rep1/ICM_input, ICM_rep2 = ICM_rep2/ICM_input,
            mor_rep1 = mor_rep1/mor_input, mor_rep2 = mor_rep2/mor_input,
            e6.5_rep1 = e6.5_rep1/e6.5_input, e6.5_rep2 = e6.5_rep2/e6.5_input, e6.5_rep3 = e6.5_rep3/e6.5_input, e7.5_rep1 = e7.5_rep1/e7.5_input, e7.5_rep2 = e7.5_rep2/e7.5_input) %>% 
  select(-contains("inp")) %>% 
  add_column(chr = embryo_filt$chr, start = embryo_filt$start, end = embryo_filt$end, .before = 1)

# join datasets
usp9x_embryo <- inner_join(usp9x_norm, embryo_norm, by = c("chr", "start", "end")) %>% select(-chr:-end)
usp9x_embryo_vm <- inner_join(usp9x_norm, embryo_norm, by = c("chr", "start", "end")) %>% 
  inner_join(vm_norm, by = c("chr", "start", "end")) %>% 
  select(-chr:-end)

##########
# ----- PCA plot

# assign PCA plot features
names <- colnames(usp9x_embryo_vm)

colors <- c("grey", "grey", "lightblue", "lightblue", "salmon1", "salmon1", 
            "red", "red", "red", 
            "orange", "orange", "orange", "gold1", "gold1", "gold1", 
            "palegreen3", "palegreen3", "darkorange", "darkorange", 
            "#1F9FB6", "#1F9FB6", "#1F9FB6", "darkblue","darkblue",
            "limegreen", "limegreen", "blue", "blue")


# ----- correct for batch effects (ESC vs embryo source)
condition <- factor(names)
batch <- factor(c(rep("1", times = 6), rep("2", times = 18), rep("1", times = 4))) # define batches
mod <- model.matrix(~1, data = condition) #no other adjustment variables aside from batch; fit an intercept matrix

usp9x_embryo_vm_m <- usp9x_embryo_vm %>% as.matrix()  
combat_cov <- sva::ComBat(dat = usp9x_embryo_vm_m, batch = batch, mod = mod, par.prior = TRUE, prior.plots = FALSE)

# -- plot PCA
pcaC <- prcomp(t(combat_cov), center = T) 
scoresC <- as.data.frame(pcaC$x)

eigs <- pcaC$sdev^2
percentVar <- round(100 * eigs/sum(eigs), 1)

ggplot(data = scoresC, aes(x = PC1, y = PC2, label = names)) +
  geom_point(color = colors, size = 2, alpha = 0.8) +
  geom_hline(yintercept = 0, color = "black") + 
  geom_vline(xintercept = 0, color= "black") + 
  #geom_text(color = colors, size = 4, hjust = -0.1, vjust = 0) + 
  #geom_text_repel(color = colors, nudge_x = -1, direction = "both", point.padding = 0.5) + 
  theme(axis.text.x = element_text(size = 12, color = "black"), axis.text.y = element_text(size = 12, color = "black")) +
  theme(axis.title.x = element_text(vjust = -0.35, size = 14, color = "black"), axis.title.y = element_text(vjust = 0.35, size = 14, color = "black")) +
  theme(panel.background = element_rect(fill = "transparent", color = "black", size = 1), plot.background = element_rect(fill = "transparent"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  labs(x = paste0("PC1 (",percentVar[1],"% of variance)"), y = paste0("PC2 (",percentVar[2],"% of variance)"), title = "PCA, batch-corrected") +
  theme(plot.title = element_text(size = 20, lineheight = 0.8, vjust = 2)) +
  xlim(c(-200, 175)) + ylim(c(-120,120)) + #create white space for labeling points
  ggsave("revision_plots/PCA_with_neg.pdf", device = "pdf", dpi = 300, width = 3, height = 3, units = "in")

# just to make sure all the sample labels are correct...
library(factoextra)
fviz_pca_ind(pcaC)

######################################
######### CUMULATIVE PLOTS ###########

# first log2 transform
usp9x_norm2 <- usp9x_norm %>% mutate(neg = (log2(neg_a + 0.5) + log2(neg_b + 0.5))/2, low_8 = (log2(low_8a + 0.5) + log2(low_8b + 0.5))/2,
                                     high_8 = (log2(high_8a + 0.5) + log2(high_8b + 0.5))/2) %>% select(-neg_a:-high_8b)

embryo_norm2 <- embryo_norm %>% select(-contains("inp")) %>% 
  mutate(x2c = (log2(rep1_2c + 0.5) + log2(rep2_2c + 0.5) + log2(rep3_2c + 0.5))/3, x4c = (log2(rep1_4c + 0.5) + log2(rep2_4c + 0.5) + log2(rep3_4c + 0.5))/3, x8c = (log2(rep1_8c + 0.5) + log2(rep2_8c + 0.5) + log2(rep3_8c + 0.5))/3,
         ICM = (log2(ICM_rep1 + 0.5) + log2(ICM_rep2 + 0.5))/2, morula = (log2(mor_rep1 + 0.5) + log2(mor_rep2 + 0.5))/2, TE = (log2(TE_rep1 + 0.5) + log2(TE_rep1 + 0.5))/2,
         e6.5 = (log2(e6.5_rep1 + 0.5) + log2(e6.5_rep2 + 0.5) + log2(e6.5_rep3 + 0.5))/3, e7.5 = (log2(e7.5_rep1 + 0.5) + log2(e7.5_rep2 + 0.5))/2) %>% 
  select(-matches("rep"))

vm_norm2 <- vm_norm %>% mutate(vm_2i = (log2(rep1_2i + 0.5) + log2(rep2_2i + 0.5))/2, vm_fbs = (log2(rep1_fbs + 0.5) + log2(rep2_fbs + 0.5))/2) %>% select(-matches("rep"))

# rearrange for plots
usp9x_norm2_plot <- usp9x_norm2 %>% select(-chr:-end) %>% gather("sample", "fc", neg, low_8, high_8) 

embryo_norm2_plot <- gather(embryo_norm2, "sample", "fc", x2c, x4c, x8c, morula, ICM, e6.5, e7.5)

vm_norm2_plot <- gather(vm_norm2, "sample", "fc", vm_2i, vm_fbs)

# stats
ks.test(usp9x_norm2$low_8, usp9x_norm2$high_8) # p-value < 2.2e-16

# define color vectors 
embryo_colors <- c("darkblue", "purple", "palegreen3", "darkorange", "red", "orange", "gold1", "darkred")
usp9x_8_colors <- c("salmon1", "lightblue", "grey")

# cumulative plots
ggplot(usp9x_norm2_plot, aes(fc, color = sample)) + 
  geom_vline(xintercept = 0, color = "black") +
  xlim(-2,2) +
  stat_ecdf(geom = "step") + xlab("log2(avg enrichment/input)") + ylab("Cumulative fraction") + 
  theme(axis.text.x = element_text(size = 12, color = "black"), axis.text.y = element_text(size = 12, color = "black"),
        panel.background = element_rect(fill = "transparent", color = "black", size = 1), plot.background = element_rect(fill = "transparent"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_color_manual(values = usp9x_8_colors)

ggplot(embryo_norm2_plot, aes(fc, color = sample)) + 
  xlab("log2(avg fold-change/input)") + ylab("Cumulative fraction") + 
  scale_color_manual(values = embryo_colors) + xlim(-2, 3) + geom_vline(xintercept = 0, color = "black") +
  stat_ecdf(geom = "step") + 
  theme(axis.text.x = element_text(size = 12, color = "black"), axis.text.y = element_text(size = 12, color = "black")) +
  #theme(axis.title.x = element_text(vjust = -0.35, size = 14, color = "black"), axis.title.y = element_text(vjust = 0.35, size = 14, color = "black")) +
  theme(panel.background = element_rect(fill = "transparent", color = "black", size = 1), plot.background = element_rect(fill = "transparent"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #+ 
  ggsave("embryo_fingerprints.pdf", device = "pdf", dpi = 300, width = 4, height = 3, units = "in")

ggplot(vm_norm2_plot, aes(fc, color = sample)) + 
  geom_vline(xintercept = 0, color = "black") +
  stat_ecdf(geom = "step") + xlab("log2(avg fold-change/input + 0.5)") + ylab("Cumulative fraction") + 
  scale_color_manual(values = c("limegreen", "blue")) + xlim(-2,2) +
  theme(axis.text.x = element_text(size = 12, color = "black"), axis.text.y = element_text(size = 12, color = "black"),
        panel.background = element_rect(fill = "transparent", color = "black", size = 1), plot.background = element_rect(fill = "transparent"), panel.grid.major = element_blank(), panel.grid.minor = element_blank())



# -----
# plot by replicates
usp9x_norm_reps <- usp9x_norm %>% transmute(neg_a = log2(neg_a + 0.5), neg_b = log2(neg_b + 0.5), 
                                            low_8a = log2(low_8a + 0.5), low_8b = log2(low_8b + 0.5), 
                                       high_8a = log2(high_8a + 0.5), high_8b = log2(high_8b + 0.5)) %>% select(-starts_with("inp"))

usp9x_norm_reps %>% gather("sample", "fc") %>% 
  ggplot(aes(fc, color = sample, linetype = sample)) + 
  geom_vline(xintercept = 0, color = "black") + xlim(-1.5, 1.5) +
  stat_ecdf(geom = "step") + xlab("log2(avg enrichment/input)") + ylab("Cumulative fraction") + 
  scale_color_manual(values = c("salmon1", "salmon1", "lightblue", "lightblue", "grey", "grey")) +
  scale_linetype_manual(values = c("solid", "dotted", "solid", "dotted", "solid", "dotted")) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(size = 12, color = "black"), axis.text.y = element_text(size = 12, color = "black")) +
  theme(axis.title.x = element_text(vjust = -0.35, size = 14), axis.title.y = element_text(vjust = 0.35, size = 14)) +
  theme(panel.background = element_rect(fill = "transparent", color = "black", size = 1), plot.background = element_rect(fill = "transparent"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggsave("~/Santos Lab Dropbox/Trisha Macrae/revisions_bioinf/revision_plots/Usp9x_fingerprints_neg.pdf", device = "pdf", dpi = 300, width = 3, height = 3, units = "in")


vm_norm_reps <- vm_norm %>% transmute(rep1_2i = log2(rep1_2i + 0.5), rep2_2i = log2(rep2_2i + 0.5), rep1_fbs = log2(rep1_fbs + 0.5), rep2_fbs = log2(rep2_fbs + 0.5)) %>% select(-starts_with("inp"))

vm_norm_reps %>% gather("sample", "fc") %>% 
  ggplot(aes(fc, color = sample, linetype = sample)) + 
  geom_vline(xintercept = 0, color = "black") + #xlim(-1.5, 1.5) +
  stat_ecdf(geom = "step") + xlab("log2(enrichment/input)") + ylab("Cumulative fraction") + 
  scale_color_manual(values = c("limegreen", "blue", "limegreen", "blue")) +
  scale_linetype_manual(values = c("solid", "solid", "dotted", "dotted")) +
  #theme(legend.position = "none") +
  theme(axis.text.x = element_text(size = 12, color = "black"), axis.text.y = element_text(size = 12, color = "black")) +
  #theme(axis.title.x = element_text(vjust = -0.35, size = 14), axis.title.y = element_text(vjust = 0.35, size = 14)) +
  theme(panel.background = element_rect(fill = "transparent", color = "black", size = 1), plot.background = element_rect(fill = "transparent"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #+
  #ggsave("fingerprints.pdf", device = "pdf", dpi = 300, width = 4, height = 3, units = "in")

# ----- STATS 
# test is Kolmogorov-Smirnov test
ks.test(usp9x_norm_reps$low_8a, usp9x_norm_reps$neg_a) #p-value < 2.2e-16
ks.test(usp9x_norm_reps$high_8a, usp9x_norm_reps$neg_a) #p-value < 2.2e-16
ks.test(usp9x_norm_reps$low_8b, usp9x_norm_reps$neg_b) #p-value < 2.2e-16
ks.test(usp9x_norm_reps$high_8b, usp9x_norm_reps$neg_b) #p-value < 2.2e-16

# old
ks.test(usp9x_norm_reps$low_8a, usp9x_norm_reps$high_8a) #p-value < 2.2e-16
ks.test(usp9x_norm_reps$low_8b, usp9x_norm_reps$high_8b) #p-value < 2.2e-16

ks.test(vm_norm_reps$rep1_2i, vm_norm_reps$rep1_fbs) #p-value < 2.2e-16
ks.test(vm_norm_reps$rep2_2i, vm_norm_reps$rep2_fbs) #p-value < 2.2e-16

embryo_norm2_stats <- embryo_norm2 %>% mutate(postimp = (e6.5 + e7.5)/2, preimp = (x2c + x4c + x8c + morula + ICM)/5) %>% select(postimp, preimp)
ks.test(embryo_norm2_stats[,1], embryo_norm2_stats[,2]) #p-value < 2.2e-16




##
# Supp Fig 4b: import bivalent (requires txt files with bivalent peaks or all K27 peaks, named as below)
bivalent <- read_tsv("bivalent_peaks_K27_coverage.txt", skip = 1, col_names = c("chr", "start", "end", "neg_a", "neg_b", "high_8a", "high_8b", "low_8a", "low_8b")) %>%
  select(-chr:-end) %>% filter_at(vars(neg_a, neg_b), all_vars(. > 1))
all_peaks <- read_tsv("all_peaks_K27_coverage.txt", skip = 1, col_names = c("chr", "start", "end", "neg_a", "neg_b", "high_8a", "high_8b", "low_8a", "low_8b")) %>%
  select(-chr:-end) %>% filter_at(vars(neg_a, neg_b), all_vars(. > 1)) 

spike2 <- c(0.64, 0.50, 0.77, 0.84, 0.53, 0.50)

#normalize, combine, and plot
bivalent_norm <- bivalent %>% do(sweep(., 2, spike2, "*")) %>% mutate(low_8 = log2((low_8a/neg_a + low_8b/neg_b)/2), high_8 = log2((high_8a/neg_a + high_8b/neg_b)/2 )) %>% 
  select(-matches("a|b")) %>% filter(low_8 != "Inf")
all_peaks_norm <- all_peaks %>% do(sweep(., 2, spike2, "*")) %>% mutate(low_8 = log2((low_8a/neg_a + low_8b/neg_b)/2), high_8 = log2((high_8a/neg_a + high_8b/neg_b)/2 )) %>% 
  select(-matches("a|b"))

set.seed(5)
sample_peaks <- sample_n(all_peaks_norm, nrow(bivalent_norm), replace = FALSE)

peaks_compare <- bind_rows("bivalent" = bivalent_norm, "xall" = all_peaks_norm, "random" = sample_peaks, .id = "peakset") 

peaks_compare %>% melt() %>% 
  mutate(variable = fct_relevel(.$variable, "high_8", "low_8")) %>% 
  ggplot(aes(peakset, value, fill = variable)) +
  scale_fill_manual(values = c("salmon1", "lightblue", "salmon1", "lightblue")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  #geom_violin(bw = 0.1, size = 0.8, width = 0.8) + geom_boxplot(width = 0.1, position = position_dodge(width = 0.8), outlier.shape = NA) +
  stat_boxplot(lwd = 0.5, geom = "errorbar", width = 0.4, position = position_dodge(width = 0.8), color = "black") +
  geom_boxplot(size = 0.5, width = 0.8, color = "black", outlier.shape = NA) +
  ylim(c(-0.75, 1.25)) +
  theme(axis.text.x = element_text(size = 12, color = "black"), axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(vjust = -0.35, size = 14), axis.title.y = element_text(vjust = 0.35, size = 14), 
        panel.background = element_rect(fill = "transparent", color = "black", size = 1), plot.background = element_rect(fill = "transparent"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #+ 
ggsave("~/Documents/***.pdf", device = "pdf", dpi = 300, width = 4.5, height = 3, units = "in")


