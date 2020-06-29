# code to produce scaled profile plots from output of deeptools computeMatrix / plotProfile --outFileNameData 
# K27 coverage over genes DE in  E8.5 Usp9x mutants (Fig. 2f, Supp Fig. 3j)

library(tidyverse)

#setwd("~/K27_invivo/")

###########################
# first genes upregulated in ko
up_genes_m <- read.delim("up_in_both_K27_bs100.txt", skip = 1, header = F)
up_genes_df <- as.data.frame(t(as.matrix(up_genes_m[,-c(1,2)])))
colnames(up_genes_df) <- c("bins", "e6.5_1", "e6.5_2", "e6.5_3", "e7.5_1", "e7.5_2", "e8.5_1", "e8.5_2", "e8.5_3")

up_genes_e6.5 <- as_tibble(up_genes_df) %>% na.omit() %>% 
  select(1:4) %>% 
  gather(key = "variable", value = "tags", 2:4) %>%
  group_by(bins) %>% 
  summarize_at(vars(tags), list(mean, min, max)) %>% 
  ungroup() %>%
  dplyr::rename(tags = fn1, ymin = fn2, ymax = fn3) 

up_genes_e7.5 <- as_tibble(up_genes_df) %>% na.omit() %>% 
  select(1,5,6) %>% 
  gather(key = "variable", value = "tags", 2:3) %>%
  group_by(bins) %>% 
  summarize_at(vars(tags), list(mean, min, max)) %>% 
  ungroup() %>%
  dplyr::rename(tags = fn1, ymin = fn2, ymax = fn3)

up_genes_e8.5 <- as_tibble(up_genes_df) %>% na.omit() %>% 
  select(1,7:9) %>% 
  gather(key = "variable", value = "tags", 2:3) %>%
  group_by(bins) %>% 
  summarize_at(vars(tags), list(mean, min, max)) %>% 
  ungroup() %>%
  dplyr::rename(tags = fn1, ymin = fn2, ymax = fn3)

up_genes_plot <-  bind_rows("e6.5" = up_genes_e6.5, "e7.5" = up_genes_e7.5, "e8.5" = up_genes_e8.5, .id = "stage")
head(up_genes_plot)

# curve fit with loess smoothing - here the se is of different points
ggplot(up_genes_plot, aes(x = bins, y = tags, color = stage, fill = stage)) +
  geom_vline(xintercept = 100, linetype = "dashed", color = "darkgrey") +
  geom_smooth(span = 0.08, alpha = 0.2) +
  scale_color_manual(values = c("#1F9FB6", "darkblue", "#904abc")) + 
  scale_fill_manual(values = c("#1F9FB6", "#8B74BD", "#663399")) + 
  ylim(c(-0.5, 1.5)) +
  theme(panel.background = element_rect(fill = "transparent", color = "black", size = 1), plot.background = element_rect(fill = "transparent"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(axis.text.x=element_text(size=12, color="black"), axis.text.y=element_text(size=12, color="black")) #+
  #ggsave("K27_up_genes.pdf", device = "pdf", dpi = 300, width = 4, height = 2, units = "in")


############
# down genes
down_genes_m <- read.delim("down_in_both_K27_bs100.txt", skip = 1, header = F)
down_genes_df <- as.data.frame(t(as.matrix(down_genes_m[,-c(1,2)])))
colnames(down_genes_df) <- c("bins", "e6.5_1", "e6.5_2", "e6.5_3", "e7.5_1", "e7.5_2", "e8.5_1", "e8.5_2", "e8.5_3")

down_genes_e6.5 <- as_tibble(down_genes_df) %>% na.omit() %>% select(1:4) %>% gather(key = "variable", value = "tags", 2:4) %>%
  group_by(bins) %>% summarize_at(vars(tags), list(mean, min, max)) %>% ungroup() %>%
  dplyr::rename(tags = fn1, ymin = fn2, ymax = fn3) 

down_genes_e7.5 <- as_tibble(down_genes_df) %>% na.omit() %>% select(1,5,6) %>% gather(key = "variable", value = "tags", 2:3) %>%
  group_by(bins) %>% summarize_at(vars(tags), list(mean, min, max)) %>% ungroup() %>%
  dplyr::rename(tags = fn1, ymin = fn2, ymax = fn3)

down_genes_e8.5 <- as_tibble(down_genes_df) %>% na.omit() %>% select(1,7:9) %>% gather(key = "variable", value = "tags", 2:3) %>%
  group_by(bins) %>% summarize_at(vars(tags), list(mean, min, max)) %>% ungroup() %>%
  dplyr::rename(tags = fn1, ymin = fn2, ymax = fn3)

down_genes_plot <-  bind_rows("e6.5" = down_genes_e6.5, "e7.5" = down_genes_e7.5, "e8.5" = down_genes_e8.5, .id = "stage")

# curve fit with loess smoothing - here the se is of different points
ggplot(down_genes_plot, aes(x = bins, y = tags, color = stage, fill = stage)) +
  geom_vline(xintercept = 100, linetype = "dashed", color = "darkgrey") +
  geom_smooth(span = 0.08, alpha = 0.2) +
  scale_color_manual(values = c("#1F9FB6", "darkblue", "#904abc")) + 
  scale_fill_manual(values = c("#1F9FB6", "#8B74BD", "#663399")) + 
  ylim(c(-0.5, 1.5)) +
  theme(panel.background = element_rect(fill = "transparent", color = "black", size = 1), plot.background = element_rect(fill = "transparent"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(axis.text.x=element_text(size=12, color="black"), axis.text.y=element_text(size=12, color="black")) #+
  #ggsave("K27_down_genes.pdf", device = "pdf", dpi = 300, width = 4, height = 2, units = "in")


