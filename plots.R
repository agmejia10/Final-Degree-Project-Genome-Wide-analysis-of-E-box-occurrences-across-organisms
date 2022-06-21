library(ggplot2)
library(data.table)
library(RColorBrewer)
library(viridis)
library(ggpubr)
library(ggridges)
library(dplyr)

all_results <- fread('all_results.txt', data.table = FALSE)

colnames(all_results) <- c('id', 'specie', 'chr', 'start', 'end', 'ebox', 'observed', 'expected', 'ratio', 'dimer', 'observed_dimer', 'expected_dimer', 'ratio_dimer')

all_results <- all_results[all_results$specie != "epsteinbarr",]

specie.names <- c("amellifera"='A.mellifera', "celegans"='C.elegans', "chimpanzee"='P.troglodytes', "chicken"='G.gallus', 
                  "drosophmelanogaster"='D.melanogaster', "dog"='C.familiaris', "macac"='M.fascicularis', 
                  "human"='H.sapiens', "xtropicalis"='X.tropicalis', "orangutan"='Pongo', "scerevisiae"='S.cerevisiae', "opossum"='Didelphidae',
                  "zebrafish"='D.rerio', "pig"='S.scrofa', "gorilla"= 'Gorilla', "mouse"='M.musculus')

specie.names[all_results$specie] -> all_results$specie


all_results$methylated <- ifelse(all_results$specie  %in% c("D.melanogaster", "C.elegans", "S.cerevisiae"),
                                          "Non-methylating", "Methylating")

#file_table <- write.table(all_results)

summary(all_results)

e <- c("CAAATG", "CAATTG", "CAACTG", "CAAGTG", "CATATG", "CATATG", "CATCTG", "CATGTG", "CACCTG", "CACGTG", "CAGCTG")

all_results_filtered <- dplyr::filter(all_results, ebox %in% e)


#ggsave('Desktop/ratio_dimerebox.pdf', plot=p3, height=9, width=16)


colors_not_methylated <- brewer.pal(3, "Blues")
names(colors_not_methylated) <- c("D.melanogaster", "C.elegans", "S.cerevisiae")

colors_methylated <- magma(13)
names(colors_methylated) <- unique(all_results_filtered$specie)[!unique(all_results_filtered$specie) %in% c("D.melanogaster", "C.elegans", "S.cerevisiae")]

clr <- c(colors_methylated, colors_not_methylated)

orden <- c('S.cerevisiae', 'C.elegans', 'D.melanogaster','A.mellifera',
           'D.rerio', 'X.tropicalis', 'G.gallus', 'Didelphidae', 'C.familiaris', 'S.scrofa', 'M.musculus', 'M.fascicularis', 'Pongo', 'Gorilla', 'P.troglodytes', 'H.sapiens')

all_results_filtered$specie <- all_results_filtered$specie %>% factor(levels = orden)
all_results_filtered$methylated <- all_results_filtered$methylated %>% factor(levels = c("Non-methylating", "Methylating"))


p1 <- ggplot(all_results_filtered, aes(x = log2(ratio), col=specie), alpha=0.2) + 
  facet_grid(ebox~.) + 
  geom_density() + 
  geom_vline(xintercept=0)+ 
  theme_light() + 
  ggtitle("Density distribution of the E-boxes across whole-genome") + 
  scale_color_manual(values = clr) + 
  theme_pubr() + 
  xlim(-4, 4) +
  # scale_x_discrete(limits = orden) +
  theme(axis.text.x = )

ggsave('../density_ebox.pdf', plot=p1, height=9, width=16)

colores = c("Non-methylating" = "#8CC7A1", "Methylating" = "#92374D")

p2 <- ggplot(all_results_filtered, aes(x = specie, y =log2(ratio), fill=methylated), alpha=0.2) + 
  facet_grid(ebox~methylated, space="free", scales="free") + 
  geom_violin() + 
  geom_hline(yintercept = 0)+ 
  theme_light() + 
  # scale_color_manual(values = clr) + 
  # scale_fill_manual(values=clr) + 
  theme_pubr() + 
  ylim(-4,4)+
  scale_fill_manual(values=colores) + 
  # scale_x_discrete(limits = orden) +
  ggtitle("Violin plots of the distribution of the E-boxes across whole-genome") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ggsave('../plotviolin_speciexebox.pdf', plot=p2, height=9, width=16)


# Quitamos para hacer una versión para el texto del Paper

results_paper <- all_results_filtered %>% 
  dplyr::filter(specie %in% c('S.cerevisiae', 'C.elegans', 'D.melanogaster', 'D.rerio','H.sapiens')) %>% droplevels()

colours_filtered <- c('S.cerevisiae' = "#283044", 'C.elegans' = "#81F4E1", 'D.melanogaster' = "#56CBF9", 'D.rerio' = "#F19953",'H.sapiens' = "#FF686B")

# para el paper 

# density
p3 <- ggplot(results_paper, aes(x = log2(ratio), col=specie), alpha=0.2) + 
  facet_grid(ebox~.) + 
  geom_density() + 
  geom_vline(xintercept=0)+ 
  theme_light() + 
  ggtitle("Density distribution of the E-boxes across whole-genome") + 
  scale_color_manual(values = colours_filtered) +
  theme_pubr() + 
  xlim(-4, 4) +
  # scale_x_discrete(limits = orden) +
  theme(axis.text.x = )

ggsave('../density_ebox_paper.pdf', plot=p3, height=9, width=7)

p4 <- ggplot(results_paper, aes(x = specie, y =log2(ratio), fill=methylated), alpha=0.2) + 
  facet_grid(ebox~methylated, space="free", scales="free") + 
  geom_violin() + 
  geom_hline(yintercept = 0)+ 
  theme_light() + 
  theme_pubr() + 
  ylim(-4,4)+
  scale_fill_manual(values=colores) + 
  # scale_x_discrete(limits = orden) +
  ggtitle("Violin plots of the distribution of the E-boxes across whole-genome") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ggsave('../plotviolin_speciexebox_paper.pdf', plot=p2, height=9, width=7)
  
  
# Dimer density of all species

all_results$ratio_dimer %>% summary()
all_results$ratio %>% summary()


p5 <- ggplot(all_results, aes(x = log2(ratio_dimer), col=specie), alpha=0.2) + 
  facet_grid(dimer~.) + 
  geom_density() + 
  geom_vline(xintercept=0)+ 
  theme_light() + 
  ggtitle("Density distribution of the dimer across whole-genome") + 
  scale_color_manual(values = clr) + 
  theme_pubr() +
  xlim(-1.5, 1.5) +
  # scale_x_discrete(limits = orden) +
  theme(axis.text.x = )

ggsave('../density_dimer.pdf', plot=p5, height=9, width=16)


## 5 especies
p6 <- ggplot(results_paper, aes(x = log2(ratio_dimer), col=specie), alpha=0.2) + 
  facet_grid(dimer~.) + 
  geom_density() + 
  geom_vline(xintercept=0)+ 
  theme_light() + 
  ggtitle("Density distribution of the dimer across whole-genome") + 
  scale_color_manual(values = colours_filtered) +
  theme_pubr() + 
  xlim(-1.5, 1.5) +
  # scale_x_discrete(limits = orden) +
  theme(axis.text.x = )

ggsave('../density_dimer_paper.pdf', plot=p6, height=9, width=7)


# # Statistical analyses
# # ANOVA (or non-parametric) to know if there are significative differences.
# 
# x <- aov(all_results_filtered$ratio ~ all_results_filtered$specie)

# TukeyHSD(x)

# # ANOVA should not be used as it is influenced by the outliers. A non-parametric method can be used.

mean(c(2,3,5))
median(c(2,3,5))

mean(c(2,3,50))
median(c(2,3,50))

# There are significative differences between two species (D.melanogaster and A.mellifera)
kruskal.test(all_results_filtered$ratio ~ all_results_filtered$specie)

x <- pairwise.wilcox.test(all_results_filtered$ratio, all_results_filtered$specie, p.adjust.method = "fdr")

x

# Median comparative of CACGTG E-box

cacgtg <- all_results_filtered %>% filter(ebox == "CACGTG")

d_melanogaster <- cacgtg %>% filter(specie == "D.melanogaster")
a_meliffera <- cacgtg %>% filter(specie == "A.mellifera")
  
median(d_melanogaster$ratio)
median(a_meliffera$ratio)

wilcox.test(d_melanogaster$ratio, a_meliffera$ratio)

# Median comparative of all hexamers

d_melanogaster <- all_results_filtered %>% filter(specie == "D.melanogaster")
a_meliffera <- all_results_filtered %>% filter(specie == "A.mellifera")

median(d_melanogaster$ratio)
median(a_meliffera$ratio)

wilcox.test(d_melanogaster$ratio, a_meliffera$ratio)

# Wilcoxon Pairwise of only cacgtg
kruskal.test(cacgtg$ratio ~ cacgtg$specie)
x <- pairwise.wilcox.test(cacgtg$ratio, cacgtg$specie, p.adjust.method = "fdr")

tabla <- x$p.value %>% log(10)
tabla <- -tabla
write.csv(tabla, "../pvalues_cacgtg_wilcoxon_pairwise.csv")

# Heatmap of -log10(pvalue)

tabla[is.infinite(tabla)] <- 300

cairo_pdf("../heatmap1 cacgtg.pdf", width = 7, height = 5)
pheatmap::pheatmap(tabla, cluster_rows = F,
                   cluster_cols = F)
dev.off()

cairo_pdf("../heatmap2 cacgtg.pdf", width = 7, height = 5)
pheatmap::pheatmap(tabla, cluster_rows = F,
                   cluster_cols = F,
                   color = c("darkblue", "red"),
                   breaks = c(0, -log10(0.05)))
dev.off()

orden <- c('S.cerevisiae', 'C.elegans', 'D.melanogaster','A.mellifera',
           'D.rerio', 'X.tropicalis', 'G.gallus', 'Didelphidae', 'C.familiaris', 'S.scrofa', 'M.musculus', 'M.fascicularis', 'Pongo', 'Gorilla', 'P.troglodytes', 'H.sapiens')

colores = c("Non-methylating" = "#8CC7A1", "Methylating" = "#92374D")

p7 <- ggplot(all_results, aes(x = specie, y = log2(ratio_dimer), fill=methylated)) + 
  geom_boxplot() + 
  ggtitle("Boxplot distribution of the dimer across whole-genome") + 
  scale_fill_manual(values = colores) + 
  theme_pubr() +
  scale_x_discrete(limits = orden) + coord_flip()
  theme(axis.text.x = )

ggsave('../boxplot_dimer.pdf', plot=p7, height=9, width=16)
