#!/usr/bin/env R
# coding utf-8

###################### Package Loading #######################
pacman::p_load(
    "Seurat", "ggplot2", "ggsci", "reshape2", "scales",
    "monocle", "viridis", "patchwork", "ggpubr", "ggrepel",
    "msigdbr", "fgsea", "tidyverse", "bseqsc", "ggsignif",
    "ggrepel", "pheatmap", "RColorBrewer", "ggpmisc", "ggplotify",
    "circlize", "zoo", "nichenetr", "tidyverse", "GSVA", "magrittr",
    "ggforce", "dplyr", "ggpointdensity", "Nebulosa", "escape",
    "clusterProfiler", "foreach"
)
setwd("/work/xiaxy/work/RA/NC")
source("/work/xiaxy/work/RA/NC/Dependent.R")
source("/work/xiaxy/work/RA/NC/Hallmarker_analysis.R")
source("/work/xiaxy/work/RA/NC/Trajactory_analysis.R")
source("/work/xiaxy/work/RA/NC/Transcription_factor_analysis.R")
#############################################################


###################### Read input ###########################
all <- readRDS("data_rds/input.Rds")
macro <- readRDS("data_rds/macrophage.Rds")
tcell <- readRDS("data_rds/tcell.Rds")
ACR20 <- readRDS("single_ACR20.Rds")
#############################################################


################## Diagnostric signature ####################
all <- AddModuleScore(all, features = list(diagnostic_signature), name = "signature") # nolint
sig_score <- aggregate(all$signature1, list(all$sampleid), mean) %>%
    mutate(group = factor(c(rep("OA-BT", 3), rep(c("RA-BT", "RA-AT"), 6)), levels = c("OA-BT", "RA-BT", "RA-AT"))) %>% # nolint
    mutate(drug = factor(c(rep("None", 3), rep("Adalimumab", 6), rep("Tofacitinib", 6)), levels = c("None", "Adalimumab", "Tofacitinib"))) %>% # nolint
    mutate(patient = c(1, 2, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9)) %>%
    filter(group %in% c("OA-BT", "RA-BT"))

p <- ggplot(sig_score, aes(x = group, y = x)) +
    geom_boxplot(aes(color = group), outlier.size = 0, outlier.colour = "white") + # nolint
    geom_point(aes(fill = group),
        position = position_jitterdodge(jitter.width = 0.4),
        color = "black", shape = 21, size = 3
    ) +
    stat_compare_means(
        comparisons =
            list(c("OA-BT", "RA-BT")),
        label.y = c(0.6)
    ) + # nolint
    scale_color_manual(values = group_color3) +
    scale_fill_manual(values = group_color3) +
    theme_classic2() +
    labs(y = "Diagnostic signature score", x = "") +
    scale_y_continuous(
        expand = c(0, 0.01), limits = c(0.1, 0.7)
    ) +
    scale_x_discrete(label = c("OA", "RA")) +
    theme(
        axis.text.x = element_text(color = "black", size = 14, angle = 45, hjust = 1), # nolint
        axis.text.y = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 16),
        legend.title = element_blank(), legend.position = "none"
    )
p
ggsave("Fig6/signature_score.pdf", p, width = 2, height = 3.5)

## p2
bulk <- read.table("input/PMID35602512.txt", header = T, sep = "\t", row.names = 1) # nolint
data_set <- list(enriched = diagnostic_signature)
result <- gsva(
    as.matrix(bulk), data_set,
    min.sz = 5,
    kcdf = "Poisson", method = "ssgsea",
    mx.diff = TRUE, verbose = FALSE, parallel.sz = 30
)
mt <- data.frame(sample = colnames(bulk), value = as.numeric(result))
mt$group <- factor(c(rep("OA", 15), rep("RA", 9)), levels = c("OA", "RA"))

p <- ggplot(mt, aes(x = group, y = value)) +
    geom_boxplot(aes(color = group), outlier.size = 0, outlier.colour = "white") + # nolint
    geom_point(aes(fill = group),
        position = position_jitterdodge(jitter.width = 0.4),
        color = "black", shape = 21, size = 3
    ) +
    stat_compare_means(
        comparisons =
            list(c("OA", "RA")),
        label.y = c(3.8)
    ) +
    scale_color_manual(values = group_color3) +
    scale_fill_manual(values = group_color3) +
    theme_classic2() +
    labs(y = "Signature score", x = "") +
    scale_y_continuous(
        expand = c(0, 0.01), limits = c(2.5, 4)
    ) +
    theme(
        axis.text.x = element_text(color = "black", size = 14, angle = 45, hjust = 1), # nolint
        axis.text.y = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 16),
        legend.title = element_blank(), legend.position = "none"
    )
p
ggsave("Fig6/PMID35602512_sig_score.pdf", p, width = 2, height = 3.5)

## p3
bulk <- read.table("input/GSE55235.txt", header = T, sep = "\t", row.names = 1) # nolint
data_set <- list(enriched = diagnostic_signature)
result <- gsva(
    as.matrix(bulk), data_set,
    min.sz = 5,
    kcdf = "Poisson", method = "ssgsea",
    mx.diff = TRUE, verbose = FALSE, parallel.sz = 30
)
mt <- data.frame(sample = colnames(bulk), value = as.numeric(result))
mt$group <- c(rep("OA", 10), rep("RA", 10))
mt$group <- factor(mt$group, levels = c("OA", "RA"))

p <- ggplot(mt, aes(x = group, y = value)) +
    geom_boxplot(aes(color = group), outlier.size = 0, outlier.colour = "white") + # nolint
    geom_point(aes(fill = group),
        position = position_jitterdodge(jitter.width = 0.4),
        color = "black", shape = 21, size = 3
    ) +
    stat_compare_means(
        comparisons =
            list(c("OA", "RA")),
        label.y = c(3.8)
    ) + # nolint
    scale_color_manual(values = group_color3) +
    scale_fill_manual(values = group_color3) +
    theme_classic2() +
    labs(y = "Signature score", x = "") +
    scale_y_continuous(
        expand = c(0, 0.01), limits = c(2.5, 4)
    ) +
    theme(
        axis.text.x = element_text(color = "black", size = 14, angle = 45, hjust = 1), # nolint
        axis.text.y = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 16),
        legend.title = element_blank(), legend.position = "none"
    )
p
ggsave("Fig6/GSE55235_sig_score.pdf", p, width = 2, height = 3.5)

## p4
bulk <- read.table("input/GSE55457.txt", header = T, sep = "\t", row.names = 1) # nolint

data_set <- list(enriched = diagnostic_signature)
result <- gsva(
    as.matrix(bulk), data_set,
    min.sz = 5,
    kcdf = "Poisson", method = "ssgsea",
    mx.diff = TRUE, verbose = FALSE, parallel.sz = 30
)
mt <- data.frame(sample = colnames(bulk), value = as.numeric(result))
mt$group <- c(rep("OA", 10), rep("RA", 13))
mt$group <- factor(mt$group, levels = c("OA", "RA"))

p <- ggplot(mt, aes(x = group, y = value)) +
    geom_boxplot(aes(color = group), outlier.size = 0, outlier.colour = "white") + # nolint
    geom_point(aes(fill = group),
        position = position_jitterdodge(jitter.width = 0.4),
        color = "black", shape = 21, size = 3
    ) +
    stat_compare_means(
        comparisons =
            list(c("OA", "RA")),
        label.y = c(2.8)
    ) +
    scale_color_manual(values = group_color3) +
    scale_fill_manual(values = group_color3) +
    theme_classic2() +
    labs(y = "Signature score", x = "") +
    scale_y_continuous(
        expand = c(0, 0.01), limits = c(1.5, 3)
    ) +
    theme(
        axis.text.x = element_text(color = "black", size = 14, angle = 45, hjust = 1), # nolint
        axis.text.y = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 16),
        legend.title = element_blank(), legend.position = "none"
    )
p
ggsave("Fig6/GSE55457_sig_score.pdf", p, width = 2, height = 3.5)

## p5
bulk <- read.table("input/GSE1919.txt", header = T, sep = "\t", row.names = 1) # nolint

data_set <- list(enriched = diagnostic_signature)
result <- gsva(
    as.matrix(bulk), data_set,
    min.sz = 5,
    kcdf = "Poisson", method = "ssgsea",
    mx.diff = TRUE, verbose = FALSE, parallel.sz = 30
)
mt <- data.frame(sample = colnames(bulk), value = as.numeric(result))
mt$group <- c(rep("OA", 5), rep("RA", 5))
mt$group <- factor(mt$group, levels = c("OA", "RA"))

p <- ggplot(mt, aes(x = group, y = value)) +
    geom_boxplot(aes(color = group), outlier.size = 0, outlier.colour = "white") + # nolint
    geom_point(aes(fill = group),
        position = position_jitterdodge(jitter.width = 0.4),
        color = "black", shape = 21, size = 3
    ) +
    stat_compare_means(
        comparisons =
            list(c("OA", "RA")),
        label.y = c(4.3)
    ) + # nolint
    scale_color_manual(values = group_color3) +
    scale_fill_manual(values = group_color3) +
    theme_classic2() +
    labs(y = "Signature score", x = "") +
    scale_y_continuous(
        expand = c(0, 0.01), limits = c(2.5, 4.5)
    ) +
    theme(
        axis.text.x = element_text(color = "black", size = 14, angle = 45, hjust = 1), # nolint
        axis.text.y = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 16),
        legend.title = element_blank(), legend.position = "none"
    )
p
ggsave("Fig6/GSE1919_sig_score.pdf", p, width = 2, height = 3.5)

## p6
bulk <- read.table("input/GSE12021.txt", header = T, sep = "\t", row.names = 1) # nolint

data_set <- list(enriched = diagnostic_signature)
result <- gsva(
    as.matrix(bulk), data_set,
    min.sz = 5,
    kcdf = "Poisson", method = "ssgsea",
    mx.diff = TRUE, verbose = FALSE, parallel.sz = 30
)
mt <- data.frame(sample = colnames(bulk), value = as.numeric(result))
mt$group <- c(
    rep("Normal", 4), "RA", rep("OA", 2),
    rep("RA", 5), "OA", "RA", "RA", rep("OA", 6),
    "RA", "RA", "RA", "OA", "RA", rep("Normal", 5)
)
mt$group <- factor(mt$group, levels = c("Normal", "OA", "RA"))
mt <- subset(mt, group %in% c("OA", "RA"))

p <- ggplot(mt, aes(x = group, y = value)) +
    geom_boxplot(aes(color = group), outlier.size = 0, outlier.colour = "white") + # nolint
    geom_point(aes(fill = group),
        position = position_jitterdodge(jitter.width = 0.4),
        color = "black", shape = 21, size = 3
    ) +
    stat_compare_means(
        comparisons =
            list(c("OA", "RA")),
        label.y = c(2.6)
    ) + # nolint
    scale_color_manual(values = group_color3) +
    scale_fill_manual(values = group_color3) +
    theme_classic2() +
    labs(y = "Signature score", x = "") +
    scale_y_continuous(
        expand = c(0, 0.01), limits = c(1, 3)
    ) +
    theme(
        axis.text.x = element_text(color = "black", size = 14, angle = 45, hjust = 1), # nolint
        axis.text.y = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 16),
        legend.title = element_blank(), legend.position = "none"
    )
p
ggsave("Fig6/GSE12021_sig_score.pdf", p, width = 2, height = 3.5)

## p7
bulk <- read.table("input/Bulk_abundance.txt", header = T, sep = "\t", row.names = 1) # nolint
data_set <- list(enriched = diagnostic_signature)
result <- gsva(
    as.matrix(bulk), data_set,
    min.sz = 5,
    kcdf = "Poisson", method = "ssgsea",
    mx.diff = TRUE, verbose = FALSE, parallel.sz = 30
)
mt <- data.frame(sample = colnames(bulk), value = as.numeric(result))
mt$group <- c(rep("OA-BT", 5), rep("RA-BT", 14), rep("RA-AT", 10))
mt$group <- factor(mt$group, levels = c("OA-BT", "RA-BT", "RA-AT"))
mt <- subset(mt, group %in% c("OA-BT", "RA-BT"))

p <- ggplot(mt, aes(x = group, y = value)) +
    geom_boxplot(aes(color = group), outlier.size = 0, outlier.colour = "white") + # nolint
    geom_point(aes(fill = group),
        position = position_jitterdodge(jitter.width = 0.4),
        color = "black", shape = 21, size = 3
    ) +
    stat_compare_means(
        comparisons =
            list(c("OA-BT", "RA-BT")),
        label.y = c(2.1)
    ) +
    scale_color_manual(values = group_color3) +
    scale_fill_manual(values = group_color3) +
    theme_classic2() +
    labs(y = "Diagnostic signature score", x = "") +
    scale_y_continuous(
        expand = c(0, 0.01), limits = c(1.6, 2.2)
    ) +
    scale_x_discrete(label = c("OA", "RA")) +
    theme(
        axis.text.x = element_text(color = "black", size = 14, angle = 45, hjust = 1), # nolint
        axis.text.y = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 16),
        legend.title = element_blank(), legend.position = "none"
    )
p
ggsave("Fig6/validation_sig_score.pdf", p, width = 2, height = 3.5)

## p8
nm <- read.table("input/E-MTAB-8322.csv", header = T, sep = ",", row.names = 1) # nolint

data_set <- list(enriched = diagnostic_signature)
result <- gsva(
    as.matrix(nm), data_set,
    min.sz = 5,
    kcdf = "Poisson", method = "ssgsea",
    mx.diff = TRUE, verbose = FALSE, parallel.sz = 30
)
mat <- data.frame(name = colnames(result), x = as.numeric(result), group = sapply(strsplit(colnames(nm), "_", fixed = T), "[", 1)) # nolint
mt <- aggregate(mat$x, list(mat$group), mean)
merge_clinical <- read.table("input/E-MTAB-8322_clinical.txt", header = T, sep = "\t") # nolint
merge_tab <- merge(merge_clinical, mt, by.x = "sample", by.y = "Group.1")
merge_tab <- merge_tab[which(merge_tab$group != "UPA"), ]
merge_tab$group <- factor(merge_tab$group, levels = c("Healthy", "Naive", "Resistant", "Remission")) # nolint
merge_tab <- subset(merge_tab, group %in% c("Healthy", "Naive"))

p <- ggplot(merge_tab, aes(x = group, y = x)) +
    geom_boxplot(aes(color = group), outlier.size = 0, outlier.colour = "white") + # nolint
    geom_point(aes(fill = group),
        position = position_jitterdodge(jitter.width = 0.4),
        color = "black", shape = 21, size = 3
    ) +
    stat_compare_means(
        comparisons =
            list(c("Healthy", "Naive")), # nolint
        label.y = c(1.95)
    ) + # nolint
    scale_color_manual(values = c(group_color3, pal_aaas(alpha = 0.8)(2))) +
    scale_fill_manual(values = c(group_color3, pal_aaas(alpha = 0.8)(2))) +
    theme_classic2() +
    labs(y = "Diagnostic signature score", x = "") +
    scale_y_continuous(
        expand = c(0, 0.01), limits = c(1.6, 2)
    ) +
    theme(
        axis.text.x = element_text(color = "black", size = 14, angle = 45, hjust = 1), # nolint
        axis.text.y = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 16),
        legend.title = element_blank(), legend.position = "none"
    )
p
ggsave("Fig6/scRNA_sig1_score.pdf", p, width = 2, height = 3.7)
#############################################################


############### Disease severity signature ##################
## p1
data_set <- list(enriched = disease_severity_signature)
nm <- read.table("input/E-MTAB-8322.csv", header = T, sep = ",", row.names = 1) # nolint
result <- gsva(
    as.matrix(nm), data_set,
    min.sz = 5,
    kcdf = "Poisson", method = "ssgsea",
    mx.diff = TRUE, verbose = FALSE, parallel.sz = 30
)
mat <- aggregate(as.numeric(result), list(sapply(strsplit(colnames(nm), "_", fixed = T), "[", 1)), mean) # nolint
merge_clinical <- read.table("input/E-MTAB-8322_clinical.txt", header = T, sep = "\t") # nolint
merge_tab <- merge(merge_clinical, mat, by.x = "sample", by.y = "Group.1")
merge_tab <- merge_tab[which(merge_tab$group != "UPA"), ]

p1 <- ggplot(merge_tab, aes(x = das28, y = x)) +
    geom_point(size = 3, color = "#D26964") +
    geom_smooth(se = FALSE, method = "lm", color = "black") +
    stat_cor(method = "pearson") +
    labs(x = "DAS28", y = "Signature score") +
    theme_classic() +
    theme(
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 12, color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 8, color = "black")
    )
p1
ggsave("Fig6/score1_DAS28_emtab8322.pdf", p1, width = 2.9, height = 2.8)

## p2
result <- gsva(
    as.matrix(all@assays$RNA@data), data_set,
    min.sz = 5,
    kcdf = "Poisson", method = "ssgsea",
    mx.diff = TRUE, verbose = FALSE, parallel.sz = 30
)
mat <- aggregate(as.numeric(result), list(all$sampleid), mean) # nolint
DAS28 <- c(NA, NA, NA, 5.29, 3.71, 4.43, 3.45, 4.15, 3.28, 5.16, 3.7, 4.22, 3.63, 6.93, 4.42) # nolint
merge_tab <- data.frame(value = mat$x, DAS28 = DAS28)

p1 <- ggplot(merge_tab[-(1:3), ], aes(x = DAS28, y = value)) +
    geom_point(size = 3, color = "#D26964") +
    geom_smooth(se = FALSE, method = "lm", color = "black") +
    stat_cor(method = "pearson") +
    labs(x = "DAS28", y = "Signature score") +
    theme_classic() +
    theme(
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 12, color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 8, color = "black")
    )
p1
ggsave("Fig6/score1_DAS28.pdf", p1, width = 2.9, height = 2.8)
#############################################################


##################### ACR20 signature #######################
## p1
ACR20 <- AddModuleScore(ACR20, features = list(ACR20_N_signature), name = "signature") # nolint
mat <- data.frame(group = ACR20$gg, value = ACR20$signature1)
ACR20a <- AddModuleScore(ACR20, features = list(ACR20_Y_signature), name = "signature") # nolint
mat$ss <- ACR20a$signature1
mat$ss_1 <- scale(mat$ss, center = -1.2)
mat$value_1 <- scale(mat$value, center = -1)
mat$tt <- mat$ss_1 / mat$value_1
mat$drug <- ACR20$drug

p <- ggplot(mat, aes(x = drug, y = tt)) +
    geom_boxplot(aes(fill = group)) +
    scale_fill_manual(values = c("#87B1C8", "#A1568E")) +
    labs(title = "", y = "Signature score\n(ACR20_Y/ACR20_N)") + # nolint
    theme_classic() +
    stat_compare_means(
        aes(group = group),
        label = "p.signif", label.y = 3.7
    ) +
    theme(
        # legend.position = "none",
        legend.title = element_blank(),
        panel.grid = element_blank(),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.x = element_text(
            angle = 45, hjust = 1,
            color = "black", size = 10
        ),
        plot.title = element_text(hjust = 0.5)
    ) +
    scale_y_continuous(limits = c(1.1 * min(mat$tt), 4)) # nolint
p
ggsave("Fig6/acr20_sig5.pdf", p, width = 4.5, height = 4.5)

## p2
p <- ggplot(mat, aes(x = group, y = tt)) +
    geom_boxplot(aes(fill = group)) +
    # geom_point(aes(color = group), position = "auto", shape = 21, size = 3) + # nolint
    scale_fill_manual(values = c("#87B1C8", "#A1568E")) +
    labs(title = "", y = "Signature score\n(ACR20_Y/ACR20_N)") + # nolint
    theme_classic() +
    stat_compare_means(
        label = "p.signif", label.y = 3.7
    ) +
    theme(
        legend.title = element_blank(),
        panel.grid = element_blank(),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.x = element_text(
            angle = 45, hjust = 1,
            color = "black", size = 10
        ),
        plot.title = element_text(hjust = 0.5)
    ) +
    scale_y_continuous(limits = c(1.1 * min(mat$tt), 4)) # nolint
p
ggsave("Fig6/acr20_sig5.pdf", p, width = 4.5, height = 4.5)
#############################################################
