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
    "clusterProfiler"
)
setwd("/work/xiaxy/work/RA/NC")
source("/work/xiaxy/work/RA/NC/Dependent.R")
source("/work/xiaxy/work/RA/NC/Hallmarker_analysis.R")
source("/work/xiaxy/work/RA/NC/Trajactory_analysis.R")
source("/work/xiaxy/work/RA/NC/Transcription_factor_analysis.R")
#############################################################


###################### Read input ###########################
data_input <- readRDS("data_rds/input.Rds")
data <- readRDS("data_rds/macrophage.Rds")
#############################################################


################## Macrophage signature #####################
## p1
subset_data <- subset(data, idents = levels(data)[1:5])
signature <- GSEABase::getGmt("~/.sc_process/database/signature/TAM_cell_signature.gmt") # nolint
escape_result <- enrichIt(
    subset_data,
    gene.sets = signature,
    groups = 1000, cores = 30, min.size = 5
)
subset_data <- Seurat::AddMetaData(subset_data, escape_result)

select <- c("Classical_TIMs", "M1_Macrophage_Polarization", "M2_Macrophage_Polarization") # nolint
mat_data <- subset_data@meta.data[, select] # nolint
mt_data <- aggregate(mat_data, list(data$celltype), median)
rownames(mt_data) <- mt_data$Group.1
mt_data <- as.data.frame(t(mt_data[, -1]))

annotation_col <- data.frame(
    celltype = colnames(mt_data)
)
rownames(annotation_col) <- colnames(mt_data)
col1 <- celltype_color[1:5]
names(col1) <- colnames(mt_data)

ann_colors <- list(
    celltype = col1
)

bk <- c(seq(-1, -0.1, by = 0.01), seq(0, 1, by = 0.01))
pdf("Fig3/f_signature.pdf", width = 3.2, height = 0.8)
pheatmap(mt_data,
    show_colnames = F, annotation_col = annotation_col,
    annotation_colors = ann_colors,
    show_rownames = T, scale = "row", cluster_cols = T,
    cluster_rows = T, color = c(
        colorRampPalette(colors = viridis(100)[1:50])(length(bk) / 2), # nolint
        colorRampPalette(colors = viridis(100)[51:100])(length(bk) / 2)
    ),
    breaks = bk, annotation_legend = TRUE,
    legend_breaks = c(-1, 0, 1),
    fontsize = 6, annotation_names_col = FALSE,
    border_color = "white", treeheight_row = 8, treeheight_col = 8
)
dev.off()

## p2
mat_data <- mat_data %>%
    mutate(group = factor(paste(data$subtype, data$drug, sep = "_"))) %>%
    mutate(sampleid = data$sampleid) %>%
    mutate(celltype = data$celltype) %>%
    mutate(new_group = paste(celltype, sampleid, sep = "_")) %>%
    mutate(subtype = data$subtype)

mycomparisons <- list(
    c("OA_BT_None", "RA_BT_Adalimumab"), c("RA_BT_Adalimumab", "RA_AT_Adalimumab"), # nolint
    c("RA_BT_Tofacitinib", "RA_AT_Tofacitinib"), c("OA_BT_None", "RA_BT_Tofacitinib") # nolint
)

select_dr <- function(x = "SPP1+ Mac", y = "M1_Macrophage_Polarization", z = "M1") { # nolint
    mt <- subset(mat_data, celltype == x)
    mt[, y] <- mt[, y] / 1000

    p <- ggplot(mt, aes(x = group, y = get(y), color = group)) +
        geom_violin() +
        scale_color_manual(values = group_color5) + # nolint
        labs(x = paste0(x, "rophage"), y = "Pathway activity", title = z) +
        stat_summary(
            fun.y = median, geom = "point",
            fill = "black", , color = "black", shape = 21, size = 1.5
        ) +
        scale_y_continuous(limits = c(ifelse(min(mt[, y]) > 0, 1.8, -0.2), 8)) + # nolint
        geom_vline(xintercept = c(1.5, 3.5), linetype = "dashed", color = "gray") + # nolint
        stat_compare_means(
            comparisons = mycomparisons, label = "p.signif",
            label.y = c(1.01, 1.01, 1.01, 1.08) * max(mt[, y])
        ) +
        theme_bw() +
        theme(
            panel.grid = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.y = element_text(size = 10, color = "black"),
            axis.title.y = element_text(size = 12, color = "black"),
            legend.position = "none",
            plot.title = element_text(hjust = 0.5)
        )
}

final <- select_dr(x = "SPP1+ Mac", y = "M1_Macrophage_Polarization", z = "M1") + # nolint
    select_dr(x = "S100A12+ Mac", y = "M1_Macrophage_Polarization", z = "M1") +
    plot_layout(nrow = 1)

ggsave("Fig3/g_select.pdf", final, width = 6.2, height = 2.5)

## p3
mat_data_acr20 <- subset(mat_data, subtype == "RA_BT")
mat_data_acr20 <- mat_data_acr20 %>%
    mutate(ACR20 = ifelse(sampleid %in% c("RA_3_1", "RA_5_1"), "ACR20_N", "ACR20_Y")) %>% # nolint
    mutate(M1 = M1_Macrophage_Polarization / 1000) %>%
    group_by(ACR20, celltype) %>%
    mutate(med = median(M1, na.rm = T))
m1_spp1 <- subset(mat_data_acr20, celltype = "SPP1+ Mac")
m1_s100a12 <- subset(mat_data_acr20, celltype = "S100A12+ Mac")

p1 <- ggplot(m1_spp1, aes(x = ACR20, y = M1)) + # nolint
    geom_point(aes(color = ACR20), position = "auto", shape = 21, size = 1) +
    scale_color_manual(values = c("#87B1C8", "#A1568E")) +
    labs(x = "SPP1+ Macrophage", y = "Pathway activity", title = "M1\nsignature") + # nolint
    theme_classic() +
    stat_compare_means(
        comparisons = list(c("ACR20_N", "ACR20_Y")),
        label.y = 1.05 * max(m1_spp1$M1),
        size = 6, bracket.size = 0.2,
        label = "p.signif"
    ) +
    geom_segment(
        aes(
            x = as.numeric(factor(ACR20)) - 0.1, y = med,
            xend = as.numeric(factor(ACR20)) + 0.1, yend = med
        ),
        color = "black"
    ) +
    scale_y_continuous(
        limits = c(0, 8)
    ) +
    theme(
        legend.position = "none",
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 12, color = "black"),
        plot.title = element_text(hjust = 0.5, color = "black", )
    )

p2 <- ggplot(m1_s100a12, aes(x = ACR20, y = M1)) + # nolint
    geom_point(aes(color = ACR20), position = "auto", shape = 21, size = 1) +
    scale_color_manual(values = c("#87B1C8", "#A1568E")) +
    labs(x = "S100A12+ Macrophage", y = "Pathway activity", title = "M1\nsignature") + # nolint
    theme_classic() +
    stat_compare_means(
        comparisons = list(c("ACR20_N", "ACR20_Y")),
        label.y = 1.05 * max(m1_s100a12$M1),
        size = 6, bracket.size = 0.2,
        label = "p.signif"
    ) +
    geom_segment(
        aes(
            x = as.numeric(factor(ACR20)) - 0.1, y = med,
            xend = as.numeric(factor(ACR20)) + 0.1, yend = med
        ),
        color = "black"
    ) +
    scale_y_continuous(
        limits = c(0.85 * min(m1_s100a12$M1), 1.2 * max(m1_s100a12$M1)) # nolint
    ) +
    theme(
        legend.position = "none",
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 12, color = "black"),
        plot.title = element_text(hjust = 0.5, color = "black", )
    )

final <- p1 + p2 + plot_layout(nrow = 1, guides = "collect")
ggsave("Fig3/Mac_M1.pdf", final, width = 4.2, height = 2.7)

## p4
bulk_data <- read.table("input/Bulk_abundance.txt", header = T, sep = "\t", row.names = 1) # nolint
data_set <- list(enriched = M1)
result <- gsva(
    as.matrix(bulk_data), data_set,
    min.sz = 5,
    kcdf = "Poisson", method = "ssgsea",
    mx.diff = TRUE, verbose = FALSE, parallel.sz = 30
)
mt <- data.frame(
    sample = colnames(bulk_data),
    value = as.numeric(result),
    group = factor(c(rep("OA-BT", 5), rep("RA-BT", 14), rep("RA-AT", 10)),
        levels = c("OA-BT", "RA-BT", "RA-AT")
    )
)
mt <- mt[-c(1, 17, 18), ] # remove extreme outliers
p <- ggplot(mt, aes(x = group, y = value)) +
    geom_boxplot(aes(color = group), outlier.size = 0, outlier.colour = "white") + # nolint
    geom_point(aes(fill = group),
        position = position_jitterdodge(jitter.width = 0.4),
        color = "black", shape = 21, size = 3
    ) +
    stat_compare_means(
        comparisons =
            list(c("OA-BT", "RA-BT"), c("RA-BT", "RA-AT")),
        label.y = c(1.8)
    ) +
    scale_color_manual(values = group_color3) +
    scale_fill_manual(values = group_color3) +
    theme_classic2() +
    scale_y_continuous(
        expand = c(0, 0.01), limits = c(1, 2.2)
    ) +
    labs(y = "Diagnostic signature score", x = "") +
    theme(
        axis.text.x = element_text(color = "black", size = 14, angle = 45, hjust = 1), # nolint
        axis.text.y = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 16),
        legend.title = element_blank(), legend.position = "none"
    )
p
ggsave("Fig3/m1_gsva.pdf", p, width = 2.3, height = 3.2)
#############################################################


################## Macrophage hallmarker ####################
RunEscape(data)

## p1
esc <- read.table("input/macrophage_hallmarker.txt", header = T, sep = "\t")
mat <- aggregate(esc, list(data$celltype), median)
rownames(mat) <- mat$Group.1
mat <- mat[, -1]
colnames(mat) <- gsub("HALLMARK_", "", colnames(mat))
mat <- as.data.frame(t(mat))[, 1:5]
annotation_col <- data.frame(
    celltype = colnames(mat)
)
rownames(annotation_col) <- colnames(mat)
col1 <- macrophage_color[1:5]
names(col1) <- colnames(mat)
ann_colors <- list(
    celltype = col1
)

bk <- c(seq(-1, -0.1, by = 0.01), seq(0, 1, by = 0.01))
pdf("Fig3/h_hall.pdf", width = 3.8, height = 4.5)
pheatmap(mat,
    show_colnames = F, annotation_col = annotation_col,
    annotation_colors = ann_colors,
    show_rownames = T, scale = "row", cluster_cols = T,
    cluster_rows = T, color = c(
        colorRampPalette(colors = brewer.pal(11, "RdYlBu")[11:6])(length(bk) / 2), # nolint
        colorRampPalette(colors = brewer.pal(11, "RdYlBu")[6:1])(length(bk) / 2)
    ),
    breaks = bk, annotation_legend = TRUE,
    legend_breaks = c(-1, 0, 1),
    fontsize = 6, annotation_names_col = FALSE,
    border_color = "white", treeheight_row = 8, treeheight_col = 8
)
dev.off()

## p2
esc <- read.table("input/macrophage_hallmarker.txt", header = T, sep = "\t")
esc <- esc[, c(4, 22:27, 35, 44, 45)]
mat <- aggregate(esc, list(paste(data$celltype, data$subtype, data$drug, sep = "_")), quantile) # nolint
colnames(mat) <- gsub("HALLMARK_", "", colnames(mat))

per_dr <- function(x = "ANGIOGENESIS", y = c("SPP1+ Mac", "S100A12+ Mac")) {
    uk <- as.data.frame(mat[, x])[, 2:4]
    uk$celltype <- sapply(strsplit(mat$Group.1, "_", fixed = T), "[", 1)
    uk$group <- paste(
        sapply(strsplit(mat$Group.1, "_", fixed = T), "[", 2),
        sapply(strsplit(mat$Group.1, "_", fixed = T), "[", 3),
        sapply(strsplit(mat$Group.1, "_", fixed = T), "[", 4),
        sep = "_"
    )
    uk$group <- factor(uk$group, levels = c(
        "OA_BT_None", "RA_BT_Adalimumab",
        "RA_AT_Adalimumab", "RA_BT_Tofacitinib", "RA_AT_Tofacitinib"
    ))
    uk <- subset(uk, celltype %in% y)
    uk$celltype <- factor(uk$celltype, levels = y)
    colnames(uk)[1:3] <- c("fr", "mi", "la")

    p <- ggplot(data = uk) +
        geom_ribbon(aes(x = group, ymin = fr, ymax = la, group = celltype, fill = celltype), alpha = 0.1, color = NA) + # nolint
        scale_fill_manual(values = macrophage_color[2:3]) +
        geom_line(aes(x = group, y = mi, group = celltype, color = celltype), linewidth = 1) + # nolint
        geom_point(aes(x = group, y = mi, color = celltype), size = 4) +
        scale_color_manual(values = macrophage_color[2:3]) +
        scale_x_discrete(
            limits = levels(uk$group), expand = c(0.04, 0.04),
            breaks = levels(uk$group), labels = c("OA-BT", "RA-BT", "RA-AT", "RA-BT", "RA-AT") # nolint
        ) + # nolint
        geom_vline(xintercept = c(1.5, 3.5), linetype = "dashed", size = 0.6, color = "lightgrey") + # nolint
        theme_bw() +
        annotate("text",
            y = 1.1 * min(uk$fr), x = 2.5, # 定义添加文本和其位置
            label = "Adalimumab",
            size = 3, color = "black"
        ) +
        annotate("text",
            y = 1.1 * min(uk$fr), x = 4.4, # 定义添加文本和其位置
            label = "Tofacitinib",
            size = 3, color = "black"
        ) +
        labs(x = "", y = "Pathway activity", title = x) +
        theme(
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 8), # nolint
            plot.title = element_text(hjust = 0.5, color = "black", size = 10),
            axis.text.y = element_text(color = "black", size = 8),
            axis.title.y = element_text(color = "black", size = 10),
            legend.title = element_blank()
        )
    return(p)
}
final <- per_dr(x = "INTERFERON_ALPHA_RESPONSE") + per_dr(x = "INTERFERON_GAMMA_RESPONSE") + # nolint
    per_dr(x = "TNFA_SIGNALING_VIA_NFKB") + per_dr(x = "INFLAMMATORY_RESPONSE") + # nolint
    per_dr(x = "ANGIOGENESIS") + per_dr(x = "HYPOXIA") +
    per_dr(x = "IL6_JAK_STAT3_SIGNALING") + per_dr(x = "IL2_STAT5_SIGNALING") +
    per_dr(x = "TGF_BETA_SIGNALING") +
    plot_layout(nrow = 3, guides = "collect")
ggsave("Fig3/j_select_hall.pdf", final, width = 10.5, height = 7.5)

## p3
esc <- read.table("input/macrophage_hallmarker.txt", header = T, sep = "\t") # nolint
colnames(esc) <- gsub("HALLMARK_", "", colnames(esc))
cell <- rownames(subset(
    data_input@meta.data,
    subtype == "RA_BT" & ACR20 %in% c("N", "Y") &
        subcell %in% c("SPP1+ Mac", "S100A12+ Mac")
))
sub <- subset(data_input, cells = cell)
pathway <- c(
    "INTERFERON_ALPHA_RESPONSE", "INTERFERON_GAMMA_RESPONSE",
    "TNFA_SIGNALING_VIA_NFKB", "INFLAMMATORY_RESPONSE",
    "ANGIOGENESIS", "HYPOXIA",
    "IL6_JAK_STAT3_SIGNALING", "TGF_BETA_SIGNALING"
)
esc_select <- esc[cell, pathway]
esc_select$group <- factor(sub$ACR20, levels = c("N", "Y"))
esc_select$celltype <- factor(sub$subcell, levels = c("SPP1+ Mac", "S100A12+ Mac")) # nolint

drawd <- function(gp) {
    p <- ggplot(esc_select, aes(x = celltype, y = get(gp))) +
        geom_boxplot(aes(fill = group),
            outlier.colour = "white", outlier.size = 0
        ) +
        theme_classic() +
        scale_fill_manual(values = c("#87B1C8", "#A1568E")) +
        labs(x = "", y = "Pathway activity", title = gp) + # nolint
        stat_compare_means(
            aes(group = group),
            label.y = max(esc_select[, gp]), size = 4,
            bracket.size = 0.2, label = "p.signif"
        ) +
        # scale_x_discrete(labels = c("ACR20_N", "ACR20_Y")) +
        scale_y_continuous(
            limits = c(min(esc_select[, gp]), 1.05 * max(esc_select[, gp])) # nolint
        ) +
        theme(
            axis.text = element_text(size = 10, color = "black"),
            axis.title = element_text(size = 12, color = "black"),
            plot.title = element_text(hjust = 0.5, color = "black", size = 10)
        )
    return(p)
}
final <- drawd(colnames(esc_select)[1]) + drawd(colnames(esc_select)[2]) +
    drawd(colnames(esc_select)[3]) + drawd(colnames(esc_select)[4]) +
    drawd(colnames(esc_select)[5]) + drawd(colnames(esc_select)[6]) +
    drawd(colnames(esc_select)[7]) + drawd(colnames(esc_select)[8]) +
    plot_layout(nrow = 2, guides = "collect")
ggsave("Fig3/mac_escape_acr20.pdf", final, width = 12, height = 6)
#############################################################


################## Trajactory analysis ######################
spp1_a <- subset(data, cells = rownames(data@meta.data)[which(data$celltype == "SPP1+ Mac" & data$drug == "Adalimumab")]) # nolint
spp1_t <- subset(data, cells = rownames(data@meta.data)[which(data$celltype == "SPP1+ Mac" & data$drug == "Tofacitinib")]) # nolint
RunMonocle(spp1_a, outputdir = "/work/xiaxy/work/RA/NC/input/", name = "spp1_a.Rds") # nolint
RunMonocle(spp1_t, outputdir = "/work/xiaxy/work/RA/NC/input/", name = "spp1_t.Rds") # nolint

## Tofacitinib
monocds <- readRDS("input/spp1_t.Rds")
monocds@phenoData@data$newstate <- ifelse(monocds@phenoData@data$State %in% 2:4, "State1", # nolint
    ifelse(monocds@phenoData@data$State == 1, "State2", "State3")
)

p <- plot_cell_trajectory(
    monocds,
    cell_size = 0.2, color_by = "newstate",
    show_tree = T, show_branch_points = F
) +
    scale_color_manual(values = c(pal_aaas()(5)[c(2, 1)], "#EFC07C")) +
    theme(
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.line.x = NULL,
        axis.line.y = NULL
    )
ggsave("Fig3/spp1_tof_mono_state.tiff", p, device = "tiff", width = 2.6, height = 2.4) # nolint

monocds@phenoData@data$acr20 <- ifelse(monocds@phenoData@data$sampleid %in% c("RA_5_1", "RA_5_2"), "N", "Y") # nolint
monocds@phenoData@data$subtype <- factor(monocds@phenoData@data$subtype, levels = c("RA_BT", "RA_AT")) # nolint
p <- plot_cell_trajectory(
    monocds,
    cell_size = 0.2, color_by = "acr20",
    show_tree = T, show_branch_points = F
) +
    facet_wrap(~subtype) +
    scale_color_manual(values = c("#87B1C8", "#A1568E")) +
    theme(
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.line.x = NULL,
        axis.line.y = NULL
    )
ggsave("Fig3/spp1_tof_mono_acr20.tiff", p, device = "tiff", width = 4.3, height = 2.7) # nolint

ss <- cbind(monocds@phenoData@data, t(monocds@reducedDimS))
ss <- subset(ss, subtype == "RA_AT")
colnames(ss)[32:33] <- c("com1", "com2")
ss$acr20 <- ifelse(ss$sampleid == "RA_5_1", "N", "Y")

p <- ggplot(ss, aes(x = com1, y = com2, color = sampleid)) +
    geom_point()

monocds@phenoData@data$group <- factor(monocds@phenoData@data$group, levels = c("RA_BT", "RA_AT")) # nolint

p <- plot_cell_trajectory(
    monocds,
    cell_size = 0.2, color_by = "acr20",
    show_tree = T, show_branch_points = F
) +
    facet_wrap(~subtype) +
    scale_color_manual(values = group_color[2:3]) +
    theme(
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.line.x = NULL,
        axis.line.y = NULL,
        legend.position = "none"
    )
ggsave("Fig3/spp1_tof_mono_group.tiff", p, device = "tiff", width = 2.6, height = 2.4) # nolint

p <- plot_cell_trajectory(
    monocds,
    cell_size = 0.2, color_by = "Pseudotime",
    show_tree = T, show_branch_points = FALSE,
    use_color_gradient = F
) +
    scale_color_gradientn(colors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)) + # nolint
    theme(
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.line.x = NULL,
        axis.line.y = NULL
    )
ggsave("Fig3/spp1_tof_mono_Pseudotime.tif", p, device = "tiff", width = 2, height = 2) # nolint
ggsave("Fig3/spp1_tof_mono_Pseudotime.pdf", p, width = 2, height = 2)

mq <- table(monocds@phenoData@data[, c("group", "newstate")]) / as.vector(table(monocds@phenoData@data$group)) # nolint
mq <- melt(mq)

p <- ggplot(mq, aes(x = newstate, y = value, fill = group)) +
    geom_bar(stat = "identity", position = "fill", width = 0.85) +
    scale_fill_manual(values = group_color[2:3]) +
    theme_classic() +
    labs(x = "", y = "Percentage") +
    scale_y_continuous(expand = c(0, 0), label = percent) +
    theme(
        panel.grid = element_blank(),
        axis.text = element_text(color = "black", size = 8),
        legend.position = "bottom",
        legend.title = element_blank()
    )

ggsave("Fig3/spp1_tof_per.pdf", p, width = 1.9, height = 3)

## density
plotdf <- pData(monocds)
plotdf1 <- subset(plotdf, newstate %in% c("State1", "State2"))
plotdf1$group <- factor(plotdf1$group, levels = rev(levels(plotdf1$group)))

library(ggridges)
p <- ggplot(
    plotdf1,
    aes(x = Pseudotime, y = group, fill = group)
) +
    geom_density_ridges(scale = 1) +
    scale_fill_manual(values = group_color[3:2]) +
    # geom_vline(xintercept = intercept, linetype = 2) +
    scale_y_discrete("") +
    theme_minimal() +
    theme(panel.grid = element_blank())

plotdf2 <- subset(plotdf, newstate %in% c("State1", "State3"))
plotdf2$group <- factor(plotdf2$group, levels = rev(levels(plotdf2$group)))
plotdf2$Pseudotime <- -1 * plotdf2$Pseudotime

p1 <- ggplot(
    plotdf2,
    aes(x = Pseudotime, y = group, fill = group)
) +
    geom_density_ridges(scale = 1) +
    scale_fill_manual(values = group_color[3:2]) +
    # geom_vline(xintercept = intercept, linetype = 2) +
    scale_y_discrete("") +
    theme_minimal() +
    theme(panel.grid = element_blank())
final <- p + p1 + plot_layout(nrow = 1, guides = "collect")
ggsave("Fig3/density.pdf", final, width = 6, height = 2)

beam_res <- BEAM(monocds, branch_point = 1, cores = 40)
beam_res <- beam_res[order(beam_res$qval), ]
beam_res <- beam_res[, c("gene_short_name", "pval", "qval")]

p <- plot_genes_branched_heatmap(
    monocds[row.names(subset(
        beam_res,
        qval < 1e-4
    )), ],
    branch_point = 1,
    num_clusters = 4,
    cores = 40,
    use_gene_short_name = T,
    show_rownames = F,
    return_heatmap = T
)
ggsave("heatmap.pdf", p$ph_res)

gene_group <- p$annotation_row
gene_group$gene <- rownames(gene_group)
write.table(gene_group, "mono_label.txt", quote = F, sep = "\t") # no

gene_group <- read.table("mono_label.txt", header = T, sep = "\t")

go <- function(x) {
    small_gene_group <- filter(gene_group, gene_group$Cluster == x)
    df_name <- bitr(
        small_gene_group$gene,
        fromType = "SYMBOL",
        toType = c("ENTREZID"), OrgDb = "org.Hs.eg.db"
    )
    go <- enrichGO(
        gene = unique(df_name$ENTREZID),
        OrgDb = org.Hs.eg.db,
        keyType = "ENTREZID",
        ont = "ALL",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2
    )
    go_res <- go@result
    return(go_res)
}

go1 <- go(1)
go4 <- go(4)
go3 <- go(3)
iterm <- c("T cell activation", "positive regulation of cytokine production", "lymphocyte proliferation") # nolint
iterm4 <- c(
    "macrophage migration", "fibroblast proliferation", "blood vessel endothelial cell migration", # nolint
    "regulation of T cell receptor signaling pathway", "cellular response to tumor necrosis factor" # nolint
) # nolint
iterm3 <- c(
    "response to interferon-gamma", "type I interferon production",
    "positive regulation of interleukin-1 production", "interleukin-12 production", # nolint
    "positive regulation of I-kappaB kinase/NF-kappaB signaling",
    "cell-cell adhesion", "toll-like receptor signaling pathway",
    "positive regulation of T cell mediated cytotoxicity" # nolint
)

regene <- function(x = go1[which(go1$Description == iterm[1]), "geneID"]) {
    entrezID <- as.vector(str_split(x, "/", simplify = T))
    SYMBOL <- bitr(entrezID, "ENTREZID", "SYMBOL", org.Hs.eg.db)[, "SYMBOL"]
    return(SYMBOL)
}
gene1 <- Reduce(union, list(
    regene(go1[which(go1$Description == iterm[1]), "geneID"]),
    regene(go1[which(go1$Description == iterm[2]), "geneID"]),
    regene(go1[which(go1$Description == iterm[3]), "geneID"])
))
gene4 <- Reduce(union, list(
    regene(go4[which(go4$Description == iterm4[1]), "geneID"]),
    regene(go4[which(go4$Description == iterm4[2]), "geneID"]),
    regene(go4[which(go4$Description == iterm4[3]), "geneID"]),
    regene(go4[which(go4$Description == iterm4[4]), "geneID"]),
    regene(go4[which(go4$Description == iterm4[5]), "geneID"])
))
gene3 <- Reduce(union, list(
    regene(go3[which(go3$Description == iterm3[1]), "geneID"]),
    regene(go3[which(go3$Description == iterm3[2]), "geneID"]),
    regene(go3[which(go3$Description == iterm3[3]), "geneID"]),
    regene(go3[which(go3$Description == iterm3[4]), "geneID"]),
    regene(go3[which(go3$Description == iterm3[5]), "geneID"]),
    regene(go3[which(go3$Description == iterm3[6]), "geneID"]),
    regene(go3[which(go3$Description == iterm3[7]), "geneID"]),
    regene(go3[which(go3$Description == iterm3[8]), "geneID"])
))

drr <- function(x) {
    p <- plot_genes_branched_heatmap(monocds[x, ],
        branch_point = 1,
        num_clusters = 1,
        cores = 40,
        use_gene_short_name = F,
        show_rownames = F,
        return_heatmap = T,
        hmcols = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(62)
    )
    return(p)
}
dr1 <- drr(gene1)
dr4 <- drr(gene4)
dr3 <- drr(gene3)
ggsave("Tactiviton.pdf", dr1$ph_res, width = 2.5, height = 1)
ggsave("macrophage_migration.pdf", dr4$ph_res, width = 2.5, height = 1)
ggsave("interferon.pdf", dr3$ph_res, width = 2.5, height = 1)

## Adalimumab
monocds <- readRDS("monocle/spp1_a.Rds")

monocds@phenoData@data$acr20 <- ifelse(monocds@phenoData@data$sampleid %in% c("RA_3_1", "RA_3_2"), "N", "Y") # nolint
monocds@phenoData@data$subtype <- factor(monocds@phenoData@data$subtype, levels = c("RA_BT", "RA_AT")) # nolint
p <- plot_cell_trajectory(
    monocds,
    cell_size = 0.2, color_by = "acr20",
    show_tree = T, show_branch_points = F
) +
    facet_wrap(~subtype) +
    scale_color_manual(values = c("#87B1C8", "#A1568E")) +
    theme(
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.line.x = NULL,
        axis.line.y = NULL
    )
ggsave("Fig3/spp1_tof_mono_acr20.tiff", p, device = "tiff", width = 4.3, height = 2.7) # nolint

p <- plot_cell_trajectory(
    monocds,
    cell_size = 0.2, color_by = "State",
    show_tree = T, show_branch_points = F
) +
    scale_color_manual(values = pal_npg()(9)) +
    theme(
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.line.x = NULL,
        axis.line.y = NULL,
        legend.position = "none"
    )
ggsave("Fig3/spp1_a_mono_state.tiff", p, device = "tiff", width = 2.6, height = 2.4) # nolint

monocds@phenoData@data$group <- factor(monocds@phenoData@data$group, levels = c("RA_BT", "RA_AT")) # nolint
p <- plot_cell_trajectory(
    monocds,
    cell_size = 0.2, color_by = "group",
    show_tree = T, show_branch_points = F
) +
    scale_color_manual(values = group_color[2:3]) +
    theme(
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.line.x = NULL,
        axis.line.y = NULL,
        legend.position = "none"
    )
ggsave("Fig3/spp1_a_mono_group.tiff", p, device = "tiff", width = 2.6, height = 2.4) # nolint

p <- plot_cell_trajectory(
    monocds,
    cell_size = 0.2, color_by = "Pseudotime",
    show_tree = T, show_branch_points = FALSE,
    use_color_gradient = F
) +
    scale_color_gradientn(colors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)) + # nolint
    theme(
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.line.x = NULL,
        axis.line.y = NULL,
        legend.position = "none"
    )
ggsave("Fig3/spp1_a_mono_Pseudotime.tif", p, device = "tiff", width = 2, height = 2) # nolint
ggsave("Fig3/spp1_tof_mono_Pseudotime.pdf", p, width = 2, height = 2)

monocds@phenoData@data$newstate <- ifelse(monocds@phenoData@data$State == 1, "State2", # nolint
    ifelse(monocds@phenoData@data$State == 2, "State3",
        ifelse(monocds@phenoData@data$State == 4, "State6",
            ifelse(monocds@phenoData@data$State == 5, "State5",
                ifelse(monocds@phenoData@data$State == 6, "State4", "State1")
            )
        )
    )
)
mq <- table(monocds@phenoData@data[, c("group", "newstate")]) / as.vector(table(monocds@phenoData@data$group)) # nolint
mq <- melt(mq)

p <- ggplot(mq, aes(x = newstate, y = value, fill = group)) +
    geom_bar(stat = "identity", position = "fill", width = 0.85) +
    scale_fill_manual(values = group_color[2:3]) +
    theme_classic() +
    labs(x = "", y = "Percentage") +
    scale_y_continuous(expand = c(0, 0), label = percent) +
    theme(
        panel.grid = element_blank(),
        axis.text = element_text(color = "black", size = 8),
        legend.position = "bottom",
        legend.title = element_blank()
    )

ggsave("Fig3/spp1_a_per.pdf", p, width = 2.7, height = 3)

## density
plotdf <- pData(monocds)
plotdf$group <- factor(plotdf$group, levels = rev(levels(plotdf$group)))

library(ggridges)
p <- ggplot(
    plotdf,
    aes(x = Pseudotime, y = group, fill = group)
) +
    geom_density_ridges(scale = 1) +
    scale_fill_manual(values = group_color[3:2]) +
    # geom_vline(xintercept = intercept, linetype = 2) +
    scale_y_discrete("") +
    theme_minimal() +
    theme(panel.grid = element_blank())

ggsave("Fig3/a_density.pdf", p, width = 3, height = 2)

my_pseudotime_de <- differentialGeneTest(
    monocds,
    fullModelFormulaStr = "~sm.ns(Pseudotime)", cores = 40
)
sig_gene_names <- row.names(subset(my_pseudotime_de, qval < 10^-20))
p <- plot_pseudotime_heatmap(
    monocds[sig_gene_names, ],
    num_clusters = 4, cores = 40,
    use_gene_short_name = TRUE, show_rownames = FALSE,
    return_heatmap = TRUE,
    hmcols = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(62)
)
ggsave("Fig3/a_mono_sig_heatmap.pdf", p, width = 6, height = 8) # nolint

gene_group <- cutree(p$tree_row, k = 4)
gene_group <- as.data.frame(gene_group)
gene_group$gene <- rownames(gene_group)
colnames(gene_group)[1] <- "Cluster"
write.table(gene_group, "Fig3/a_sig_label.txt", quote = F, sep = "\t") # nolint

gene_group <- read.table("Fig3/a_sig_label.txt", header = T, sep = "\t")

go <- function(x) {
    small_gene_group <- filter(gene_group, gene_group$Cluster == x)
    df_name <- bitr(
        small_gene_group$gene,
        fromType = "SYMBOL",
        toType = c("ENTREZID"), OrgDb = "org.Hs.eg.db"
    )
    go <- enrichGO(
        gene = unique(df_name$ENTREZID),
        OrgDb = org.Hs.eg.db,
        keyType = "ENTREZID",
        ont = "ALL",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2
    )
    go_res <- go@result
    return(go_res)
}

go1 <- go(1)
go4 <- go(4)
go3 <- go(3)
#############################################################


#################### SCENIC analysis ########################
scenic <- read.table("input/macrophage_scenic.txt", header = T, sep = "\t")
tr <- c(
    "STAT1 (222g)", "SPI1 (377g)", "KLF2 (80g)", "CEBPD_extended (15g)",
    "EGR1 (19g)", "FOSB (12g)", "NFKB1 (181g)", "BACH1_extended (22g)",
    "HIF1A (13g)", "CEBPB (310g)", "JUNB (63g)"
)

scenic <- scenic[tr, ]
mat <- data@meta.data[colnames(scenic), ]
mat <- subset(mat, celltype %in% c("SPP1+ Mac", "S100A12+ Mac"))
scenic <- scenic[, rownames(mat)]
mat$ggr <- paste(mat$subtype, mat$drug, mat$celltype, sep = "_")

mt <- matrix(NA, nrow = 11, ncol = 10)
for (i in 1:nrow(mt)) {
    mt[i, ] <- aggregate(as.vector(t(scenic)[, i]), list(mat$ggr), median)$x
}
colnames(mt) <- aggregate(as.vector(t(scenic)[, 1]), list(mat$ggr), median)$Group.1 # nolint
rownames(mt) <- rownames(scenic)
mat$tt <- paste(mat$subtype, mat$celltype, sep = "_")
mtt <- matrix(NA, nrow = 11, ncol = 6)
for (i in 1:nrow(mtt)) {
    mtt[i, ] <- aggregate(as.vector(t(scenic)[, i]), list(mat$tt), median)$x
}
colnames(mtt) <- aggregate(as.vector(t(scenic)[, 1]), list(mat$tt), median)$Group.1 # nolint
rownames(mtt) <- rownames(scenic)

mk <- data.frame(
    tr = rownames(mtt), spp1_a = mtt[, 6] / mtt[, 2],
    spp1_b = mt[, 8] / mt[, 4], spp1_c = mt[, 10] / mt[, 6],
    s100a_a = mtt[, 5] / mtt[, 1], s100a_b = mt[, 7] / mt[, 3],
    s100a_c = mt[, 9] / mt[, 5]
)
for (i in 1:nrow(mk)) {
    for (j in 2:ncol(mk)) {
        mk[i, j] <- log(mk[i, j], 2)
    }
}
for (i in 1:nrow(mk)) {
    for (j in 2:ncol(mk)) {
        if (mk[i, j] %in% c("-Inf", "Inf", "NaN")) {
            mk[i, j] <- 0
        }
    }
}
ma <- matrix(NA, nrow = nrow(mk), ncol = 6)
for (i in 1:nrow(mk)) {
    ma[i, 1] <- wilcox.test(as.numeric(scenic[i, which(mat$tt == "RA_BT_SPP1+ Mac" & mat$celltype == "SPP1+ Mac")]), as.numeric(scenic[i, which(mat$tt == "OA_BT_SPP1+ Mac" & mat$celltype == "SPP1+ Mac")]))$p.value # nolint
    ma[i, 2] <- wilcox.test(as.numeric(scenic[i, which(mat$ggr == "RA_BT_Adalimumab_SPP1+ Mac" & mat$celltype == "SPP1+ Mac")]), as.numeric(scenic[i, which(mat$ggr == "RA_AT_Adalimumab_SPP1+ Mac" & mat$celltype == "SPP1+ Mac")]))$p.value # nolint
    ma[i, 3] <- wilcox.test(as.numeric(scenic[i, which(mat$ggr == "RA_BT_Tofacitinib_SPP1+ Mac" & mat$celltype == "SPP1+ Mac")]), as.numeric(scenic[i, which(mat$ggr == "RA_AT_Tofacitinib_SPP1+ Mac" & mat$celltype == "SPP1+ Mac")]))$p.value # nolint
    ma[i, 4] <- wilcox.test(as.numeric(scenic[i, which(mat$tt == "RA_BT_S100A12+ Mac" & mat$celltype == "S100A12+ Mac")]), as.numeric(scenic[i, which(mat$tt == "OA_BT_S100A12+ Mac" & mat$celltype == "S100A12+ Mac")]))$p.value # nolint
    ma[i, 5] <- wilcox.test(as.numeric(scenic[i, which(mat$ggr == "RA_BT_Adalimumab_S100A12+ Mac" & mat$celltype == "S100A12+ Mac")]), as.numeric(scenic[i, which(mat$ggr == "RA_AT_Adalimumab_S100A12+ Mac" & mat$celltype == "S100A12+ Mac")]))$p.value # nolint
    ma[i, 6] <- wilcox.test(as.numeric(scenic[i, which(mat$ggr == "RA_BT_Tofacitinib_S100A12+ Mac" & mat$celltype == "S100A12+ Mac")]), as.numeric(scenic[i, which(mat$ggr == "RA_AT_Tofacitinib_S100A12+ Mac" & mat$celltype == "S100A12+ Mac")]))$p.value # nolint
}
colnames(ma) <- colnames(mk)[2:7] # nolint
rownames(ma) <- rownames(scenic)
mk <- melt(mk)
mk$pvalue <- melt(ma)$value
mk$value <- ifelse(mk$value > 1, 1, ifelse(mk$value < -1, -1, mk$value))

p <- ggplot(mk, aes(x = variable, y = tr)) +
    geom_point(aes(size = -log(pvalue, 10), color = value)) +
    scale_color_gradientn(colors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)) + # nolint
    geom_vline(xintercept = 3.5, color = "gray", linetype = "dashed") +
    theme_bw() +
    labs(x = "", y = "") +
    scale_x_discrete(labels = rep(c("RA-BT vs. OA-BT", "RA-BT vs. RA-AT\n(Adalimumab)", "RA-BT vs. RA-AT\n(Tofacitinib)"), 2)) + # nolint
    theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, color = "black", hjust = 1),
        axis.text.y = element_text(color = "black", size = 12)
    )
p
ggsave("Fig3/tf.pdf", p, width = 5.9, height = 4.5)
#############################################################
