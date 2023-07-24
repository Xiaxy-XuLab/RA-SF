#!/usr/bin/env R
# coding utf-8

###################### Package Loading ###########################
pacman::p_load(
    "Seurat", "ggplot2", "ggsci", "reshape2", "scales",
    "monocle", "viridis", "patchwork", "ggpubr", "ggrepel",
    "msigdbr", "fgsea", "tidyverse", "bseqsc", "ggsignif",
    "ggrepel", "pheatmap", "RColorBrewer", "ggpmisc", "ggplotify",
    "circlize", "zoo", "nichenetr", "tidyverse", "GSVA", "magrittr",
    "ggforce", "dplyr"
)
setwd("/work/xiaxy/work/RA/")
source("/work/xiaxy/work/RA/Dependent.R")
#############################################################


###################### Read input ###########################
data_input <- readRDS("data_rds/input.Rds")
#############################################################


####################### Tsne Dotplot ########################
## tsne
p <- DimPlot(
    data_input,
    label = FALSE, reduction = "tsne"
) +
    scale_color_manual(values = celltype_color) + # nolint
    theme(
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank()
    )
ggsave("Fig1/a_umap.pdf", p, width = 3, height = 5)

p <- DimPlot(data_input,
    label = FALSE, group.by = "subtype",
    reduction = "tsne", shuffle = T
) +
    scale_color_manual(values = group_color3) + # nolint
    theme(
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank()
    )
ggsave("Fig1/b_umap.pdf", p, width = 3, height = 3)

## dotplot
pos_calculate <- function(x) {
    return(100 * sum(x > 0) / length(x))
}

genes <- c(
    "CD68", "VCAN", "CD14", "CD1C", "FCER1A", "CLEC10A",
    "CD3D", "CD3E", "CD2", "NKG7", "KLRD1", "PRF1", "CD79A",
    "MS4A1", "CD19", "JCHAIN", "TCF4", "IGKC", "CTSK",
    "MMP9", "ACP5", "COL1A2", "DCN", "PRG4"
)
exp_pos <- matrix(NA, nrow = length(genes), ncol = length(unique(data_input$celltype))) # nolint
exp_mean <- matrix(NA, nrow = length(genes), ncol = length(unique(data_input$celltype))) # nolint

for (i in 1:length(genes)) {
    exp_mean[i, ] <- aggregate(data_input@assays$RNA@data[genes[i], ], list(data_input$celltype), mean)$x # nolint
    exp_pos[i, ] <- aggregate(data_input@assays$RNA@data[genes[i], ], list(data_input$celltype), pos_calculate)$x # nolint
}

mat <- data.frame(
    gene = factor(rep(genes, times = length(unique(data_input$celltype))), levels = rev(genes)), # nolint
    celltype = rep(aggregate(data_input@assays$RNA@data[genes[1], ], list(data_input$celltype), pos_calculate)$Group.1, each = length(genes)), # nolint
    exp_pos = melt(exp_pos)$value,
    exp_mean = melt(exp_mean)$value
)
mat$exp_pos <- ifelse(mat$exp_pos < 1, NA, mat$exp_pos)

p <- ggplot(mat) +
    geom_point(aes(
        x = celltype, y = gene, fill = exp_mean,
        size = exp_pos
    ), color = "black", shape = 21) +
    scale_fill_gradientn(colors = dotplot_color) +
    theme_bw() +
    scale_y_discrete(position = "right") +
    scale_x_discrete(position = "bottom") +
    labs(fill = "Mean expression", size = "Percent expressed") +
    geom_hline(yintercept = c(3.5, 6.5, 9.5, 12.5, 15.5, 18.5, 21.5), color = "lightgrey", linetype = "dashed") + # nolint
    theme(
        axis.title = element_blank(),
        legend.position = "bottom",
        axis.text.x = element_text(angle = 270, color = "black", size = 10, hjust = 0, vjust = 0.5), # nolint
        axis.text.y = element_text(size = 10, color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(linewidth = 1.2, color = "black")
    ) +
    guides(
        fill = guide_legend(override.aes = list(size = 2), title.position = "top", title.hjust = 0.5, nrow = 1, keywidth = 0.1), # nolint
        size = guide_legend(title.position = "top", title.hjust = 0.5, nrow = 1, keywidth = 0.1) # nolint
    )
ggsave("Fig1/c_dotplot.pdf", p, width = 3.8, height = 6.5)
#############################################################


################### Percentage analysis #####################
## p1
mat <- melt(table(data_input@meta.data[, c("sampleid", "celltype")]) / as.vector(table(data_input$sampleid))) # nolint
mat$patient <- paste(
    sapply(strsplit(as.character(mat$sampleid), "_", fixed = T), "[", 1),
    sapply(strsplit(as.character(mat$sampleid), "_", fixed = T), "[", 2),
    sep = "_"
)
mat$value <- ifelse(
    sapply(strsplit(as.character(mat$sampleid), "_", fixed = T), "[", 3) == 2,
    mat$value, -1 * mat$value
)
mat$patient <- factor(mat$patient, levels = rev(unique(mat$patient)))
mat$celltype <- factor(mat$celltype, levels = rev(levels(mat$celltype)))

p <- ggplot(mat, aes(x = value, y = patient, fill = celltype)) +
    geom_bar(position = "fill", stat = "identity", width = 0.85) +
    scale_fill_manual(values = rev(celltype_color)) + # nolint
    theme_bw() +
    guides(fill = guide_legend(title = NULL)) +
    labs(y = "", x = "Proportion of clusters") +
    theme(
        legend.title = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.text.y = element_text(color = "black", size = 10), # nolint
        axis.text.x = element_text(color = "black", size = 10),
        axis.title = element_text(color = "black", size = 12)
    ) +
    scale_x_continuous(
        expand = c(0, 0.01),
        labels = percent, position = "top"
    ) +
    guides(fill = guide_legend(reverse = TRUE))
ggsave("Fig1/d_percent.pdf", p, width = 4.3, height = 3)

## p2
mat <- melt(table(data_input@meta.data[, c("subtype", "celltype")]) / as.vector(table(data_input$subtype))) # nolint
mat$subtype <- factor(mat$subtype, levels = c("OA_BT", "RA_BT", "RA_AT"))

p <- ggplot(mat, aes(x = celltype, y = value, fill = subtype)) +
    geom_bar(position = "fill", stat = "identity") +
    scale_fill_manual(values = group_color3) +
    theme_bw() +
    labs(x = "", y = "Proportion of groups") +
    theme(
        legend.title = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 10), # nolint
        axis.text.y = element_text(color = "black", size = 10),
        axis.title = element_text(color = "black", size = 12)
    ) +
    scale_y_continuous(expand = c(0, 0.01), labels = percent)
ggsave("Fig1/e_percent.pdf", p, width = 4.3, height = 2.7)

## p3
mat <- melt(table(data_input@meta.data[, c("sampleid", "celltype")]) / as.vector(table(data_input$sampleid))) # nolint
mat %<>% mutate(., patient = paste(
    sapply(strsplit(as.character(mat$sampleid), "_", fixed = T), "[", 1),
    sapply(strsplit(as.character(mat$sampleid), "_", fixed = T), "[", 2),
    sep = "_"
)) %<>% mutate(.,
    group = factor(
        rep(
            c(rep("OA_BT", 3), rep(c("RA_BT", "RA_AT"), 6)),
            length(unique(mat$celltype))
        ),
        levels = c("OA_BT", "RA_BT", "RA_AT")
    )
) %<>% mutate(.,
    drug = factor(
        rep(c(rep("None", 3), rep(c("Adalimumab", "Tofacitinib"), each = 6)), length(unique(mat$celltype))), # nolint
        levels = c("None", "Adalimumab", "Tofacitinib")
    )
) %<>% mutate(.,
    value = 100 * value
)

draw <- function(x) {
    data_set <- subset(mat, celltype == x)
    p <- ggplot(data_set, aes(group, value)) +
        stat_summary(fun.y = mean, geom = "bar", color = "black", fill = "white", width = .6) + # nolint
        stat_summary(fun.data = mean_se, geom = "errorbar", color = "black", width = .2) + # nolint
        geom_point(aes(color = group, shape = drug), size = 5) +
        geom_line(aes(group = patient), linewidth = 0.6, colour = "#9C9C9C") +
        scale_shape_manual(values = c(15, 1, 19)) +
        scale_color_manual(values = group_color3) +
        theme_classic() +
        stat_compare_means(comparisons = list(c("OA_BT", "RA_BT"), c("OA_BT", "RA_AT"), c("RA_BT", "RA_AT"))) + # nolint
        scale_y_continuous(expand = c(0, 0), limits = c(0, 50, 100), breaks = c(0, 50, 100), labels = c("0%", "50%", "100%")) + # nolint
        labs(y = "Proportion", title = x) +
        theme(
            axis.title.x = element_blank(),
            plot.title = element_text(hjust = 0.5, color = "black", size = 16),
            axis.title.y = element_text(size = 16, color = "black"),
            axis.text.x = element_text(size = 14, color = "black", angle = 45, hjust = 1), # nolint
            axis.text.y = element_text(size = 14, color = "black"),
            axis.ticks.length = unit(0.4, "lines"),
            legend.position = "top",
            legend.title = element_blank(),
            legend.text = element_text(size = 14, color = "black")
        ) +
        guides(
            shape = guide_legend(order = 1),
            color = "none"
        )
    return(p)
}
final <- draw("Macrophage") + draw("DC") + draw("T cell") +
    plot_layout(guides = "collect") &
    theme(legend.position = "top")
ggsave("Fig1/f_percent.pdf", final, width = 7.2, height = 5)
#############################################################


################# Clinical index analysis ###################
## p1
clinical_index <- read.table("input/validation_cohort_clinical_index.txt", header = T, sep = "\t") # nolint
clinical_index %<>% mutate(., Group = factor(Group, levels = c("RA-BT", "RA-AT"))) # nolint

clinical_graph <- function(x) {
    p <- ggplot(clinical_index, aes(Group, get(x))) +
        stat_summary(fun.y = mean, geom = "bar", color = "black", fill = "white", width = .6) + # nolint
        stat_summary(fun.data = mean_se, geom = "errorbar", color = "black", width = .2) + # nolint
        geom_point(aes(color = Group, shape = Drug), size = 5) +
        scale_shape_manual(values = c(1, 15, 19)) +
        scale_color_manual(values = group_color3[2:3]) +
        theme_classic() +
        geom_line(aes(group = PatientID), size = 0.6, colour = "#9C9C9C") +
        stat_compare_means(comparisons = list(c("RA-BT", "RA-AT")), label.y = 1.05 * max(clinical_index[, x], na.rm = T)) + # nolint
        scale_y_continuous(expand = c(0, 0), limits = c(0, floor(max(clinical_index[, x], na.rm = T) * 1.2))) + # nolint
        labs(y = x) +
        theme(
            axis.title.x = element_blank(),
            plot.title = element_text(hjust = 0.5, color = "black", size = 16),
            axis.title.y = element_text(size = 16, color = "black"),
            axis.text.x = element_text(size = 14, color = "black", angle = 45, hjust = 1), # nolint
            axis.text.y = element_text(size = 14, color = "black"),
            axis.ticks.length = unit(0.4, "lines"),
            legend.position = "top",
            legend.title = element_blank(),
            legend.text = element_text(size = 14, color = "black")
        ) +
        guides(
            shape = guide_legend(order = 1),
            color = "none"
        )
    return(p)
}

final <- clinical_graph("DAS28") + clinical_graph("SDAI") +
    plot_layout(nrow = 1, guides = "collect")
ggsave("Fig1/vali_info.pdf", final, width = 7, height = 4)

## p2
p <- ggplot(clinical_index_single, aes(x = DAS28, y = SDAI)) +
    geom_point(color = "#C76862", size = 3) +
    geom_smooth(aes(x = DAS28, y = SDAI), se = FALSE, method = "lm", color = "black") + # nolint
    stat_cor(method = "pearson") +
    labs(y = "SDAI", x = "DAS28") +
    theme_classic() +
    theme(
        panel.grid = element_blank(),
        axis.text = element_text(size = 8, color = "black")
    )
ggsave("Fig1/das28_SDAI.pdf", p, width = 3, height = 3)

## p3
clinical_index_ACR20 <- read.table("input/validation_cohort_clinical_index_ACR20.txt", header = T, sep = "\t") # nolint
clinical_graph_ACR20 <- function(x) {
    p <- ggplot(clinical_index_ACR20, aes(ACR20, get(x))) +
        stat_summary(fun.y = mean, geom = "bar", color = "black", fill = "white", width = .6) + # nolint
        stat_summary(fun.data = mean_se, geom = "errorbar", color = "black", width = .2) + # nolint
        geom_point(aes(color = ACR20, shape = Drug), size = 5) +
        scale_shape_manual(values = c(1, 19)) +
        scale_color_manual(values = c("#87B1C8", "#A1568E")) +
        theme_classic() +
        stat_compare_means(comparisons = list(c("N", "Y")), label.y = 1.05 * max(clinical_index_ACR20[, x], na.rm = T)) + # nolint
        scale_y_continuous(expand = c(0, 0), limits = c(0, round(max(clinical_index_ACR20[, x], na.rm = T) * 1.2, digits = 2))) + # nolint
        labs(y = x) +
        theme(
            axis.title.x = element_blank(),
            plot.title = element_text(hjust = 0.5, color = "black", size = 16),
            axis.title.y = element_text(size = 16, color = "black"),
            axis.text.x = element_text(size = 14, color = "black", angle = 45, hjust = 1), # nolint
            axis.text.y = element_text(size = 14, color = "black"),
            axis.ticks.length = unit(0.4, "lines"),
            legend.position = "top",
            legend.title = element_blank(),
            legend.text = element_text(size = 14, color = "black")
        ) +
        guides(
            shape = guide_legend(order = 1),
            color = "none"
        )
    return(p)
}
final <- clinical_graph_ACR20("DAS28") + clinical_graph_ACR20("SDAI") +
    plot_layout(nrow = 1, guides = "collect")
ggsave("Fig1/vali_ACR20.pdf", final, width = 6.5, height = 4)
#############################################################


################# Flow cytometry analysis ###################
flow_cytometry <- read.table("input/Flow_cytometry_CD64_CD11b.txt", header = T, sep = "\t") # nolint
p <- ggplot(flow_cytometry, aes(x = group, y = value)) +
    geom_boxplot() +
    geom_point(aes(fill = group),
        position = position_auto(scale = FALSE, jitter.width = 0.4),
        color = "black", shape = 21, size = 3
    ) +
    labs(y = "% CD11b+CD64+\namong live cells", x = "") +
    scale_fill_manual(values = group_color3[1:2]) +
    theme_classic() +
    scale_y_continuous(limits = c(0, 100)) +
    stat_compare_means(comparisons = list(c("OA", "RA")), label.y = 90) +
    theme(
        legend.position = "none",
        axis.text = element_text(color = "black")
    )
ggsave("Fig1/flow_cytometry.pdf", p, width = 3.5, height = 4)
#############################################################


################## Deconvolution analysis ###################
## p1
bulk <- read.table("input/Bulk_abundance.txt", header = T, sep = "\t", row.names = 1) # nolint
bseqsc_config("/work/xiaxy/CIBERSORT.R")

count <- as.matrix(data_input@assays$RNA@counts)
celltype <- data.frame(
    cell = rownames(data_input@meta.data),
    celltype = as.character(data_input$celltype)
)
markerlist <- list(
    Macrophage = c("CD68", "FCGR1A", "CD163", "VCAN", "CD80", "CD86"),
    DC = c("FCER1A", "CD1C", "CD1E", "PLD4"),
    "T cell" = c("CD3D", "LTB", "IL7R", "CD3G"),
    "NK cell" = c("NKG7", "KLRD1", "KLRC1", "KLRB1", "KLRF1"),
    "B cell" = c("CD79A", "MS4A1", "CD19"),
    "Plasma cell" = c("JCHAIN", "TCF4", "PLD4", "BCL11A"),
    Osteoclast = c("MMP9", "CTSK", "APOE", "MMP14"),
    Fibroblast = c("PRG4", "CLU", "MMP3", "LUM", "COL1A2", "COL3A1")
)

bseqsc_result <- bseqsc_basis(
    count, markerlist,
    clusters = celltype$celltype,
    samples = data_input$sampleid, ct.scale = TRUE
)
fit <- bseqsc_proportions(bulk, bseqsc_result, verbose = TRUE)

result_data <- melt(as.matrix(fit$coefficients)) %>%
    arrange(Var1) %>%
    mutate(
        group = factor(rep(c(rep("OA-BT", 5), rep("RA-BT", 14), rep("RA-AT", 10)), length(markerlist)), # nolint
            levels = c("OA-BT", "RA-BT", "RA-AT")
        )
    )

p <- ggplot(result_data, aes(x = Var1, y = value)) +
    geom_boxplot(aes(color = group), outlier.size = 0, outlier.colour = "white") + # nolint
    geom_point(aes(fill = group),
        position = position_jitterdodge(jitter.width = 0.4),
        color = "black", shape = 21
    ) +
    scale_color_manual(values = group_color3) +
    scale_fill_manual(values = group_color3) +
    theme_classic2() +
    labs(y = "Deconvolution score", x = "") +
    scale_y_continuous(
        expand = c(0, 0.01), limits = c(0, 1)
    ) +
    theme(
        axis.text.x = element_text(color = "black", size = 14, angle = 45, hjust = 1), # nolint
        axis.text.y = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 16),
        legend.position = c(0.15, 0.9), legend.title = element_blank()
    )
p
ggsave("Fig1/m_bulk-decon1.pdf", p, width = 5, height = 4)

## p2
gse <- read.table("input/GSE55235_Deconvolution_result.txt", header = T, sep = "\t", row.names = 1) # nolint

mat <- data.frame(
    macro = as.numeric(gse[1, ]),
    group = c(rep("OA", 10), rep("RA", 10))
)
mat$group <- factor(mat$group, levels = c("OA", "RA"))

p <- ggplot(mat, aes(x = group, y = macro)) +
    geom_boxplot(aes(color = group), outlier.size = 0, outlier.colour = "white") + # nolint
    geom_point(aes(fill = group),
        position = position_jitterdodge(jitter.width = 0.7),
        color = "black", shape = 21, size = 4
    ) +
    scale_color_manual(values = group_color3) +
    scale_fill_manual(values = group_color3) +
    stat_compare_means(comparisons = list(c("OA", "RA")), label.y = 0.9) +
    theme_classic2() +
    labs(y = "Deconvolution score", x = "", title = "GSE55235\nMacrophage") +
    scale_y_continuous(
        expand = c(0, 0.01), limits = c(0, 1)
    ) +
    theme(
        axis.text.x = element_text(color = "black", size = 14), # nolint
        axis.text.y = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 16),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5)
    )
p
ggsave("Fig1/n_bulk-decon1.pdf", p, width = 2.4, height = 3)
#############################################################
