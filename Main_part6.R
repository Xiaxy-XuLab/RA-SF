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
data_input <- readRDS("data_rds/input.Rds")
data <- readRDS("data_rds/tcell.Rds")
#############################################################


###################### DEG number ###########################
DEG <- function(x) {
    matr <- matrix(NA, nrow = 10, ncol = 2)
    k <- 1
    for (i in unique(x$celltype)) {
        print(i)
        sub <- subset(x, cells = rownames(x@meta.data)[which(x$celltype == i)]) # nolint
        matr[k, 1] <- i
        if (table(x@active.ident)[1] > 10 & table(x@active.ident)[2] > 10) {
            marker <- FindAllMarkers(sub, only.pos = T, min.pct = 0.1, logfc.threshold = 0.25) # nolint
            if (nrow(marker) > 0) {
                am <- subset(marker, p_val_adj < 0.05)
                if (nrow(am) > 0) {
                    matr[k, 2] <- nrow(am)
                } else {
                    matr[k, 2] <- NA
                }
            }
        }
        k <- k + 1
    }
    colnames(matr) <- c("celltype", "num")
    matr <- as.data.frame(matr)
    dmat <- as.data.frame(cbind(data@meta.data, data@reductions$tsne@cell.embeddings)) # nolint
    dmat_mer <- merge(dmat, matr, by.x = "celltype", by.y = "celltype")
    return(dmat_mer)
}

deg_graph <- function(x, y) {
    p1 <- ggplot(x, aes(x = tSNE_1, y = tSNE_2, color = as.numeric(num))) + # nolint
        geom_point(size = 0.1) +
        scale_color_gradientn(colors = rev(colorRampPalette(colors = deg_color)(100))) + # nolint
        theme_classic() +
        labs(color = "", title = "") +
        theme(
            axis.title.x = element_blank(),
            axis.title.y = element_blank(), # nolint
            axis.line = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            plot.margin = margin(0, 0, 0, 0, "cm"),
            plot.title = element_text(hjust = 0.5, vjust = 0, color = "black", size = 14), # nolint
            legend.position = y
        ) +
        guides(color = guide_colourbar(title.position = "top", title.hjust = 1)) # nolint
    return(p1)
}

sub1 <- subset(data, cells = rownames(data@meta.data)[which(data$subtype %in% c("OA_BT", "RA_BT"))]) # nolint
sub1@active.ident <- factor(sub1$subtype)
sub2 <- subset(data, cells = rownames(data@meta.data)[which(data$subtype %in% c("RA_BT", "RA_AT"))]) # nolint
sub2@active.ident <- factor(sub2$subtype)
ddr1 <- DEG(sub1)
ddr2 <- DEG(sub2)

sub4 <- subset(data, cells = rownames(data@meta.data)[
    which(data$subtype %in% c("RA_BT", "RA_AT") & data$sampleid %in% unique(data$sampleid)[4:9]) # nolint
])
sub4@active.ident <- factor(sub4$subtype)
sub5 <- subset(data, cells = rownames(data@meta.data)[
    which(data$subtype %in% c("RA_BT", "RA_AT") & data$sampleid %in% unique(data$sampleid)[10:15]) # nolint
])
sub5@active.ident <- factor(sub5$subtype)

ddr4 <- DEG(sub4)
ddr5 <- DEG(sub5)

final <- deg_graph(ddr1, "bottom") + deg_graph(ddr2, "bottom") +
    deg_graph(ddr4, "bottom") + deg_graph(ddr5, "bottom") +
    plot_layout(nrow = 2)
ggsave("Fig4/d_deg4.pdf", final, width = 8.2, height = 6)

## p2
## CXCL13+ CD4+ T
marker_a <- read.table("input/cxcl13_deg_OA_BT-RA_BT.txt", header = T, sep = "\t") # nolint
marker_b <- read.table("input/cxcl13_deg_RA_BT-RA_AT.txt", header = T, sep = "\t") # nolint
marker_c <- read.table("input/cxcl13_deg_RA_BT-RA_AT_ada.txt", header = T, sep = "\t") # nolint
marker_d <- read.table("input/cxcl13_deg_RA_BT-RA_AT_tof.txt", header = T, sep = "\t") # nolint

vol <- function(x, y, z) {
    x[, "avg_log2FC"] <- ifelse(x[, "cluster"] == y, -1 * x[, "avg_log2FC"], x[, "avg_log2FC"]) # nolint
    x[, "p_val_adj"] <- ifelse(x[, "p_val_adj"] == 0,
        min(x[which(x[, "p_val_adj"] != 0), "p_val_adj"]), x[, "p_val_adj"]
    )
    x[, "label"] <- ifelse(x[, "gene"] %in% z, x[, "gene"], NA)
    x[, "group"] <- ifelse(x[, "p_val_adj"] >= 0.05, "N.S.",
        ifelse(x[, "avg_log2FC"] > 0, "Up", "Down")
    )
    x[, "group"] <- factor(x[, "group"], levels = c("N.S.", "Down", "Up")) # nolint
    x <- x[order(x$group), ]

    p <- ggplot(x, aes(
        x = avg_log2FC, y = -log(p_val_adj, 10),
        color = group
    )) +
        geom_point(size = 1) +
        scale_color_manual(
            values = c(
                "lightgrey",
                "#77A270", "#E3B66B"
            )
        ) +
        # geom_text_repel(aes(label = label), color = "#3750A1", max.overlaps = 30, size = 2.5) + # nolint
        theme_classic() +
        labs(title = "", y = "", x = "") + # nolint
        theme(
            plot.title = element_text(hjust = 0.5, size = 12),
            legend.position = "none",
            axis.text = element_text(size = 10, color = "black"),
            axis.title = element_text(size = 12, color = "black")
        )
    return(p)
}

label1 <- c("LAG3", "ITGA4", "JUN", "JUND", "JUNB", "IRF1", "UBE2K", "SPP1", "STAT3", "STAT5B", "H3F3B", "JAK3", "HLA-DRB1", "XCL1", "FCER1A") # nolint
label2 <- c("JAK3", "IFI16", "IRF1", "STAT1", "STAT3", "JAK1", "STAT5B", "IL6R", "S100A10", "RGCC", "S100A6") # nolint
label3 <- c("S100A4", "S100A10", "S100A6", "STAT5B", "JAK3", "IRF1", "JUND", "IFI16", "CCL4", "CXCL8", "JUNB", "IL6R") # nolint
label4 <- c("SPP1", "CXCR6", "JAK3", "IL32", "STAT3", "IFI16", "STAT1", "IRF1", "IL16", "CCR5", "ISG20", "S100A10", "RGCC", "JUND") # nolint
final <- vol(marker_a, "OA_BT", label1) + vol(marker_b, "RA_AT", label2) +
    vol(marker_c, "RA_AT", label3) + vol(marker_d, "RA_AT", label4) +
    plot_layout(nrow = 1) &
    theme(
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank()
    ) &
    labs(title = "")
ggsave("Fig4/m_vol.tif", final, device = "tiff", width = 14, height = 3.5)
#############################################################


################## Hallmarker analysis ######################
## p1
esc <- read.table("input/tcell_hall.txt", header = T, sep = "\t")
mat <- aggregate(esc, list(data$celltype), mean)
rownames(mat) <- mat$Group.1
mat <- mat[, -1]
colnames(mat) <- gsub("HALLMARK_", "", colnames(mat))
mat <- as.data.frame(t(mat))
mat1 <- mat[c(4, 22:27, 35, 44, 45), ]

annotation_col <- data.frame(
    celltype = colnames(mat1)
)
rownames(annotation_col) <- colnames(mat1)
col1 <- celltype_color
names(col1) <- colnames(mat1)

ann_colors <- list(
    celltype = col1
)

bk <- c(seq(-1, -0.1, by = 0.01), seq(0, 1, by = 0.01))
pdf("Fig2/b_hall.pdf", width = 4.2, height = 3)
pheatmap(mat1,
    show_colnames = F, annotation_col = annotation_col,
    annotation_colors = ann_colors,
    show_rownames = T, scale = "row", cluster_cols = T,
    cluster_rows = T, color = c(
        colorRampPalette(colors = brewer.pal(11, "RdYlBu")[11:6])(length(bk) / 2), # nolint
        colorRampPalette(colors = brewer.pal(11, "RdYlBu")[6:1])(length(bk) / 2)
    ),
    breaks = bk, annotation_legend = TRUE,
    clustering_method = "ward.D",
    legend_breaks = c(-1, 0, 1),
    fontsize = 6, annotation_names_col = FALSE,
    border_color = "white", treeheight_row = 8, treeheight_col = 8
)
dev.off()

## p2
esc <- read.table("input/tcell_hall.txt", header = T, sep = "\t")
mat <- aggregate(esc, list(data$celltype), mean)
rownames(mat) <- mat$Group.1
mat <- mat[, -1]
colnames(mat) <- gsub("HALLMARK_", "", colnames(mat))
mat <- as.data.frame(t(mat))
mat1 <- mat[-c(4, 22:27, 35, 44, 45), ]

annotation_col <- data.frame(
    celltype = colnames(mat1)
)
rownames(annotation_col) <- colnames(mat1)
col1 <- celltype_color
names(col1) <- colnames(mat1)

ann_colors <- list(
    celltype = col1
)

bk <- c(seq(-1, -0.1, by = 0.01), seq(0, 1, by = 0.01))
pdf("Fig2/c_hall.pdf", width = 4.2, height = 3.8)
pheatmap(mat1,
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

## p3
esc <- read.table("input/tcell_hall.txt", header = T, sep = "\t")
colnames(esc) <- gsub("HALLMARK_", "", colnames(esc))
esc <- esc[rownames(data@meta.data)[which(data$celltype %in% c("GZMK+ CD8+ Tem1", "Treg", "CD4+ Tex"))], c(4, 22:27, 35, 44, 45)] # nolint
sub <- subset(data, cells = rownames(esc))

esc <- esc[, c("TNFA_SIGNALING_VIA_NFKB", "INTERFERON_GAMMA_RESPONSE", "IL6_JAK_STAT3_SIGNALING")] # nolint
esc$group <- paste(sub$subtype, sub$drug, sep = "_")
esc$group <- factor(esc$group, levels = unique(esc$group))
esc$celltype <- sub$celltype
my_comparisons <- list(
    c("OA_BT_None", "RA_BT_Adalimumab"), c("OA_BT_None", "RA_BT_Tofacitinib"),
    c("RA_BT_Adalimumab", "RA_AT_Adalimumab"), c("RA_BT_Tofacitinib", "RA_AT_Tofacitinib") # nolint
)

p1 <- ggplot(esc, aes(x = group, y = TNFA_SIGNALING_VIA_NFKB, fill = group)) +
    geom_violin() +
    geom_boxplot(width = 0.2, position = position_dodge(0.9), outlier.colour = NA, fill = "white") + # nolint
    theme_bw() +
    labs(title = "TNFA_SIGNALING_VIA_NFKB", y = "Pathway activity") +
    facet_wrap(~celltype) +
    theme(
        legend.position = "none",
        axis.text.x = element_text(hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(hjust = 0.5, vjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = rel(1)),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black", size = 1),
        plot.title = element_text(hjust = 0.5, color = "black", size = 12)
    ) +
    scale_fill_manual(values = group_color5) + # nolint
    stat_compare_means(comparisons = my_comparisons, label = "p.signif")

p2 <- ggplot(esc, aes(x = group, y = INTERFERON_GAMMA_RESPONSE, fill = group)) +
    geom_violin() +
    geom_boxplot(width = 0.2, position = position_dodge(0.9), outlier.colour = NA, fill = "white") + # nolint
    theme_bw() +
    labs(title = "INTERFERON_GAMMA_RESPONSE", y = "Pathway activity") +
    facet_wrap(~celltype) +
    theme(
        legend.position = "none",
        axis.text.x = element_text(hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(hjust = 0.5, vjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = rel(1)),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black", size = 1),
        plot.title = element_text(hjust = 0.5, color = "black", size = 12)
    ) +
    scale_fill_manual(values = group_color5) + # nolint
    stat_compare_means(comparisons = my_comparisons, label = "p.signif")

p3 <- ggplot(esc, aes(x = group, y = IL6_JAK_STAT3_SIGNALING, fill = group)) +
    geom_violin() +
    geom_boxplot(width = 0.2, position = position_dodge(0.9), outlier.colour = NA, fill = "white") + # nolint
    theme_bw() +
    labs(title = "IL6_JAK_STAT3_SIGNALING", y = "Pathway activity") +
    facet_wrap(~celltype) +
    theme(
        legend.position = "none",
        axis.text.x = element_text(hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(hjust = 0.5, vjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = rel(1)),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black", size = 1),
        plot.title = element_text(hjust = 0.5, color = "black", size = 12)
    ) +
    scale_fill_manual(values = group_color5) + # nolint
    stat_compare_means(comparisons = my_comparisons, label = "p.signif")

final <- p1 + p2 + p3 + plot_layout(nrow = 1, guides = "collect")
ggsave("t.pdf", final, width = 15, height = 3.2)
#############################################################


#################### Scenic analysis ########################
## p1
scenic <- read.table("input/Tcell_scenic.txt", header = T, sep = "\t")

tr <- c(
    "IRF1 (121g)", "RUNX3 (25g)", "EOMES (81g)",
    "KLF2_extended (70g)", "JUND_extended (1429g)", "JUN (19g)",
    "JUNB (20g)", "ETV7_extended (44g)", "BATF (81g)",
    "CREM (43g)", "NR3C1 (20g)", "ETS1_extended (182g)", "FOXP3_extended (30g)"
)
scenic <- scenic[tr, ]
mat <- data@meta.data[colnames(scenic), ]
mat <- subset(mat, celltype %in% c("GZMK+ CD8+ Tem1", "Treg", "CD4+ Tex"))
scenic <- scenic[, rownames(mat)]
mat$ggr <- paste(mat$subtype, mat$drug, mat$celltype, sep = "_")

mt <- matrix(NA, nrow = 13, ncol = 15)
for (i in 1:nrow(mt)) {
    mt[i, ] <- aggregate(as.vector(t(scenic)[, i]), list(mat$ggr), median)$x
}
colnames(mt) <- aggregate(as.vector(t(scenic)[, 1]), list(mat$ggr), median)$Group.1 # nolint
rownames(mt) <- rownames(scenic)
mat$tt <- paste(mat$subtype, mat$celltype, sep = "_")
mtt <- matrix(NA, nrow = 13, ncol = 9)
for (i in 1:nrow(mtt)) {
    mtt[i, ] <- aggregate(as.vector(t(scenic)[, i]), list(mat$tt), median)$x
}
colnames(mtt) <- aggregate(as.vector(t(scenic)[, 1]), list(mat$tt), median)$Group.1 # nolint
rownames(mtt) <- rownames(scenic)

mk <- data.frame(
    tr = rownames(mtt), tem1_a = mtt[, 8] / mtt[, 2],
    tem1_b = mt[, 11] / mt[, 5], tem1_c = mt[, 14] / mt[, 8],
    treg_a = mtt[, 9] / mtt[, 3], treg_b = mt[, 12] / mt[, 6],
    treg_c = mt[, 15] / mt[, 9], tex_a = mtt[, 7] / mtt[, 1],
    tex_b = mt[, 10] / mt[, 4], tex_c = mt[, 13] / mt[, 7]
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
ma <- matrix(NA, nrow = nrow(mk), ncol = 9)
for (i in 1:nrow(mk)) {
    ma[i, 1] <- wilcox.test(as.numeric(scenic[i, which(mat$tt == "RA_BT_GZMK+ CD8+ Tem1" & mat$celltype == "GZMK+ CD8+ Tem1")]), as.numeric(scenic[i, which(mat$tt == "OA_BT_GZMK+ CD8+ Tem1" & mat$celltype == "GZMK+ CD8+ Tem1")]))$p.value # nolint
    ma[i, 2] <- wilcox.test(as.numeric(scenic[i, which(mat$ggr == "RA_BT_Adalimumab_GZMK+ CD8+ Tem1" & mat$celltype == "GZMK+ CD8+ Tem1")]), as.numeric(scenic[i, which(mat$ggr == "RA_AT_Adalimumab_GZMK+ CD8+ Tem1" & mat$celltype == "GZMK+ CD8+ Tem1")]))$p.value # nolint
    ma[i, 3] <- wilcox.test(as.numeric(scenic[i, which(mat$ggr == "RA_BT_Tofacitinib_GZMK+ CD8+ Tem1" & mat$celltype == "GZMK+ CD8+ Tem1")]), as.numeric(scenic[i, which(mat$ggr == "RA_AT_Tofacitinib_GZMK+ CD8+ Tem1" & mat$celltype == "GZMK+ CD8+ Tem1")]))$p.value # nolint
    ma[i, 4] <- wilcox.test(as.numeric(scenic[i, which(mat$tt == "OA_BT_Treg" & mat$celltype == "Treg")]), as.numeric(scenic[i, which(mat$tt == "RA_BT_Treg" & mat$celltype == "Treg")]))$p.value # nolint
    ma[i, 5] <- wilcox.test(as.numeric(scenic[i, which(mat$ggr == "RA_BT_Adalimumab_Treg" & mat$celltype == "Treg")]), as.numeric(scenic[i, which(mat$ggr == "RA_AT_Adalimumab_Treg" & mat$celltype == "Treg")]))$p.value # nolint
    ma[i, 6] <- wilcox.test(as.numeric(scenic[i, which(mat$ggr == "RA_BT_Tofacitinib_Treg" & mat$celltype == "Treg")]), as.numeric(scenic[i, which(mat$ggr == "RA_AT_Tofacitinib_Treg" & mat$celltype == "Treg")]))$p.value # nolint
    ma[i, 7] <- wilcox.test(as.numeric(scenic[i, which(mat$tt == "RA_BT_CD4+ Tex" & mat$celltype == "CD4+ Tex")]), as.numeric(scenic[i, which(mat$tt == "OA_BT_CD4+ Tex" & mat$celltype == "CD4+ Tex")]))$p.value # nolint
    ma[i, 8] <- wilcox.test(as.numeric(scenic[i, which(mat$ggr == "RA_BT_Adalimumab_CD4+ Tex" & mat$celltype == "CD4+ Tex")]), as.numeric(scenic[i, which(mat$ggr == "RA_AT_Adalimumab_CD4+ Tex" & mat$celltype == "CD4+ Tex")]))$p.value # nolint
    ma[i, 9] <- wilcox.test(as.numeric(scenic[i, which(mat$ggr == "RA_BT_Tofacitinib_CD4+ Tex" & mat$celltype == "CD4+ Tex")]), as.numeric(scenic[i, which(mat$ggr == "RA_AT_Tofacitinib_CD4+ Tex" & mat$celltype == "CD4+ Tex")]))$p.value # nolint
}
colnames(ma) <- colnames(mk)[2:10] # nolint
rownames(ma) <- rownames(scenic)
mk <- melt(mk)
mk$pvalue <- melt(ma)$value
mk$value <- ifelse(mk$value > 1, 1, ifelse(mk$value < -1, -1, mk$value))

p <- ggplot(mk, aes(x = variable, y = tr)) +
    geom_point(aes(size = -log(pvalue, 10), color = value)) +
    scale_color_gradientn(colors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)) + # nolint
    geom_vline(xintercept = c(3.5, 6.5), color = "gray", linetype = "dashed") +
    theme_bw() +
    labs(x = "", y = "") +
    scale_x_discrete(labels = rep(c("RA-BT vs. OA-BT", "RA-BT vs. RA-AT\n(Adalimumab)", "RA-BT vs. RA-AT\n(Tofacitinib)"), 3)) + # nolint
    theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, color = "black", hjust = 1),
        axis.text.y = element_text(color = "black", size = 12)
    )
p
ggsave("Fig2/tf.pdf", p, width = 6.7, height = 4.5)

## p2
sub1 <- subset(data, ident = "GZMK+ CD8+ Tem1")
sub2 <- subset(data, ident = "Treg")
sub3 <- subset(data, ident = "CD4+ Tex")

tr <- c(
    "IRF1", "RUNX3", "EOMES", "KLF2", "JUND", "JUN",
    "JUNB", "ETV7", "BATF", "CREM", "NR3C1", "ETS1", "FOXP3"
)
corre1 <- function(x, y) {
    sdai_rr <- c()
    sdai_pp <- c()
    das28_rr <- c()
    das28_pp <- c()
    for (i in x) {
        print(i)
        mat <- as.data.frame(aggregate(as.numeric(y@assays$RNA@data[i, ]), list(y$sampleid), mean)) # nolint
        mat <- mat[-(1:3), ]
        result_sdai <- cor.test(mat$x, clinical_index_single$sdai_index[-(1:3)])
        sdai_r <- as.numeric(result_sdai$estimate)
        sdai_p <- as.numeric(result_sdai$p.value)
        result_das28 <- cor.test(mat$x, clinical_index_single$das28_index[-(1:3)]) # nolint
        das28_r <- as.numeric(result_das28$estimate)
        das28_p <- as.numeric(result_das28$p.value)
        sdai_rr <- append(sdai_rr, sdai_r)
        sdai_pp <- append(sdai_pp, sdai_p)
        das28_rr <- append(das28_rr, das28_r)
        das28_pp <- append(das28_pp, das28_p)
    }
    mt <- data.frame(sdai_r = sdai_rr, sdai_p = sdai_pp, das28_r = das28_rr, das28_p = das28_pp) # nolint
    rownames(mt) <- tr
    mt$gene <- tr
    return(mt)
}
tem1 <- corre1(tr, sub1)
treg <- corre1(tr, sub2)
cxcl13 <- corre1(tr, sub3)
sdai <- data.frame(row.names = tr, tem1 = tem1$sdai_r, treg = treg$sdai_r, cxcl13 = cxcl13$sdai_r) # nolint
das28 <- data.frame(row.names = tr, tem1 = tem1$das28_r, treg = treg$das28_r, cxcl13 = cxcl13$das28_r) # nolint

bk <- c(seq(-1, -0.1, by = 0.01), seq(0, 1, by = 0.01))
pdf("Fig2/sce_heat_sdai.pdf", width = 4.2, height = 1.5)
pheatmap(t(sdai),
    show_colnames = T,
    show_rownames = T, cluster_cols = F,
    cluster_rows = F, color = c(
        colorRampPalette(colors = brewer.pal(11, "RdBu")[11:6])(length(bk) / 2), # nolint
        colorRampPalette(colors = brewer.pal(11, "RdBu")[6:1])(length(bk) / 2)
    ),
    breaks = bk, annotation_legend = TRUE,
    legend_breaks = c(-1, 0, 1),
    fontsize = 6, annotation_names_col = FALSE,
    border_color = "white", treeheight_row = 8, treeheight_col = 8
)
dev.off()

corre <- function(x, y) {
    mat <- aggregate(as.numeric(y@assays$RNA@data[x, ]), list(y$sampleid), mean)
    mat <- cbind(mat, clinical_index_single)
    mat <- mat[-(1:3), ]

    p1 <- ggplot(mat) +
        geom_point(aes(x = x, y = das28_index, color = subtype), size = 3) +
        scale_color_manual(values = group_color3[2:3]) +
        geom_smooth(aes(x = x, y = das28_index), se = FALSE, method = "lm", color = "black") + # nolint
        stat_cor(aes(x = x, y = das28_index), method = "pearson") +
        labs(x = paste0("Expression of ", x), y = "DAS28") +
        theme_classic() +
        theme(
            axis.text = element_text(size = 10, color = "black"),
            axis.title = element_text(size = 12, color = "black"),
            legend.title = element_blank(),
            legend.text = element_text(size = 8, color = "black")
        )

    p2 <- ggplot(mat) +
        geom_point(aes(x = x, y = sdai_index, color = subtype), size = 3) +
        scale_color_manual(values = group_color3[2:3]) +
        geom_smooth(aes(x = x, y = sdai_index), se = FALSE, method = "lm", color = "black") + # nolint
        stat_cor(aes(x = x, y = sdai_index), method = "pearson") +
        labs(x = paste0("Expression of ", x), y = "SDAI") +
        theme_classic() +
        theme(
            axis.text = element_text(size = 10, color = "black"),
            axis.title = element_text(size = 12, color = "black"),
            legend.title = element_blank(),
            legend.text = element_text(size = 8, color = "black")
        )
    final <- p2 + p1 + plot_layout(nrow = 1, guides = "collect")
    final
    # ggsave("Fig3/STAT1_clinical_cor.pdf", final, width = 3.6, height = 5) # nolint
    ggsave("Fig4/clinical_cor.pdf", final, width = 5.5, height = 2.5)
}
corre("IRF1", sub2)
corre("IRF1", sub3)
corre("ETV7", sub3)
#############################################################


################## cellphonedb analysis #####################
## p1
meanvalue <- read.table("input/t_interaction_means.txt", header = T, sep = "\t", check.names = FALSE) # nolint
pvalue <- read.table("input/t_interaction_pvalues.txt", header = T, sep = "\t", check.names = FALSE) # nolint

group1 <- c(
    "GZMK+ CD8+ Tem1_OA_BT_None|Treg_OA_BT_None", "GZMK+ CD8+ Tem1_RA_BT_Adalimumab|Treg_RA_BT_Adalimumab", # nolint
    "GZMK+ CD8+ Tem1_RA_AT_Adalimumab|Treg_RA_AT_Adalimumab", "GZMK+ CD8+ Tem1_RA_BT_Tofacitinib|Treg_RA_BT_Tofacitinib", # nolint
    "GZMK+ CD8+ Tem1_RA_AT_Tofacitinib|Treg_RA_AT_Tofacitinib",
    "GZMK+ CD8+ Tem1_OA_BT_None|CD4+ Tex_OA_BT_None", "GZMK+ CD8+ Tem1_RA_BT_Adalimumab|CD4+ Tex_RA_BT_Adalimumab", # nolint
    "GZMK+ CD8+ Tem1_RA_AT_Adalimumab|CD4+ Tex_RA_AT_Adalimumab", "GZMK+ CD8+ Tem1_RA_BT_Tofacitinib|CD4+ Tex_RA_BT_Tofacitinib", # nolint
    "GZMK+ CD8+ Tem1_RA_AT_Tofacitinib|CD4+ Tex_RA_AT_Tofacitinib"
)
group2 <- c(
    "Treg_OA_BT_None|GZMK+ CD8+ Tem1_OA_BT_None", "Treg_RA_BT_Adalimumab|GZMK+ CD8+ Tem1_RA_BT_Adalimumab", # nolint
    "Treg_RA_AT_Adalimumab|GZMK+ CD8+ Tem1_RA_AT_Adalimumab", "Treg_RA_BT_Tofacitinib|GZMK+ CD8+ Tem1_RA_BT_Tofacitinib", # nolint
    "Treg_RA_AT_Tofacitinib|GZMK+ CD8+ Tem1_RA_AT_Tofacitinib",
    "Treg_OA_BT_None|CD4+ Tex_OA_BT_None", "Treg_RA_BT_Adalimumab|CD4+ Tex_RA_BT_Adalimumab", # nolint
    "Treg_RA_AT_Adalimumab|CD4+ Tex_RA_AT_Adalimumab", "Treg_RA_BT_Tofacitinib|CD4+ Tex_RA_BT_Tofacitinib", # nolint
    "Treg_RA_AT_Tofacitinib|CD4+ Tex_RA_AT_Tofacitinib"
)
group3 <- c(
    "CD4+ Tex_OA_BT_None|GZMK+ CD8+ Tem1_OA_BT_None", "CD4+ Tex_RA_BT_Adalimumab|GZMK+ CD8+ Tem1_RA_BT_Adalimumab", # nolint
    "CD4+ Tex_RA_AT_Adalimumab|GZMK+ CD8+ Tem1_RA_AT_Adalimumab", "CD4+ Tex_RA_BT_Tofacitinib|GZMK+ CD8+ Tem1_RA_BT_Tofacitinib", # nolint
    "CD4+ Tex_RA_AT_Tofacitinib|GZMK+ CD8+ Tem1_RA_AT_Tofacitinib",
    "CD4+ Tex_OA_BT_None|Treg_OA_BT_None", "CD4+ Tex_RA_BT_Adalimumab|Treg_RA_BT_Adalimumab", # nolint
    "CD4+ Tex_RA_AT_Adalimumab|Treg_RA_AT_Adalimumab", "CD4+ Tex_RA_BT_Tofacitinib|Treg_RA_BT_Tofacitinib", # nolint
    "CD4+ Tex_RA_AT_Tofacitinib|Treg_RA_AT_Tofacitinib"
)

lr1 <- c(
    "CCL5_CCR4", "CCL5_CCR5", "CCL4_CCR5", "CD8 receptor_LCK",
    "TNFSF12_TNFRSF25", "FAM3C_CLEC2D", "LAMP1_FAM3C", "SPP1_a4b1 complex"
)
lr2 <- c(
    "SELL_SELPLG", "TNFRSF1B_GRN", "CD27_CD70", "CD55_ADGRE5",
    "SPP1_PTGER4", "SPP1_CD44", "SPP1_a4b1 complex",
    "LGALS9_SORL1", "LGALS9_HAVCR2", "CCR6_CCL20"
)
lr3 <- c(
    "CD2_CD58", "SPP1_PTGER4", "CD55_ADGRE5", "PDCD1_FAM3C",
    "BTLA_TNFRSF14", "SELL_SELPLG", "NCR3_BAG6",
    "IGFL2_IGFLR1", "TNFSF12_TNFRSF25"
)

interac <- function(group = c(group2, group3),
                    lr = c(lr2, lr3),
                    lev = c(2, 8, 4, 10, 6, 1, 7, 3, 9, 5)) { # nolint
    means_g1 <- meanvalue[which(meanvalue$interacting_pair %in% lr), c(2, which(colnames(meanvalue) %in% group))] # nolint
    pvalue_g1 <- pvalue[which(pvalue$interacting_pair %in% lr), c(2, which(colnames(pvalue) %in% group))] # nolint
    pvalue_g1 <- melt(pvalue_g1)
    means_g1 <- melt(means_g1)
    pvalue_g1$value <- -log(pvalue_g1$value + 0.0001)
    pvalue_g1$value <- ifelse(pvalue_g1$value < 0, 0, pvalue_g1$value)
    means_g1$pvalue <- pvalue_g1$value
    means_g1$variable <- factor(means_g1$variable, levels = levels(means_g1$variable)[lev]) # nolint
    means_g1 <- means_g1[order(means_g1$variable), ] # nolint
    means_g1$interacting_pair <- factor(means_g1$interacting_pair, levels = rev(unique(lr))) # nolint

    p <- ggplot(means_g1, aes(x = variable, y = interacting_pair)) +
        geom_point(aes(size = pvalue, fill = value), color = "black", shape = 21) + # nolint
        scale_fill_gradient(low = "white", high = "#FF772F") +
        theme_bw() +
        scale_size(range = c(0, 8)) +
        theme(
            axis.text.x = element_blank(),
            axis.text.y = element_text(color = "black", size = 12),
            axis.ticks.x = element_blank(),
            axis.title = element_blank(),
            panel.grid = element_blank()
        )
    ggsave("Fig4/j_treg_tex_inte1.pdf", p, width = 6, height = 3.2)
}

interac(
    group = group1,
    lr = lr1,
    lev = c(2, 8, 4, 10, 6, 1, 7, 3, 9, 5)
)
interac(
    group = group2,
    lr = lr2,
    lev = c(2, 8, 4, 10, 6, 1, 7, 3, 9, 5)
)
interac(
    group = group3,
    lr = lr3,
    lev = c(1, 7, 3, 9, 5, 2, 8, 4, 10, 6)
)

## p2
count <- read.table("input/tcell_inter_counts.txt", header = T, sep = "\t")
count$group1 <- paste(
    sapply(strsplit(count$SOURCE, "_", fixed = T), "[", 2),
    sapply(strsplit(count$SOURCE, "_", fixed = T), "[", 3),
    sep = "_"
)
count$group2 <- paste(
    sapply(strsplit(count$TARGET, "_", fixed = T), "[", 2),
    sapply(strsplit(count$TARGET, "_", fixed = T), "[", 3),
    sep = "_"
)
count$SOURCE <- sapply(strsplit(count$SOURCE, "_", fixed = T), "[", 1)
count$TARGET <- sapply(strsplit(count$TARGET, "_", fixed = T), "[", 1)
oabt <- count[which(count$group1 == "OA_BT" & count$group2 == "OA_BT"), ]
rabt <- count[which(count$group1 == "RA_BT" & count$group2 == "RA_BT"), ]
raat <- count[which(count$group1 == "RA_AT" & count$group2 == "RA_AT"), ]

rabt <- dcast(rabt[, 1:3], SOURCE ~ TARGET)
rownames(rabt) <- rabt$SOURCE
rabt <- rabt[, -1]

raat <- dcast(raat[, 1:3], SOURCE ~ TARGET)
rownames(raat) <- raat$SOURCE
raat <- raat[, -1]

oabt <- dcast(oabt[, 1:3], SOURCE ~ TARGET)
rownames(oabt) <- oabt$SOURCE
oabt <- oabt[, -1]

type <- c("GZMK+ CD8+ Tem1", "Treg", "CD4+ Tex")
oabt <- oabt[type, type]
rabt <- rabt[type, type]
raat <- raat[type, type]

bk <- c(seq(0, 60, by = 1))
pdf("Fig4/q1.pdf", width = 2.6, height = 2)
pheatmap(oabt,
    cluster_rows = F, cluster_cols = F,
    scale = "none", border = "white",
    color = c(colorRampPalette(colors = brewer.pal(11, "RdBu")[11:1])(length(bk))), # nolint
    breaks = bk, treeheight_row = 0, treeheight_col = 0
)
dev.off()
pdf("Fig4/q2.pdf", width = 2.6, height = 2)
pheatmap(rabt,
    cluster_rows = F, cluster_cols = F,
    scale = "none", border = "white",
    color = c(colorRampPalette(colors = brewer.pal(11, "RdBu")[11:1])(length(bk))), # nolint
    breaks = bk, treeheight_row = 0, treeheight_col = 0
)
dev.off()
pdf("Fig4/q3.pdf", width = 2.6, height = 2)
pheatmap(raat,
    cluster_rows = F, cluster_cols = F,
    scale = "none", border = "white",
    color = c(colorRampPalette(colors = brewer.pal(11, "RdBu")[11:1])(length(bk))), # nolint
    breaks = bk, treeheight_row = 0, treeheight_col = 0
)
dev.off()

## p3
count <- read.table("input/tcell_drug_inter_counts1.txt", header = T, sep = "\t") # nolint
count$group1 <- paste(
    sapply(strsplit(count$SOURCE, "_", fixed = T), "[", 2),
    sapply(strsplit(count$SOURCE, "_", fixed = T), "[", 3),
    sapply(strsplit(count$SOURCE, "_", fixed = T), "[", 4),
    sep = "_"
)
count$group2 <- paste(
    sapply(strsplit(count$TARGET, "_", fixed = T), "[", 2),
    sapply(strsplit(count$TARGET, "_", fixed = T), "[", 3),
    sapply(strsplit(count$TARGET, "_", fixed = T), "[", 4),
    sep = "_"
)
count$SOURCE <- sapply(strsplit(count$SOURCE, "_", fixed = T), "[", 1)
count$TARGET <- sapply(strsplit(count$TARGET, "_", fixed = T), "[", 1)
oabt <- count[which(count$group1 == "OA_BT_None" & count$group2 == "OA_BT_None"), ] # nolint
rabt_a <- count[which(count$group1 == "RA_BT_Adalimumab" & count$group2 == "RA_BT_Adalimumab"), ] # nolint
rabt_t <- count[which(count$group1 == "RA_BT_Tofacitinib" & count$group2 == "RA_BT_Tofacitinib"), ] # nolint
raat_a <- count[which(count$group1 == "RA_AT_Adalimumab" & count$group2 == "RA_AT_Adalimumab"), ] # nolint
raat_t <- count[which(count$group1 == "RA_AT_Tofacitinib" & count$group2 == "RA_AT_Tofacitinib"), ] # nolint

rabt_a <- dcast(rabt_a[, 1:3], SOURCE ~ TARGET)
rownames(rabt_a) <- rabt_a$SOURCE
rabt_a <- rabt_a[, -1]
rabt_t <- dcast(rabt_t[, 1:3], SOURCE ~ TARGET)
rownames(rabt_t) <- rabt_t$SOURCE
rabt_t <- rabt_t[, -1]

raat_a <- dcast(raat_a[, 1:3], SOURCE ~ TARGET)
rownames(raat_a) <- raat_a$SOURCE
raat_a <- raat_a[, -1]
raat_t <- dcast(raat_t[, 1:3], SOURCE ~ TARGET)
rownames(raat_t) <- raat_t$SOURCE
raat_t <- raat_t[, -1]

oabt <- dcast(oabt[, 1:3], SOURCE ~ TARGET)
rownames(oabt) <- oabt$SOURCE
oabt <- oabt[, -1]

type <- c("GZMK+ CD8+ Tem1", "Treg", "CD4+ Tex")
oabt <- oabt[type, type]
rabt_a <- rabt_a[type, type]
rabt_t <- rabt_t[type, type]
raat_a <- raat_a[type, type]
raat_t <- raat_t[type, type]
ada <- rabt_a / raat_a
tof <- rabt_t / raat_t
for (i in 1:nrow(ada)) {
    for (j in 1:nrow(ada)) {
        ada[i, j] <- log(ada[i, j], 2)
    }
}
for (i in 1:nrow(ada)) {
    for (j in 1:nrow(ada)) {
        tof[i, j] <- log(tof[i, j], 2)
    }
}
mt <- cbind(ada, tof)

bk <- c(seq(-2, -0.01, by = 0.01), seq(0, 2, by = 0.01))
pdf("Fig4/p11.pdf", width = 3.6, height = 2)
pheatmap(mt,
    cluster_cols = FALSE, cluster_rows = FALSE, treeheight_row = 0,
    treeheight_col = 0, border_color = "white",
    color = c(
        colorRampPalette(colors = brewer.pal(11, "RdBu")[11:6])(length(bk) / 2),
        colorRampPalette(colors = brewer.pal(11, "RdBu")[6:1])(length(bk) / 2)
    ),
    breaks = bk
)
dev.off()
#############################################################


#################### DEG ACR20 analysis #####################
x <- "GZMK+ CD8+ Tem1"
sub <- subset(data_input, cells = rownames(data_input@meta.data)[which(data_input$subcell == x & data_input$ACR20 %in% c("N", "Y"))]) # nolint
sub@active.ident <- factor(sub$ACR20)
marker <- FindAllMarkers(sub, only.pos = T, min.pct = 0.25)
print(head(marker[which(marker$cluster == "N"), "gene"], 200))
print(head(marker[which(marker$cluster == "Y"), "gene"], 200))

genes <- c("CXCR4", "JUND", "MIF", "IFITM1", "ISG15")
genes <- c("MIF", "IFITM1")
genes <- c("CXCL13", "IFITM1", "MIF", "TIGIT", "GZMA", "NFKBIA", "IRF1")

mat <- as.data.frame(t(sub@assays$RNA@data[genes, ]))
mat$group <- sub$ACR20

plot_list <- foreach::foreach(gp = colnames(mat)[1:(ncol(mat) - 1)]) %do% {
    ggplot(mat, aes(x = group, y = get(gp), fill = group)) +
        geom_boxplot(outlier.colour = "white", outlier.size = 0) +
        scale_fill_manual(values = c("#87B1C8", "#A1568E")) +
        theme_classic() +
        labs(x = "", y = "Expression level", title = gp) +
        scale_y_continuous(expand = c(0, 0), limits = c(0, 1.2 * max(mat[, gp]))) + # nolint
        stat_compare_means(
            comparisons = list(
                c("Y", "N")
            ), size = 5, label = "p.signif"
        ) +
        scale_x_discrete(labels = c("ACR20_N", "ACR20_Y")) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 8), # nolint
            plot.title = element_text(hjust = 0.5, color = "black", size = 8, face = "bold"), # nolint
            axis.title.y = element_text(color = "black", size = 8),
            axis.text.y = element_text(color = "black", size = 6),
            axis.title.x = element_text(color = "black", size = 7)
        ) # nolint
}
final <- wrap_plots(plot_list, ncol = length(plot_list), guides = "collect")
ggsave("Fig4/CXCL13_deg.pdf", final, width = 5, height = 2)
#############################################################
