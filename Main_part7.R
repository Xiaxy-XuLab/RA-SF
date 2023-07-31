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
#############################################################


################## cellphonedb analysis #####################
## p1
meanvalue <- read.table("input/macro_tcell_inter_drug_means.txt", header = T, sep = "\t", check.names = FALSE) # nolint
pvalue <- read.table("input/macro_tcell_inter_drug_pvalues.txt", header = T, sep = "\t", check.names = FALSE) # nolint

group1 <- c(
    "SPP1+ Mac_OA_BT_None|GZMK+ CD8+ Tem1_OA_BT_None", "SPP1+ Mac_RA_BT_Adalimumab|GZMK+ CD8+ Tem1_RA_BT_Adalimumab", # nolint
    "SPP1+ Mac_RA_AT_Adalimumab|GZMK+ CD8+ Tem1_RA_AT_Adalimumab", "SPP1+ Mac_RA_BT_Tofacitinib|GZMK+ CD8+ Tem1_RA_BT_Tofacitinib", # nolint
    "SPP1+ Mac_RA_AT_Tofacitinib|GZMK+ CD8+ Tem1_RA_AT_Tofacitinib",
    "SPP1+ Mac_OA_BT_None|Treg_OA_BT_None", "SPP1+ Mac_RA_BT_Adalimumab|Treg_RA_BT_Adalimumab", # nolint
    "SPP1+ Mac_RA_AT_Adalimumab|Treg_RA_AT_Adalimumab", "SPP1+ Mac_RA_BT_Tofacitinib|Treg_RA_BT_Tofacitinib", # nolint
    "SPP1+ Mac_RA_AT_Tofacitinib|Treg_RA_AT_Tofacitinib",
    "SPP1+ Mac_OA_BT_None|CD4+ Tex_OA_BT_None", "SPP1+ Mac_RA_BT_Adalimumab|CD4+ Tex_RA_BT_Adalimumab", # nolint
    "SPP1+ Mac_RA_AT_Adalimumab|CD4+ Tex_RA_AT_Adalimumab", "SPP1+ Mac_RA_BT_Tofacitinib|CD4+ Tex_RA_BT_Tofacitinib", # nolint
    "SPP1+ Mac_RA_AT_Tofacitinib|CD4+ Tex_RA_AT_Tofacitinib"
)
group2 <- paste(sapply(strsplit(group1, "|", fixed = T), "[", 2), sapply(strsplit(group1, "|", fixed = T), "[", 1), sep = "|") # nolint

group3 <- c(
    "S100A12+ Mac_OA_BT_None|GZMK+ CD8+ Tem1_OA_BT_None", "S100A12+ Mac_RA_BT_Adalimumab|GZMK+ CD8+ Tem1_RA_BT_Adalimumab", # nolint
    "S100A12+ Mac_RA_AT_Adalimumab|GZMK+ CD8+ Tem1_RA_AT_Adalimumab", "S100A12+ Mac_RA_BT_Tofacitinib|GZMK+ CD8+ Tem1_RA_BT_Tofacitinib", # nolint
    "S100A12+ Mac_RA_AT_Tofacitinib|GZMK+ CD8+ Tem1_RA_AT_Tofacitinib",
    "S100A12+ Mac_OA_BT_None|Treg_OA_BT_None", "S100A12+ Mac_RA_BT_Adalimumab|Treg_RA_BT_Adalimumab", # nolint
    "S100A12+ Mac_RA_AT_Adalimumab|Treg_RA_AT_Adalimumab", "S100A12+ Mac_RA_BT_Tofacitinib|Treg_RA_BT_Tofacitinib", # nolint
    "S100A12+ Mac_RA_AT_Tofacitinib|Treg_RA_AT_Tofacitinib",
    "S100A12+ Mac_OA_BT_None|CD4+ Tex_OA_BT_None", "S100A12+ Mac_RA_BT_Adalimumab|CD4+ Tex_RA_BT_Adalimumab", # nolint
    "S100A12+ Mac_RA_AT_Adalimumab|CD4+ Tex_RA_AT_Adalimumab", "S100A12+ Mac_RA_BT_Tofacitinib|CD4+ Tex_RA_BT_Tofacitinib", # nolint
    "S100A12+ Mac_RA_AT_Tofacitinib|CD4+ Tex_RA_AT_Tofacitinib"
)
group4 <- paste(sapply(strsplit(group3, "|", fixed = T), "[", 2), sapply(strsplit(group3, "|", fixed = T), "[", 1), sep = "|") # nolint

lr1 <- c(
    "SPP1_CD44", "SPP1_PTGER4", "PLXNB2_SEMA4D",
    "CD74_COPA", "TNF_FAS", "ICAM1_SPN", "ICAM1_ITGAL",
    "LGALS9_SORL1", "TNF_VSIR", "ADRB2_VEGFB"
)
lr2 <- c(
    "TIGIT_NECTIN2", "TNFRSF1B_GRN", "CD2_CD58",
    "CD28_CD86", "CD6_ALCAM",
    "TFRC_TNFSF13B", "CCR5_CCL7",
    "IL15 receptor_IL15", "SPN_SIGLEC1"
)
lr3 <- c(
    "SELL_SELPLG", "C5AR1_RPS19", "TNF_FAS", "CD55_ADGRE5",
    "ICAM1_SPN", "TNFRSF10B_FASLG", "FLT1_VEGFB", "CCL4L2_VSIR",
    "CCL4_CCR5", "TNFRSF10A_TNFSF10", "TNFRSF10B_TNFSF10",
    "TNFSF14_TNFRSF14"
)
lr4 <- c(
    "TIGIT_NECTIN2", "HLA-F_LILRB2", "CXCR6_CXCL16",
    "CD2_CD58", "CD28_CD86", "CD6_ALCAM", "CD55_ADGRE5",
    "TNFRSF1B_GRN", "CCL5_CCR5", "CCR5_CCL7",
    "CCL4_SLC7A1", "CCL4_CCR5", "CD46_JAG1", "NOTCH1_JAG1",
    "CCR6_CCL20", "CXCR3_CCL20", "CD226_NECTIN2"
)

interac <- function(group = group4,
                    lr = lr4) {
    means_g1 <- meanvalue[which(meanvalue$interacting_pair %in% lr), c(2, which(colnames(meanvalue) %in% group))] # nolint
    pvalue_g1 <- pvalue[which(pvalue$interacting_pair %in% lr), c(2, which(colnames(pvalue) %in% group))] # nolint
    pvalue_g1 <- melt(pvalue_g1)
    means_g1 <- melt(means_g1)
    pvalue_g1$value <- -log(pvalue_g1$value + 0.0001)
    pvalue_g1$value <- ifelse(pvalue_g1$value < 0, 0, pvalue_g1$value)
    means_g1$pvalue <- pvalue_g1$value
    means_g1$variable <- factor(means_g1$variable, levels = group) # nolint
    means_g1 <- means_g1[order(means_g1$variable), ] # nolint
    means_g1$interacting_pair <- factor(means_g1$interacting_pair, levels = rev(lr)) # nolint
    means_g1$value <- ifelse(means_g1$value > 4, 4, means_g1$value)

    p <- ggplot(means_g1, aes(x = variable, y = interacting_pair)) +
        geom_point(aes(size = pvalue, fill = value), color = "black", shape = 21) + # nolint
        scale_fill_gradient(low = "white", high = "#FF772F") +
        theme_bw() +
        geom_vline(xintercept = c(5.5, 10.5), color = "gray", linetype = "dashed") + # nolint
        scale_size(range = c(0, 7)) +
        theme(
            axis.text.x = element_blank(),
            axis.text.y = element_text(color = "black", size = 12),
            axis.ticks.x = element_blank(),
            axis.title = element_blank(),
            panel.grid = element_blank()
        )
    ggsave("Fig5/interaction.pdf", p, width = 6.5, height = 4.6)
}
interac(group1, lr1)
interac(group2, lr2)
interac(group3, lr3)
interac(group4, lr4)

## p2
data <- read.table("input/macro_tcell_inter_count_network.txt", header = T, sep = "\t") # nolint
data$cellfrom <- sapply(strsplit(data$SOURCE, "_", fixed = T), "[", 1)
data$cellto <- sapply(strsplit(data$TARGET, "_", fixed = T), "[", 1)
data$from <- paste(
    sapply(strsplit(data$SOURCE, "_", fixed = T), "[", 2),
    sapply(strsplit(data$SOURCE, "_", fixed = T), "[", 3),
    sep = "-"
)
data$to <- paste(
    sapply(strsplit(data$TARGET, "_", fixed = T), "[", 2),
    sapply(strsplit(data$TARGET, "_", fixed = T), "[", 3),
    sep = "-"
)
data$fromgroup <- ifelse(data$cellfrom %in% c("SPP1+ Mac", "S100A12+ Mac", "C1QC+ Mac", "IEG Mac", "Mac", "CCL4+ Mac"), "Macrophage", "Tcell") # nolint
data$togroup <- ifelse(data$cellto %in% c("SPP1+ Mac", "S100A12+ Mac", "C1QC+ Mac", "IEG Mac", "Mac", "CCL4+ Mac"), "Macrophage", "Tcell") # nolint
data <- data[-which(data$fromgroup == data$togroup), ]

ddr <- function(x) {
    sub <- subset(data, fromgroup == x)
    dat <- dcast(sub[, 1:3], SOURCE ~ TARGET)
    rownames(dat) <- dat$SOURCE
    dat <- dat[, -1]
    if (x == "Tcell") {
        dat <- as.data.frame(t(dat))
    }
    dat <- dat[
        paste(rep(levels(macro), 3), rep(unique(macro$subtype), each = length(levels(macro))), sep = "_"), # nolint
        paste(rep(levels(tcell), 3), rep(unique(tcell$subtype), each = length(levels(tcell))), sep = "_"), # nolint
    ]
    dat_oa <- dat[1:6, 1:10]
    dat_oa <- dat_oa[1:5, 1:8]
    dat_bt <- dat[7:12, 11:20]
    dat_bt <- dat_bt[1:5, 1:8]
    dat_at <- dat[13:18, 21:30]
    dat_at <- dat_at[1:5, 1:8]

    mk <- (dat_bt - dat_at) / dat_bt
    bk <- c(seq(-0.4, -0.01, by = 0.01), seq(0, 0.6, by = 0.01))
    p0 <- pheatmap(
        mk,
        scale = "none",
        color = c(
            colorRampPalette(colors = brewer.pal(11, "RdBu")[11:6])(length(bk) * 2 / 5), # nolint
            colorRampPalette(colors = brewer.pal(11, "RdBu")[6:1])(length(bk) * 3 / 5) # nolint
        ),
        cluster_rows = F, cluster_cols = F,
        show_rownames = F, show_colnames = F,
        breaks = bk, border = "white"
    )
    ggsave("Fig5/remission.pdf", p0, width = 3.8, height = 1.7)

    bk <- c(seq(0, max(dat_oa, dat_bt, dat_at), by = 1))
    p1 <- pheatmap(
        dat_oa,
        scale = "none",
        color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(83),
        cluster_rows = F, cluster_cols = F,
        show_rownames = F, show_colnames = F,
        breaks = bk, border = "white"
    )

    p2 <- pheatmap(
        dat_bt,
        scale = "none",
        color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(83),
        cluster_rows = F, cluster_cols = F,
        show_rownames = F, show_colnames = F,
        breaks = bk, border = "white"
    )

    p3 <- pheatmap(
        dat_at,
        scale = "none",
        color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(83),
        cluster_rows = F, cluster_cols = F,
        show_rownames = F, show_colnames = F,
        breaks = bk, border = "white"
    )
    final <- as.ggplot(p1) + as.ggplot(p2) + as.ggplot(p3) + plot_layout(nrow = 3, guides = "collect") # nolint
    return(final)
}
final_draw <- ddr("Macrophage")
ggsave("Fig5/Macrophage_tcell.pdf", final, width = 3.8, height = 5)

pdf("Fig5/row_label.pdf", width = 3.6, height = 0.5)
show_col(tcell_color, ncol = 10, labels = FALSE, borders = "white")
dev.off()

pdf("Fig5/col_label.pdf", width = 0.5, height = 1.5)
show_col(macrophage_color, ncol = 1, labels = FALSE, borders = "white")
dev.off()

## p3
data <- read.table("input/ACR20_count.txt", header = T, sep = "\t")
data$cellfrom <- sapply(strsplit(data$SOURCE, "-", fixed = T), "[", 1)
data$cellto <- sapply(strsplit(data$TARGET, "-", fixed = T), "[", 1)
data$from <- sapply(strsplit(data$SOURCE, "-", fixed = T), "[", 2)
data$to <- sapply(strsplit(data$TARGET, "-", fixed = T), "[", 2)
data <- subset(
    data,
    cellfrom %in% unique(data$cellfrom)[c(1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12, 17, 24)] & # nolint
        cellto %in% unique(data$cellfrom)[c(1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12, 17, 24)] # nolint
)

data$fromgroup <- ifelse(data$cellfrom %in% c("SPP1+ Mac", "S100A12+ Mac", "C1QC+ Mac", "IEG Mac", "Mac"), "Macrophage", "Tcell") # nolint
data$togroup <- ifelse(data$cellto %in% c("SPP1+ Mac", "S100A12+ Mac", "C1QC+ Mac", "IEG Mac", "Mac"), "Macrophage", "Tcell") # nolint
data <- data[-which(data$fromgroup == data$togroup), ]
ddr <- function(x) {
    sub <- subset(data, fromgroup == x)
    dat <- dcast(sub[, 1:3], SOURCE ~ TARGET)
    rownames(dat) <- dat$SOURCE
    dat <- dat[, -1]
    if (x == "Tcell") {
        dat <- as.data.frame(t(dat))
    }
    dat <- dat[
        paste(rep(levels(macro)[1:5], 2), rep(c("N", "Y"), each = 5), sep = "-"), # nolint
        paste(rep(levels(tcell)[1:8], 2), rep(c("N", "Y"), each = 8), sep = "-"), # nolint
    ]
    dat_oa <- dat[1:5, 1:8]
    dat_bt <- dat[6:10, 9:16]

    bk <- c(seq(0, 83, by = 1))
    p1 <- pheatmap(
        dat_oa,
        scale = "none",
        color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(length(bk)), # nolint
        cluster_rows = F, cluster_cols = F,
        show_rownames = F, show_colnames = F,
        breaks = bk, border = "white"
    )

    p2 <- pheatmap(
        dat_bt,
        scale = "none",
        color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(length(bk)), # nolint
        cluster_rows = F, cluster_cols = F,
        show_rownames = F, show_colnames = F,
        breaks = bk, border = "white"
    )

    final <- as.ggplot(p1) + as.ggplot(p2) + plot_layout(nrow = 2, guides = "collect") # nolint
    return(final)
}
final_draw <- ddr("Macrophage")
ggsave("Fig5/Macrophage_tcell.pdf", final_draw, width = 3.8, height = 3.4)

## p4
meanvalue <- read.table("input/spp1_tcell_inter_drug_acr20_means.txt", header = T, sep = "\t", check.names = FALSE) # nolint
pvalue <- read.table("input/spp1_tcell_inter_drug_acr20_pvalues.txt", header = T, sep = "\t", check.names = FALSE) # nolint

group1 <- c(
    "SPP1+ Mac-N|GZMK+ CD8+ Tem1-N", "SPP1+ Mac-Y|GZMK+ CD8+ Tem1-Y",
    "SPP1+ Mac-N|Treg-N", "SPP1+ Mac-Y|Treg-Y",
    "SPP1+ Mac-N|CD4+ Tex-N", "SPP1+ Mac-Y|CD4+ Tex-Y"
)
group2 <- paste(sapply(strsplit(group1, "|", fixed = T), "[", 2), sapply(strsplit(group1, "|", fixed = T), "[", 1), sep = "|") # nolint

group3 <- c(
    "S100A12+ Mac-N|GZMK+ CD8+ Tem1-N", "S100A12+ Mac-Y|GZMK+ CD8+ Tem1-Y",
    "S100A12+ Mac-N|Treg-N", "S100A12+ Mac-Y|Treg-Y",
    "S100A12+ Mac-N|CD4+ Tex-N", "S100A12+ Mac-Y|CD4+ Tex-Y"
)
group4 <- paste(sapply(strsplit(group3, "|", fixed = T), "[", 2), sapply(strsplit(group3, "|", fixed = T), "[", 1), sep = "|") # nolint

lr1 <- c(
    "CD74_MIF", "CD74_COPA",
    "CD47_SIRPG", "TNFRSF10B_TNFSF10", "TNF_RIPK1",
    "IL1B_ADRB2", "TNFSF12_TNFRSF25",
    "CSF1R_CSF1", "VEGFA_FLT1", "FLT1_VEGFB",
    "TNF_VSIR", "CD40_CD40LG",
    "TNFSF14_TNFRSF14",
    "NRP1_VEGFB", "CCL4_CCR5",
    "TNF_ICOS", "CCL4L2_VSIR", "IL7 receptor_IL7"
)
lr2 <- c(
    "CD2_CD58", "CD28_CD86", "CD6_ALCAM",
    "TIGIT_NECTIN2", "PGRMC2_CCL4L2", "TNFRSF1B_GRN",
    "HLA-C_FAM3C", "TNFRSF1A_GRN", "TFRC_TNFSF13B",
    "ADRB2_VEGFB", "IL15 receptor_IL15", "CCR5_CCL7",
    "SPN_SIGLEC1", "CD52_SIGLEC10", "CSF1_SIRPA", "SELL_SELPLG"
)
lr3 <- c(
    "CD55_ADGRE5",
    "FLT1_VEGFB",
    "TNF_RIPK1", "SIRPA_CD47",
    "TGFB1_TGFbeta receptor1",
    "IL1B_ADRB2", "CD47_SIRPG", "CD40_CD40LG",
    "VEGFA_FLT1", "SEMA4A_PLXND1",
    "TNF_VSIR", "TNF_NOTCH1", "ICOSLG_ICOS",
    "TNF_ICOS", "TNFSF12_TNFRSF25", "LGALS9_HAVCR2"
)
lr4 <- c(
    "CD52_SIGLEC10", "LGALS9_LRP1", "CD99_PILRA", "TNFSF14_LTBR",
    "IFNG_Type II IFNR", "SPN_SIGLEC1", "PDCD1_FAM3C",
    "HLA-C_FAM3C", "ADRB2_VEGFB", "TNFSF14_TNFRSF14", "CSF1_SLC7A1",
    "PGRMC2_CCL4L2", "CXCR3_CCL20", "CTLA4_CD86",
    "CD46_JAG1", "CCR6_CCL20", "MIF_TNFRSF14", "TNF_TNFRSF1A",
    "TNF_TNFRSF1B", "TNF_NOTCH1"
)

interac <- function(group = group4,
                    lr = lr4) {
    means_g1 <- meanvalue[which(meanvalue$interacting_pair %in% lr), c(2, which(colnames(meanvalue) %in% group))] # nolint
    pvalue_g1 <- pvalue[which(pvalue$interacting_pair %in% lr), c(2, which(colnames(pvalue) %in% group))] # nolint
    pvalue_g1 <- melt(pvalue_g1)
    means_g1 <- melt(means_g1)
    pvalue_g1$value <- -log(pvalue_g1$value + 0.0001)
    pvalue_g1$value <- ifelse(pvalue_g1$value < 0, 0, pvalue_g1$value)
    means_g1$pvalue <- pvalue_g1$value
    means_g1$variable <- factor(means_g1$variable, levels = group) # nolint
    means_g1 <- means_g1[order(means_g1$variable), ] # nolint
    means_g1$interacting_pair <- factor(means_g1$interacting_pair, levels = rev(lr)) # nolint
    means_g1$value <- ifelse(means_g1$value > 4, 4, means_g1$value)

    p <- ggplot(means_g1, aes(x = variable, y = interacting_pair)) +
        geom_point(aes(size = pvalue, fill = value), color = "black", shape = 21) + # nolint
        scale_fill_gradient(low = "white", high = "#FF772F") +
        theme_bw() +
        geom_vline(xintercept = c(2.5, 4.5), color = "gray", linetype = "dashed") + # nolint
        scale_size(range = c(0, 7)) +
        theme(
            axis.text.x = element_blank(),
            axis.text.y = element_text(color = "black", size = 12),
            axis.ticks.x = element_blank(),
            axis.title = element_blank(),
            panel.grid = element_blank()
        )
    ggsave("Fig5/ddd.pdf", p, width = 5, height = 5.4)
}

interac(group1, lr1)
interac(group2, lr2)

meanvalue <- read.table("input/s100a12_tcell_inter_drug_acr20_means.txt", header = T, sep = "\t", check.names = FALSE) # nolint
pvalue <- read.table("input/s100a12_tcell_inter_drug_acr20_pvalues.txt", header = T, sep = "\t", check.names = FALSE) # nolint
interac(group3, lr3)
interac(group4, lr4)
#############################################################


#################### circle analysis ########################
da <- read.table("input/macro_tcell_inter_drug_count_network1.txt", header = T, sep = "\t") # nolint
da$sg <- paste(
    sapply(strsplit(da$SOURCE, "_", fixed = T), "[", 2),
    sapply(strsplit(da$SOURCE, "_", fixed = T), "[", 3),
    sapply(strsplit(da$SOURCE, "_", fixed = T), "[", 4),
    sep = "_"
)
da$tg <- paste(
    sapply(strsplit(da$TARGET, "_", fixed = T), "[", 2),
    sapply(strsplit(da$TARGET, "_", fixed = T), "[", 3),
    sapply(strsplit(da$TARGET, "_", fixed = T), "[", 4),
    sep = "_"
)
da <- da[which(da$sg == da$tg), ]
da$source_cell <- sapply(strsplit(da$SOURCE, "_", fixed = T), "[", 1)
da$target_cell <- sapply(strsplit(da$TARGET, "_", fixed = T), "[", 1)
da <- subset(da, source_cell == "SPP1+ Mac" & target_cell %in% c("GZMK+ CD8+ Tem1", "Treg", "CD4+ Tex")) # nolint
da <- subset(da, source_cell == "S100A12+ Mac" & target_cell %in% c("GZMK+ CD8+ Tem1", "Treg", "CD4+ Tex")) # nolint

dt1 <- dcast(da[, 1:3], SOURCE ~ TARGET)
rownames(dt1) <- dt1$SOURCE
dt1 <- dt1[, -1]
dtt1 <- na.fill(dt1, 0)
rownames(dtt1) <- rownames(dt1)
dtt1 <- dtt1[c(1, 4, 2, 3), c(5, 8, 6, 7, 9, 12, 10, 11, 1, 4, 2, 3)] # nolint
mt1 <- dtt1[, 1:4]
mt2 <- dtt1[, 5:8]
mt3 <- dtt1[, 9:12]

group <- structure(
    c(
        rep("SPP1+ Mac", 4), rep("GZMK+ CD8+ Tem1", 4),
        rep("Treg", 4), rep("CD4+ Tex", 4)
    ),
    names = c(rownames(dtt1), colnames(dtt1))
)
group <- structure(
    c(
        rep("S100A12+ Mac", 4), rep("GZMK+ CD8+ Tem1", 4),
        rep("Treg", 4), rep("CD4+ Tex", 4)
    ),
    names = c(rownames(dtt1), colnames(dtt1))
)

grid_col <- structure(
    rep(pal_lancet(alpha = 0.5)(9)[c(4, 2, 6, 5)], 4),
    names = c(
        rownames(dtt1), colnames(dtt1)
    )
)

circos.clear()
pdf("Fig5/circle_S100A12.pdf", width = 4, height = 4)
# 最外层添加一个空白轨迹
chordDiagram(
    dtt1,
    group = group, grid.col = grid_col,
    annotationTrack = c("grid", "axis"),
    preAllocateTracks = list(
        track.height = mm_h(4),
        track.margin = c(mm_h(4), 0)
    )
)

# 将扇形标签放置在格子中
circos.track(
    track.index = 2,
    panel.fun = function(x, y) {
        sector.index <- get.cell.meta.data("sector.index")
        xlim <- get.cell.meta.data("xlim")
        ylim <- get.cell.meta.data("ylim")
        circos.text(
            mean(xlim), mean(ylim),
            sector.index,
            cex = 0.6,
            niceFacing = TRUE
        )
    },
    bg.border = NA
)
# 高亮最外层的扇形格子
highlight.sector(
    rownames(mt1),
    track.index = 1, col = "#F1B26E", # "#C6B4D3",
    text = "S100A12+ Mac", cex = 0.8, text.col = "white",
    niceFacing = TRUE
)
highlight.sector(
    colnames(mt1),
    track.index = 1, col = "#E64B35E5",
    text = "GZMK+ CD8+ Tem1", cex = 0.8, text.col = "white",
    niceFacing = TRUE
)
highlight.sector(
    colnames(mt2),
    track.index = 1, col = "#91D1C2E5",
    text = "Treg", cex = 0.8, text.col = "white",
    niceFacing = TRUE
)
highlight.sector(
    colnames(mt3),
    track.index = 1, col = "#00A087E5",
    text = "CD4+ Tex", cex = 0.8, text.col = "white",
    niceFacing = TRUE
)
dev.off()
#############################################################


#################### nichenet analysis ######################
## nichenet
ligand_target_matrix <- readRDS("/work/xiaxy/work/RA/database/nichenet/ligand_target_matrix.rds") # nolint
lr_network <- readRDS("/work/xiaxy/work/RA/database/nichenet/lr_network.rds") # nolint
weighted_networks <- readRDS("/work/xiaxy/work/RA/database/nichenet/weighted_networks.rds") # nolint
weighted_networks_lr <- weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from, to), by = c("from", "to")) # nolint

marker <- function(x) {
    sub <- subset(tcell, ident = x)
    sub@active.ident <- factor(sub$subtype)
    mark <- rownames(FindMarkers(sub, ident.1 = "RA_BT", ident.2 = "OA_BT", only.pos = T, min.pct = 0.25)) # nolint
    mark <- mark[!grepl("^MT-|XIST", mark)]
    return(mark[1:100])
}

receiver_genes <- unique(c(marker("GZMK+ CD8+ Tem1"), marker("Treg"), marker("CD4+ Tex"))) # nolint

tcell@active.ident <- factor(paste(tcell$celltype, tcell$subtype, sep = "_"))
names(tcell@active.ident) <- rownames(tcell@meta.data)
receiver <- c("GZMK+ CD8+ Tem1_RA_BT", "Treg_RA_BT", "CD4+ Tex_RA_BT")
receobject <- subset(tcell, ident = receiver)
expressed_genes_receiver <- receiver %>%
    unique() %>%
    lapply(get_expressed_genes, receobject, 0.10)
expressed_genes_receiver <- expressed_genes_receiver %>%
    unlist() %>%
    unique()
background_expressed_genes <- expressed_genes_receiver %>%
    .[. %in% rownames(ligand_target_matrix)]

sender_celltypes <- c("GZMK+ CD8+ Tem1_RA_BT", "Treg_RA_BT", "CD4+ Tex_RA_BT")
tcell@active.ident <- factor(paste(tcell$celltype, tcell$subtype, sep = "_"))
names(tcell@active.ident) <- rownames(tcell@meta.data)
sendobject <- subset(tcell, ident = sender_celltypes)
list_expressed_genes_sender <- sender_celltypes %>%
    unique() %>%
    lapply(get_expressed_genes, sendobject, 0.10)
expressed_genes_sender <- list_expressed_genes_sender %>%
    unlist() %>%
    unique()

geneset_oi <- receiver_genes
geneset_oi <- geneset_oi %>%
    .[. %in% rownames(ligand_target_matrix)]

ligands <- lr_network %>%
    pull(from) %>%
    unique()
receptors <- lr_network %>%
    pull(to) %>%
    unique()

expressed_ligands <- intersect(ligands, expressed_genes_sender)
expressed_receptors <- intersect(receptors, expressed_genes_receiver)

potential_ligands <- lr_network %>%
    filter(from %in% expressed_ligands & to %in% expressed_receptors) %>%
    pull(from) %>%
    unique()
ligand_activities <- predict_ligand_activities(
    geneset = geneset_oi,
    background_expressed_genes = background_expressed_genes,
    ligand_target_matrix = ligand_target_matrix,
    potential_ligands = potential_ligands
)

ligand_activities <- ligand_activities %>%
    arrange(-pearson) %>%
    mutate(rank = rank(desc(pearson)))

best_upstream_ligands <- ligand_activities %>%
    top_n(60, pearson) %>%
    arrange(-pearson) %>%
    pull(test_ligand) %>%
    unique()

active_ligand_target_links_df <- best_upstream_ligands %>%
    lapply(get_weighted_ligand_target_links,
        geneset = geneset_oi,
        ligand_target_matrix = ligand_target_matrix, n = 200
    ) %>%
    bind_rows() %>%
    drop_na()
active_ligand_target_links <- prepare_ligand_target_visualization(
    ligand_target_df = active_ligand_target_links_df,
    ligand_target_matrix = ligand_target_matrix, cutoff = 0.33
)
order_ligands <- intersect(
    best_upstream_ligands,
    colnames(active_ligand_target_links)
) %>%
    rev() %>%
    make.names()
order_targets <- active_ligand_target_links_df$target %>%
    unique() %>%
    intersect(rownames(active_ligand_target_links)) %>%
    make.names()
rownames(active_ligand_target_links) <- rownames(active_ligand_target_links) %>% # nolint
    make.names() # make.names() for heatmap visualization of genes like H2-T23
colnames(active_ligand_target_links) <- colnames(active_ligand_target_links) %>%
    make.names() # make.names() for heatmap visualization of genes like H2-T23

vis_ligand_target <- active_ligand_target_links[order_targets, order_ligands] %>% t() # nolint
vis_ligand_target <- vis_ligand_target[-c(3, 4, 5, 6, 7, 8, 10, 11, 20, 22, 24, 26, 28), ] # nolint
rownames(vis_ligand_target) <- gsub("[.]", "-", rownames(vis_ligand_target))

bk <- c(seq(0, 0.009, by = 0.00001))
colo <- colorRampPalette(brewer.pal(n = 9, name = "OrRd"))(100)
p2 <- pheatmap(vis_ligand_target,
    cluster_cols = FALSE, cluster_rows = FALSE, treeheight_row = 0,
    treeheight_col = 0, border_color = "black",
    color = rev(c(
        colorRampPalette(colors = colo[100:81])(length(bk) / 2), # noilnt
        colorRampPalette(colors = colo[80:1])(length(bk) / 2) # nolint
    ))
)
ggsave("Fig5/Tcell_Macrophage_lt-net.pdf", p2, width = 7.8, height = 5.2)

gge <- rownames(vis_ligand_target)
gge <- gsub("[.]", "-", gge)
t_tcell <- tcell
t_tcell@active.ident <- factor(paste(t_tcell$celltype, t_tcell$subtype, t_tcell$drug, sep = "_")) # nolint
names(t_tcell@active.ident) <- rownames(t_tcell@meta.data)
sub <- subset(t_tcell, ident = c(
    "GZMK+ CD8+ Tem1_OA_BT_None", "GZMK+ CD8+ Tem1_RA_BT_Adalimumab",
    "GZMK+ CD8+ Tem1_RA_AT_Adalimumab", "GZMK+ CD8+ Tem1_RA_BT_Tofacitinib",
    "GZMK+ CD8+ Tem1_RA_AT_Tofacitinib", "Treg_OA_BT_None",
    "Treg_RA_BT_Adalimumab", "Treg_RA_AT_Adalimumab",
    "Treg_RA_BT_Tofacitinib", "Treg_RA_AT_Tofacitinib",
    "CD4+ Tex_OA_BT_None", "CD4+ Tex_RA_BT_Adalimumab",
    "CD4+ Tex_RA_AT_Adalimumab", "CD4+ Tex_RA_BT_Tofacitinib",
    "CD4+ Tex_RA_AT_Tofacitinib"
))
sub$kk <- sub@active.ident
mat <- sub@assays$RNA@data[gge, ]
mt <- aggregate(t(mat), list(sub$kk), mean)
rownames(mt) <- mt$Group.1
mt <- mt[, -1]
tt <- as.data.frame(t(mt))
tt <- tt[, c(6, 9, 7, 10, 8, 11, 14, 12, 15, 13, 1, 4, 2, 5, 3)]

bk <- c(seq(-1.5, -0.01, by = 0.01), seq(0, 1.5, by = 0.01))
colo2 <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)
p1 <- pheatmap(tt,
    scale = "row",
    cluster_cols = FALSE, cluster_rows = FALSE, treeheight_row = 0,
    treeheight_col = 0, border_color = "black",
    color = rev(c(
        colorRampPalette(colors = colo2[100:51])(length(bk) / 2),
        colorRampPalette(colors = colo2[50:1])(length(bk) / 2)
    )),
    breaks = bk, gaps_col = c(5, 10)
)
ggsave("Fig5/Tcell_Macrophage_ligand_exp.pdf", p1, width = 2.5, height = 6)

ggt <- colnames(vis_ligand_target)
ggt <- gsub("[.]", "-", ggt)
t_macro <- macro
t_macro@active.ident <- factor(paste(t_macro$celltype, t_macro$subtype, t_macro$drug, sep = "_")) # nolint
names(t_macro@active.ident) <- rownames(t_macro@meta.data)
sub <- subset(t_macro, ident = c(
    "SPP1+ Mac_OA_BT_None", "SPP1+ Mac_RA_BT_Adalimumab",
    "SPP1+ Mac_RA_AT_Adalimumab", "SPP1+ Mac_RA_BT_Tofacitinib",
    "SPP1+ Mac_RA_AT_Tofacitinib", "S100A12+ Mac_OA_BT_None",
    "S100A12+ Mac_RA_BT_Adalimumab", "S100A12+ Mac_RA_AT_Adalimumab",
    "S100A12+ Mac_RA_BT_Tofacitinib", "S100A12+ Mac_RA_AT_Tofacitinib"
))
sub$kk <- sub@active.ident
mat <- sub@assays$RNA@data[ggt, ]
mt <- aggregate(t(mat), list(sub$kk), mean)
rownames(mt) <- mt$Group.1
mt <- mt[, -1]
tt <- as.data.frame(t(mt))
tt <- tt[, c(6, 9, 7, 10, 8, 1, 4, 2, 5, 3)]

bk <- c(seq(-1.5, -0.01, by = 0.01), seq(0, 1.5, by = 0.01))
colo2 <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)
p1 <- pheatmap(t(tt),
    scale = "column",
    cluster_cols = FALSE, cluster_rows = FALSE, treeheight_row = 0,
    treeheight_col = 0, border_color = "black",
    color = rev(c(
        colorRampPalette(colors = colo2[100:51])(length(bk) / 2),
        colorRampPalette(colors = colo2[50:1])(length(bk) / 2)
    )),
    breaks = bk, gaps_row = c(5)
)
ggsave("Fig5/Tcell_Macrophage_receiver_exp.pdf", p1, width = 6, height = 2.1)

###
ligand_target_matrix <- readRDS("/work/xiaxy/work/RA/database/nichenet/ligand_target_matrix.rds") # nolint
lr_network <- readRDS("/work/xiaxy/work/RA/database/nichenet/lr_network.rds") # nolint
weighted_networks <- readRDS("/work/xiaxy/work/RA/database/nichenet/weighted_networks.rds") # nolint
weighted_networks_lr <- weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from, to), by = c("from", "to")) # nolint

marker <- function(x) {
    sub <- subset(macro, ident = x)
    sub@active.ident <- factor(sub$subtype)
    mark <- rownames(FindMarkers(sub, ident.1 = "RA_BT", ident.2 = "OA_BT", only.pos = T, min.pct = 0.25)) # nolint
    mark <- mark[!grepl("^MT-|XIST", mark)]
    return(mark[1:100])
}

receiver_genes <- unique(c(marker("SPP1+ Mac"), marker("S100A12+ Mac"))) # nolint

macro@active.ident <- factor(paste(macro$celltype, macro$subtype, sep = "_"))
names(macro@active.ident) <- rownames(macro@meta.data)
receiver <- c("SPP1+ Mac_RA_BT", "S100A12+ Mac_RA_BT")
receobject <- subset(macro, ident = receiver)
expressed_genes_receiver <- receiver %>%
    unique() %>%
    lapply(get_expressed_genes, receobject, 0.10)
expressed_genes_receiver <- expressed_genes_receiver %>%
    unlist() %>%
    unique()
background_expressed_genes <- expressed_genes_receiver %>%
    .[. %in% rownames(ligand_target_matrix)]

sender_celltypes <- c("SPP1+ Mac_RA_BT", "S100A12+ Mac_RA_BT")
macro@active.ident <- factor(paste(macro$celltype, macro$subtype, sep = "_"))
names(macro@active.ident) <- rownames(macro@meta.data)
sendobject <- subset(macro, ident = sender_celltypes)
list_expressed_genes_sender <- sender_celltypes %>%
    unique() %>%
    lapply(get_expressed_genes, sendobject, 0.10)
expressed_genes_sender <- list_expressed_genes_sender %>%
    unlist() %>%
    unique()

geneset_oi <- receiver_genes
geneset_oi <- geneset_oi %>%
    .[. %in% rownames(ligand_target_matrix)]

ligands <- lr_network %>%
    pull(from) %>%
    unique()
receptors <- lr_network %>%
    pull(to) %>%
    unique()

expressed_ligands <- intersect(ligands, expressed_genes_sender)
expressed_receptors <- intersect(receptors, expressed_genes_receiver)

potential_ligands <- lr_network %>%
    filter(from %in% expressed_ligands & to %in% expressed_receptors) %>%
    pull(from) %>%
    unique()
ligand_activities <- predict_ligand_activities(
    geneset = geneset_oi,
    background_expressed_genes = background_expressed_genes,
    ligand_target_matrix = ligand_target_matrix,
    potential_ligands = potential_ligands
)

ligand_activities <- ligand_activities %>%
    arrange(-pearson) %>%
    mutate(rank = rank(desc(pearson)))

best_upstream_ligands <- ligand_activities %>%
    top_n(60, pearson) %>%
    arrange(-pearson) %>%
    pull(test_ligand) %>%
    unique()

active_ligand_target_links_df <- best_upstream_ligands %>%
    lapply(get_weighted_ligand_target_links,
        geneset = geneset_oi,
        ligand_target_matrix = ligand_target_matrix, n = 200
    ) %>%
    bind_rows() %>%
    drop_na()
active_ligand_target_links <- prepare_ligand_target_visualization(
    ligand_target_df = active_ligand_target_links_df,
    ligand_target_matrix = ligand_target_matrix, cutoff = 0.33
)
order_ligands <- intersect(
    best_upstream_ligands,
    colnames(active_ligand_target_links)
) %>%
    rev() %>%
    make.names()
order_targets <- active_ligand_target_links_df$target %>%
    unique() %>%
    intersect(rownames(active_ligand_target_links)) %>%
    make.names()
rownames(active_ligand_target_links) <- rownames(active_ligand_target_links) %>% # nolint
    make.names() # make.names() for heatmap visualization of genes like H2-T23
colnames(active_ligand_target_links) <- colnames(active_ligand_target_links) %>%
    make.names() # make.names() for heatmap visualization of genes like H2-T23

vis_ligand_target <- active_ligand_target_links[order_targets, order_ligands] %>% t() # nolint
vis_ligand_target <- vis_ligand_target[-c(1, 3, 4, 5, 6, 9, 11, 12, 13, 17, 18, 20, 24, 31), ] # nolint
rownames(vis_ligand_target) <- gsub("[.]", "-", rownames(vis_ligand_target))

bk <- c(seq(0, 0.006, by = 0.00001))
colo <- colorRampPalette(brewer.pal(n = 9, name = "OrRd"))(100)
p2 <- pheatmap(vis_ligand_target,
    cluster_cols = FALSE, cluster_rows = FALSE, treeheight_row = 0,
    treeheight_col = 0, border_color = "black",
    color = rev(c(
        colorRampPalette(colors = colo[100:81])(length(bk) / 2), # noilnt
        colorRampPalette(colors = colo[80:1])(length(bk) / 2) # nolint
    ))
)
ggsave("Fig5/Macrophage_Tcell_lt-net.pdf", p2, width = 7, height = 5)

gge <- rownames(vis_ligand_target)
gge <- gsub("[.]", "-", gge)
t_macro <- macro
t_macro@active.ident <- factor(paste(t_macro$celltype, t_macro$subtype, t_macro$drug, sep = "_")) # nolint
names(t_macro@active.ident) <- rownames(t_macro@meta.data)
sub <- subset(t_macro, ident = c(
    "SPP1+ Mac_OA_BT_None", "SPP1+ Mac_RA_BT_Adalimumab",
    "SPP1+ Mac_RA_AT_Adalimumab", "SPP1+ Mac_RA_BT_Tofacitinib",
    "SPP1+ Mac_RA_AT_Tofacitinib", "S100A12+ Mac_OA_BT_None",
    "S100A12+ Mac_RA_BT_Adalimumab", "S100A12+ Mac_RA_AT_Adalimumab",
    "S100A12+ Mac_RA_BT_Tofacitinib", "S100A12+ Mac_RA_AT_Tofacitinib"
))
sub$kk <- sub@active.ident
mat <- sub@assays$RNA@data[gge, ]
mt <- aggregate(t(mat), list(sub$kk), mean)
rownames(mt) <- mt$Group.1
mt <- mt[, -1]
tt <- as.data.frame(t(mt))
tt <- tt[, c(6, 9, 7, 10, 8, 1, 4, 2, 5, 3)]

bk <- c(seq(-1.5, -0.01, by = 0.01), seq(0, 1.5, by = 0.01))
colo2 <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)
p1 <- pheatmap(tt,
    scale = "row",
    cluster_cols = FALSE, cluster_rows = FALSE, treeheight_row = 0,
    treeheight_col = 0, border_color = "black",
    color = rev(c(
        colorRampPalette(colors = colo2[100:51])(length(bk) / 2),
        colorRampPalette(colors = colo2[50:1])(length(bk) / 2)
    )),
    breaks = bk, gaps_col = 5
)
ggsave("Fig5/ligand_exp.pdf", p1, width = 2.5, height = 6)

ggt <- colnames(vis_ligand_target)
ggt <- gsub("[.]", "-", ggt)
t_tcell <- tcell
t_tcell@active.ident <- factor(paste(t_tcell$celltype, t_tcell$subtype, t_tcell$drug, sep = "_")) # nolint
names(t_tcell@active.ident) <- rownames(t_tcell@meta.data)
sub <- subset(t_tcell, ident = c(
    "GZMK+ CD8+ Tem1_OA_BT_None", "GZMK+ CD8+ Tem1_RA_BT_Adalimumab",
    "GZMK+ CD8+ Tem1_RA_AT_Adalimumab", "GZMK+ CD8+ Tem1_RA_BT_Tofacitinib",
    "GZMK+ CD8+ Tem1_RA_AT_Tofacitinib", "Treg_OA_BT_None",
    "Treg_RA_BT_Adalimumab", "Treg_RA_AT_Adalimumab",
    "Treg_RA_BT_Tofacitinib", "Treg_RA_AT_Tofacitinib",
    "CD4+ Tex_OA_BT_None", "CD4+ Tex_RA_BT_Adalimumab",
    "CD4+ Tex_RA_AT_Adalimumab", "CD4+ Tex_RA_BT_Tofacitinib",
    "CD4+ Tex_RA_AT_Tofacitinib"
))
sub$kk <- sub@active.ident
mat <- sub@assays$RNA@data[ggt, ]
mt <- aggregate(t(mat), list(sub$kk), mean)
rownames(mt) <- mt$Group.1
mt <- mt[, -1]
tt <- as.data.frame(t(mt))
tt <- tt[, c(6, 9, 7, 10, 8, 11, 14, 12, 15, 13, 1, 4, 2, 5, 3)]

bk <- c(seq(-1.5, -0.01, by = 0.01), seq(0, 1.5, by = 0.01))
colo2 <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)
p1 <- pheatmap(t(tt),
    scale = "column",
    cluster_cols = FALSE, cluster_rows = FALSE, treeheight_row = 0,
    treeheight_col = 0, border_color = "black",
    color = rev(c(
        colorRampPalette(colors = colo2[100:51])(length(bk) / 2),
        colorRampPalette(colors = colo2[50:1])(length(bk) / 2)
    )),
    breaks = bk, gaps_row = c(5, 10)
)
ggsave("Fig5/receiver_exp.pdf", p1, width = 6, height = 2.1)
#############################################################


##################### italk analysis ########################
## p1
mat <- data@meta.data[, c("subcell", "subtype", "drug")]
mat <- mat[which(mat$subcell %in% c("SPP1+ Mac", "S100A12+ Mac", "Treg", "GZMK+ CD8+ Tem1", "CD4+ Tex")), ] # nolint
mat$subcell <- as.character(mat$subcell)
mat$subtype <- as.character(mat$subtype)
mat$group <- paste(mat$subcell, mat$subtype, sep = "_")

iTalk_data <- as.data.frame(t(data@assays$RNA@counts[, rownames(mat)]))
iTalk_data$cell_type <- mat$group
iTalk_data$compare_group <- mat$subtype
colo <- c(pal_npg(alpha = 0.9)(9), pal_aaas(alpha = 0.7)(9), pal_lancet(alpha = 0.5)(9)) # nolint
my10colors <- my36colors <- colo
highly_exprs_genes <- rawParse(iTalk_data, top_genes = 50, stats = "mean")

comm_list <- c("cytokine")
cell_types <- unique(iTalk_data$cell_type)
cell_col <- structure(
    rep(c("#C6B4D3", "#F1B26E", "#E64B35E5", "#91D1C2E5", "#00A087E5"), each = 3), # nolint
    names = cell_types[c(1, 6, 11, 2, 7, 15, 3, 10, 14, 5, 8, 12, 4, 9, 13)]
)

iTalk_res <- NULL
for (comm_type in comm_list) {
    res_cat <- FindLR(highly_exprs_genes, datatype = "mean count", comm_type = comm_type) # nolint
    iTalk_res <- rbind(iTalk_res, res_cat)
}
iTalk_res <- iTalk_res[order(iTalk_res$cell_from_mean_exprs * iTalk_res$cell_to_mean_exprs, decreasing = T), ] # nolint

iTalk_res$group_from <- paste(sapply(strsplit(iTalk_res$cell_from, "_", fixed = T), "[", 2), sapply(strsplit(iTalk_res$cell_from, "_", fixed = T), "[", 3), sep = "_") # nolint
iTalk_res$group_to <- paste(sapply(strsplit(iTalk_res$cell_to, "_", fixed = T), "[", 2), sapply(strsplit(iTalk_res$cell_to, "_", fixed = T), "[", 3), sep = "_") # nolint
iTalk_res <- iTalk_res[which(iTalk_res$group_from == iTalk_res$group_to), ]
iTalk_res$ex_from <- sapply(strsplit(iTalk_res$cell_from, "_", fixed = T), "[", 1) # nolint
iTalk_res$ex_to <- sapply(strsplit(iTalk_res$cell_to, "_", fixed = T), "[", 1)
saveres <- iTalk_res

iTalk_res <- iTalk_res[order(iTalk_res$cell_from_mean_exprs * iTalk_res$cell_to_mean_exprs, decreasing = T), ][1:120, ] # nolint
iTalk_res$color <- ifelse(iTalk_res$ex_from == "SPP1+ Mac", as.character(cell_col[1]), # nolint
    ifelse(iTalk_res$ex_from == "S100A12+ Mac", as.character(cell_col[4]), # nolint
        ifelse(iTalk_res$ex_from == "GZMK+ CD8+ Tem1", as.character(cell_col[7]), # nolint
            ifelse(iTalk_res$ex_from == "Treg", as.character(cell_col[10]), as.character(cell_col[13])) # nolint
        )
    )
)

iTalk_res$cell_from <- gsub("RA_BT", "RA_0BT", iTalk_res$cell_from)
iTalk_res$cell_to <- gsub("RA_BT", "RA_0BT", iTalk_res$cell_to)
names(cell_col) <- gsub("RA_BT", "RA_0BT", names(cell_col))

pdf("lpplot_cytokine-group2.pdf", width = 6, height = 6)
LRPlot(iTalk_res[1:120, ], datatype = "mean count", cell_col = cell_col, link.arr.lwd = iTalk_res$cell_from_mean_exprs[1:120], link.arr.width = iTalk_res$cell_to_mean_exprs[1:120], link.arr.col = iTalk_res$color) # nolint
dev.off()

## p2
group_color <- c("#87B1C8", "#A1568E")
mat <- data@meta.data[, c("subcell", "ACR20")]
mat <- mat[which(mat$subcell %in% c("SPP1+ Mac", "S100A12+ Mac", "Treg", "GZMK+ CD8+ Tem1", "CD4+ Tex")), ] # nolint
mat <- subset(mat, ACR20 %in% c("N", "Y"))
mat$subcell <- as.character(mat$subcell)
mat$group <- paste(mat$subcell, mat$ACR20, sep = "_")

iTalk_data <- as.data.frame(t(data@assays$RNA@counts[, rownames(mat)]))
iTalk_data$cell_type <- mat$group
iTalk_data$compare_group <- mat$ACR20
my10colors <- my36colors <- colo
highly_exprs_genes <- rawParse(iTalk_data, top_genes = 50, stats = "mean")

comm_list <- c("cytokine")
cell_types <- unique(iTalk_data$cell_type)
cell_col <- structure(
    rep(c("#C6B4D3", "#F1B26E", "#E64B35E5", "#91D1C2E5", "#00A087E5"), each = 2), # nolint
    names = cell_types[c(1, 7, 2, 6, 5, 9, 3, 10, 4, 8)]
)

iTalk_res <- NULL
for (comm_type in comm_list) {
    res_cat <- FindLR(highly_exprs_genes, datatype = "mean count", comm_type = comm_type) # nolint
    iTalk_res <- rbind(iTalk_res, res_cat)
}
iTalk_res <- iTalk_res[order(iTalk_res$cell_from_mean_exprs * iTalk_res$cell_to_mean_exprs, decreasing = T), ] # nolint

iTalk_res$group_from <- sapply(strsplit(iTalk_res$cell_from, "_", fixed = T), "[", 2) # nolint
iTalk_res$group_to <- sapply(strsplit(iTalk_res$cell_to, "_", fixed = T), "[", 2) # nolint
iTalk_res <- iTalk_res[which(iTalk_res$group_from == iTalk_res$group_to), ]
iTalk_res$ex_from <- sapply(strsplit(iTalk_res$cell_from, "_", fixed = T), "[", 1) # nolint
iTalk_res$ex_to <- sapply(strsplit(iTalk_res$cell_to, "_", fixed = T), "[", 1)
saveres <- iTalk_res

iTalk_res <- iTalk_res[order(iTalk_res$cell_from_mean_exprs * iTalk_res$cell_to_mean_exprs, decreasing = T), ][1:120, ] # nolint
iTalk_res$color <- ifelse(iTalk_res$ex_from == "SPP1+ Mac", as.character(cell_col[1]), # nolint
    ifelse(iTalk_res$ex_from == "S100A12+ Mac", as.character(cell_col[4]), # nolint
        ifelse(iTalk_res$ex_from == "GZMK+ CD8+ Tem1", as.character(cell_col[7]), # nolint
            ifelse(iTalk_res$ex_from == "Treg", as.character(cell_col[10]), as.character(cell_col[13])) # nolint
        )
    )
)

iTalk_res$cell_from <- gsub("RA_BT", "RA_0BT", iTalk_res$cell_from)
iTalk_res$cell_to <- gsub("RA_BT", "RA_0BT", iTalk_res$cell_to)
names(cell_col) <- gsub("RA_BT", "RA_0BT", names(cell_col))

pdf("lpplot_cytokine-group2.pdf", width = 6, height = 6)
LRPlot(iTalk_res[1:120, ], datatype = "mean count", cell_col = cell_col, link.arr.lwd = iTalk_res$cell_from_mean_exprs[1:120], link.arr.width = iTalk_res$cell_to_mean_exprs[1:120], link.arr.col = iTalk_res$color) # nolint
dev.off()
#############################################################
