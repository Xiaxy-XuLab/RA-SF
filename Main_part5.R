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
data <- readRDS("data_rds/tcell.Rds")
#############################################################


###################### Tsne plot ############################
p <- DimPlot(data, label = F, reduction = "tsne") + NoLegend() +
    scale_color_manual(
        values = celltype_color
    ) +
    theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank()
    )
ggsave("Fig4/a_umap.tif", p, width = 4, height = 4, device = "tiff")
#############################################################


################# Percentage analysis #######################
mat_data <- melt(table(data@meta.data[, c("sampleid", "celltype")]) / as.vector(table(data_input$sampleid))) # nolint
mat_data <- mat_data %>%
    mutate(group = factor(
        rep(
            c(rep("OA-BT", 3), rep(c("RA-BT", "RA-AT"), 6)),
            length(unique(mat_data$celltype))
        ),
        levels = c("OA-BT", "RA-BT", "RA-AT")
    )) %>%
    mutate(gg = ifelse(group == "OA-BT", "OA", "RA")) %>%
    mutate(patient = gsub("_[0-9]$", "", sampleid)) %>%
    mutate(treat = sapply(strsplit(as.character(group), "-", fixed = T), "[", 2)) %>% # nolint
    mutate(patient = rep(c(1:3, rep(4:9, each = 2)), length(unique(mat_data$celltype)))) %>% # nolint
    mutate(drug = factor(rep(
        c(rep("None", 3), rep(c("Adalimumab", "Tofacinitib"), each = 6)),
        length(unique(mat_data$celltype))
    ), levels = c("None", "Adalimumab", "Tofacinitib")))

type <- "CD4+ Tex"
mt_data <- mat_data[which(mat_data$celltype == type), ]
p <- ggplot(mt_data, aes(group, value)) +
    stat_summary(fun.y = mean, geom = "bar", color = "black", fill = "white", width = .6) + # nolint
    stat_summary(fun.data = mean_se, geom = "errorbar", color = "black", width = .2) + # nolint
    geom_point(aes(color = group, shape = drug), size = 5) +
    geom_line(aes(group = patient), linewidth = 0.6, colour = "#9C9C9C") +
    scale_shape_manual(values = c(15, 1, 19)) +
    scale_color_manual(values = group_color3) +
    theme_classic() +
    stat_compare_means(comparisons = list(c("OA-BT", "RA-BT"), c("OA-BT", "RA-AT"), c("RA-BT", "RA-AT"))) + # nolint
    scale_y_continuous(expand = c(0, 0)) + # nolint
    labs(y = "Proportion of CXCL13+ CD4+ T", title = "") +
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
ggsave("Fig4/e_percent.pdf", p, width = 2.4, height = 4)

## p2
decon <- function(x = bulk) {
    bseqsc_config("/work/xiaxy/CIBERSORT.R")
    count <- as.matrix(data@assays$RNA@counts)
    celltype <- data.frame(
        cell = rownames(data@meta.data),
        celltype = as.character(data$celltype)
    )
    bb <- bseqsc_basis(
        count, tmarkerlist,
        clusters = celltype$celltype,
        samples = data$orig.ident, ct.scale = TRUE
    )
    fit <- bseqsc_proportions(x, bb, verbose = TRUE)
    return(melt(as.matrix(fit$coefficients)))
}

bulk <- read.table("input/Bulk_abundance.txt", header = T, sep = "\t", row.names = 1) # nolint
mt_data <- decon(bulk)
mt_data <- mt_data %>%
    arrange(Var1) %>%
    mutate(group = factor(rep(
        c(rep("OA-BT", 5), rep("RA-BT", 14), rep("RA-AT", 10)),
        length(unique(Var1))
    ), levels = c("OA-BT", "RA-BT", "RA-AT")))

tt_data <- subset(mt_data, Var1 == "CD4+ Tex")
bulk_clinical <- read.table("input/Bulk_clinical.txt", header = T, sep = "\t")
ttb <- merge(tt_data, bulk_clinical, by.x = "Var2", by.y = "sample")
ttb$patient <- paste(
    sapply(strsplit(as.character(ttb$Var2), "_", fixed = T), "[", 1),
    sapply(strsplit(as.character(ttb$Var2), "_", fixed = T), "[", 3),
    sep = "_"
)
ttb$drug <- factor(ttb$drug, levels = c("None", "Adalimumab", "Tofacitinib"))

p <- ggplot(ttb, aes(group, value)) +
    stat_summary(fun.y = mean, geom = "bar", color = "black", fill = "white", width = .6) + # nolint
    stat_summary(fun.data = mean_se, geom = "errorbar", color = "black", width = .2) + # nolint
    geom_point(aes(color = group, shape = drug), size = 5) +
    geom_line(aes(group = patient), linewidth = 0.6, colour = "#9C9C9C") +
    scale_shape_manual(values = c(15, 1, 19)) +
    scale_color_manual(values = group_color3) +
    theme_classic() +
    stat_compare_means(comparisons = list(c("OA-BT", "RA-BT"), c("OA-BT", "RA-AT"), c("RA-BT", "RA-AT"))) + # nolint
    # scale_y_continuous(expand = c(0, 0), limits = c(0, 0.25), breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25), labels = c("0.0%", "5%", "10%", "15%", "20%", "25%")) + # nolint
    labs(y = "Proportion of CXCL13+ CD4+ T", title = "") +
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
p
ggsave("Fig4/k_percent.pdf", p, width = 2.4, height = 4)

## p3
bulk <- read.table("input/GSE12021.txt", header = T, sep = "\t", row.names = 1) # nolint
mt_data <- decon(bulk)
mt_data <- mt_data[order(mt_data$Var1), ]
mt_data$Var1 <- factor(mt_data$Var1, levels = levels(mt_data$Var1))
mt_data$group <- c(
    rep("Normal", 4), "RA", rep("OA", 2),
    rep("RA", 5), "OA", "RA", "RA", rep("OA", 6),
    "RA", "RA", "RA", "OA", "RA", rep("Normal", 5)
)
mt_data <- subset(mt_data, group %in% c("OA", "RA"))
mt_data$group <- factor(mt_data$group, levels = c("OA", "RA"))
mt_data <- subset(mt_data, Var1 == "CD4+ Tex")

p <- ggplot(mt_data, aes(x = group, y = value)) +
    geom_boxplot(aes(color = group), outlier.size = 0, outlier.colour = "white") + # nolint
    geom_point(aes(fill = group),
        position = position_jitterdodge(jitter.width = 0.4),
        color = "black", shape = 21, size = 5
    ) +
    stat_compare_means(
        comparisons =
            list(c("OA", "RA")),
        label.y = c(0.38)
    ) +
    scale_color_manual(values = c("#87B1C8", "#A1568E")) +
    scale_fill_manual(values = c("#87B1C8", "#A1568E")) +
    theme_classic2() +
    labs(y = "Signature score", x = "") +
    scale_y_continuous(
        expand = c(0, 0), limits = c(0, 0.42)
    ) +
    theme(
        axis.text.x = element_text(color = "black", size = 14, angle = 45, hjust = 1), # nolint
        axis.text.y = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 16),
        legend.title = element_blank(), legend.position = "none"
    )
p
ggsave("Fig4/GSE12021.pdf", p, width = 2, height = 3.5)

## p4
bulk <- read.table("input/PMID35602512.txt", header = T, sep = "\t", row.names = 1) # nolint
mt_data <- decon(bulk)
mt_data <- mt_data[order(mt_data$Var1), ]
mt_data$Var1 <- factor(mt_data$Var1, levels = levels(mt_data$Var1))
mt_data$group <- c(rep("OA", 15), rep("RA", 9))
mt_data$group <- factor(mt_data$group, levels = c("OA", "RA"))

p <- ggplot(mt_data, aes(x = group, y = value)) +
    geom_boxplot(aes(color = group), outlier.size = 0, outlier.colour = "white") + # nolint
    geom_point(aes(fill = group),
        position = position_jitterdodge(jitter.width = 0.4),
        color = "black", shape = 21, size = 5
    ) +
    stat_compare_means(
        comparisons =
            list(c("OA", "RA")),
        label.y = c(0.38)
    ) +
    scale_color_manual(values = c("#87B1C8", "#A1568E")) +
    scale_fill_manual(values = c("#87B1C8", "#A1568E")) +
    theme_classic2() +
    labs(y = "Signature score", x = "") +
    scale_y_continuous(
        expand = c(0, 0), limits = c(0, 0.42)
    ) +
    theme(
        axis.text.x = element_text(color = "black", size = 14, angle = 45, hjust = 1), # nolint
        axis.text.y = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 16),
        legend.title = element_blank(), legend.position = "none"
    )
p
ggsave("Fig4/PMID35602512.pdf", p, width = 2, height = 3.5)

## p5
bulk <- read.table("input/Bulk_abundance.txt", header = T, sep = "\t", row.names = 1) # nolint
patients <- colnames(bulk)[c(6:12, 14:19)]
bulk <- bulk[, patients]
mt_data <- decon(bulk)
mt_data <- mt_data[order(mt_data$Var1), ]
mt_data$Var1 <- factor(mt_data$Var1, levels = levels(mt_data$Var1))
mt_data$group <- factor(c(
    rep("ACR20_N", 3), rep("ACR20_Y", 2),
    "ACR20_N", "ACR20_Y", "ACR20_N", "ACR20_Y",
    rep("ACR20_N", 3), "ACR20_Y"
))

p <- ggplot(mt_data, aes(x = Var1, y = value)) +
    geom_boxplot(aes(color = group), outlier.size = 0, outlier.colour = "white") + # nolint
    geom_point(aes(fill = group),
        position = position_jitterdodge(jitter.width = 0.4),
        color = "black", shape = 21, size = 3
    ) +
    scale_color_manual(values = c("#87B1C8", "#A1568E")) +
    scale_fill_manual(values = c("#87B1C8", "#A1568E")) +
    theme_classic2() +
    labs(y = "Deconvolution score", x = "") +
    scale_y_continuous(
        expand = c(0, 0.01), limits = c(0, 1)
    ) +
    theme(
        axis.text.x = element_text(color = "black", size = 14, angle = 45, hjust = 1), # nolint
        axis.text.y = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 16),
        legend.position = c(0.65, 0.9), legend.title = element_blank()
    )
p
ggsave("Fig4/k_bulk-decon.pdf", p, width = 5.8, height = 4)

## p6
mat <- melt(table(data@meta.data[, c("sampleid", "celltype")]) / as.vector(table(data_input$sampleid))) # nolint
mat <- mat %>%
    mutate(group = factor(rep(
        c(rep("OA-BT", 3), rep(c("RA-BT", "RA-AT"), 6)),
        length(unique(data$celltype))
    ), levels = c("OA-BT", "RA-BT", "RA-AT"))) %>%
    mutate(gg = ifelse(group == "OA-BT", "OA", "RA")) %>%
    mutate(patient = gsub("_[0-9]$", "", sampleid)) %>%
    mutate(treat = sapply(strsplit(as.character(group), "-", fixed = T), "[", 2)) %>% # nolint
    mutate(drug = factor(rep(
        c(rep("None", 3), rep(c("Adalimumab", "Tofacinitib"), each = 6)),
        length(unique(data$celltype))
    ), levels = c("None", "Adalimumab", "Tofacinitib")))

type <- "GZMK+ CD8+ Tem1"
ma <- mat[which(mat$celltype == type), ]
ma_cli <- cbind(ma, clinical_index_single)
ma_cli <- ma_cli[-(1:3), 1:12]
ma_cli$value <- 100 * ma_cli$value

p <- ggplot(ma_cli, aes(x = das28_index, y = value)) +
    geom_point(aes(shape = drug, color = group), size = 3) +
    scale_shape_manual(values = c(1, 19)) +
    scale_color_manual(values = group_color3[2:3]) +
    theme_classic() +
    geom_smooth(method = "lm", se = F, color = "black") +
    labs(y = "Proportion of CD8+ Tem1") +
    stat_cor() +
    theme(
        axis.text = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 10, color = "black"),
        legend.title = element_blank(),
        legend.position = "top"
    )
p
ggsave("Fig4/f_correlation1.pdf", p, width = 3, height = 3.3)
#############################################################


################# Percentage analysis #######################
## p1
sub <- subset(data, idents = "CD4+ Tex")
sub$gg <- paste(sub$subtype, sub$drug, sep = "_")
sub$gg <- factor(sub$gg, levels = unique(sub$gg))
exp <- sub@assays$RNA@data
mat <- data.frame(t(exp[c("CXCL13", "LAG3", "PDCD1", "CTLA4"), ]), group = sub$gg) # nolint
rownames(mat) <- NULL
mt <- melt(mat)

p <- ggplot(mt, aes(x = variable, y = value, fill = group)) +
    geom_boxplot(outlier.size = 0, outlier.alpha = 0) +
    labs(y = "Expression level in CXCL13+ CD4+ T") +
    stat_compare_means(comparisons = list(
        c("OA_BT_None", "RA_BT_Adalimumab"), c("OA_BT_None", "RA_BT_Tofacitinib"), # nolint
        c("RA_BT_Adalimumab", "RA_AT_Adalimumab"), c("RA_BT_Tofacitinib", "RA_AT_Tofacitinib") # nolint
    )) +
    scale_fill_manual(values = group_color5) + # nolint
    theme_classic() +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 10)) +
    theme(
        axis.title.x = element_blank(),
        legend.title = element_blank()
    )
p
ggsave("Fig4/Dys_marker_expression.pdf", p, width = 5.5, height = 2)

## p2
bulk <- read.table("input/Bulk_abundance.txt", header = T, sep = "\t", row.names = 1) # nolint
genelist <- c("CXCL13", "CTLA4", "LAG3", "PDCD1", "TIGIT", "HAVCR2")

data_set <- list(enriched = genelist)
result <- gsva(
    as.matrix(bulk), data_set,
    min.sz = 5,
    kcdf = "Poisson", method = "ssgsea",
    mx.diff = TRUE, verbose = FALSE, parallel.sz = 30
)
mat <- data.frame(name = colnames(result), value = as.numeric(result))
mat$group <- factor(c(rep("OA-BT", 5), rep("RA-BT", 14), rep("RA-AT", 10)))
bulk_clinical <- read.table("input/Bulk_clinical.txt", header = T, sep = "\t")
ttb <- merge(mat, bulk_clinical, by.x = "name", by.y = "sample")
ttb$patient <- paste(
    sapply(strsplit(as.character(ttb$name), "_", fixed = T), "[", 1),
    sapply(strsplit(as.character(ttb$name), "_", fixed = T), "[", 3),
    sep = "_"
)
ttb$drug <- factor(ttb$drug, levels = c("None", "Adalimumab", "Tofacitinib"))
ttb$drug <- factor(ttb$drug, levels = c("Adalimumab", "Tofacitinib"))
ttb$group <- factor(ttb$group, levels = c("OA-BT", "RA-BT", "RA-AT"))

p <- ggplot(ttb, aes(group, value)) +
    stat_summary(fun.y = mean, geom = "bar", color = "black", fill = "white", width = .6) + # nolint
    stat_summary(fun.data = mean_se, geom = "errorbar", color = "black", width = .2) + # nolint
    geom_point(aes(color = group, shape = drug), size = 5) +
    geom_line(aes(group = patient), size = 0.6, colour = "#9C9C9C") +
    scale_shape_manual(values = c(15, 1, 19)) +
    scale_color_manual(values = group_color) +
    theme_classic() +
    stat_compare_means(comparisons = list(c("OA-BT", "RA-BT"), c("OA-BT", "RA-AT"), c("RA-BT", "RA-AT"))) + # nolint
    scale_y_continuous(expand = c(0, 0), limits = c(0, 2)) + # nolint
    labs(y = "Dysfunction scores", x = "", title = "") +
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
p
ggsave("Fig4/dys.pdf", p, width = 2.8, height = 4)

ss <- subset(ttb, patient %in% ttb$patient[which(duplicated(ttb$patient))])
val <- data.frame(drug = ss$drug[2:10], value = log(ss$value[2:10] / ss$value[12:20], 2)) # nolint
p <- ggplot(val, aes(x = drug, y = value)) +
    geom_boxplot() +
    geom_point(aes(shape = drug),
        position = position_auto(scale = FALSE), size = 3
    ) +
    scale_shape_manual(values = c(1, 19)) +
    stat_compare_means(comparisons = list(c("Adalimumab", "Tofacitinib"))) +
    theme_classic() +
    labs(x = "", y = "Log2 FC (Dysfunction scores)") +
    theme(
        axis.text.x = element_text(
            angle = 30, hjust = 1, color = "black", size = 10
        ),
        axis.text.y = element_text(color = "black", size = 10),
        axis.title.y = element_text(color = "black", size = 12),
        legend.position = "none"
    )
p
ggsave("Fig4/dys1.pdf", p, height = 4, width = 2.8)

vp_case1 <- function(cell, gene_signature, file_name, ncol, cc) {
    sub <- subset(data_input, cells = cell)
    sub <- subset(sub, cells = rownames(sub@meta.data)[which(sub$subcell == cc)]) # nolint
    sub@active.ident <- factor(sub$ACR20)

    plot_case1 <- function(signature) {
        VlnPlot(sub,
            features = signature,
            pt.size = 0,
            y.max = 1.2 * max(as.numeric(sub@assays$RNA@data[signature, ])),
            cols = c("#87B1C8", "#A1568E")
        ) +
            geom_boxplot(width = .1, col = "black", fill = "white", outlier.size = 0) + # nolint
            stat_compare_means(
                comparisons = list(levels(sub)),
                label = "p.signif",
                label.y = 1.1 * max(as.numeric(sub@assays$RNA@data[signature, ])) # nolint
            ) +
            labs(y = "", x = "") +
            scale_x_discrete(labels = c("ACR20_N", "ACR20_Y")) &
            theme(
                legend.position = "none",
                axis.text.x = element_text(angle = 45, color = "black", hjust = 1, size = 11), # nolint
                axis.text.y = element_text(color = "black", size = 14)
            )
    }

    purrr::map(gene_signature, plot_case1) %>% cowplot::plot_grid(plotlist = ., ncol = ncol) # nolint
    file_name <- paste0(file_name, ".pdf")
    ggsave(file_name, width = 2 * ncol, height = 3 * length(gene_signature) / ncol) # nolint
}
genes <- c("CXCL13", "CTLA4", "LAG3", "PDCD1")
vp_case1(
    cell = rownames(data_input@meta.data)[which(data_input$ACR20 %in% c("Y", "N"))], # nolint
    gene_signature = genes, file_name = "Fig4/ACR20_Tex", ncol = 4, cc = "CD4+ Tex" # nolint
)
#############################################################
