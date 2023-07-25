#!/usr/bin/env R
# coding utf-8

###################### Package Loading #######################
pacman::p_load(
    "Seurat", "ggplot2", "ggsci", "reshape2", "scales",
    "monocle", "viridis", "patchwork", "ggpubr", "ggrepel",
    "msigdbr", "fgsea", "tidyverse", "bseqsc", "ggsignif",
    "ggrepel", "pheatmap", "RColorBrewer", "ggpmisc", "ggplotify",
    "circlize", "zoo", "nichenetr", "tidyverse", "GSVA", "magrittr",
    "ggforce", "dplyr", "ggpointdensity", "Nebulosa"
)
setwd("/work/xiaxy/work/RA/NC")
source("/work/xiaxy/work/RA/NC/Dependent.R")
#############################################################


###################### Read input ###########################
data_input <- readRDS("data_rds/input.Rds")
data <- readRDS("data_rds/macrophage.Rds")
#############################################################


###################### Tsne plot ############################
p <- DimPlot(data, label = F, reduction = "tsne") + NoLegend() +
    scale_color_manual(
        values = macrophage_color
    ) +
    theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank()
    )
ggsave("Fig3/a_umap.tif", p, width = 3, height = 3, device = "tiff")
#############################################################


################# Percentage analysis #######################
## p1
data_matrix <- melt(table(data@meta.data[, c("sampleid", "celltype")]) / as.vector(table(data_input$sampleid))) # nolint
data_matrix <- data_matrix %>%
    mutate(group = factor(
        rep(
            c(rep("OA-BT", 3), rep(c("RA-BT", "RA-AT"), 6)),
            length(unique(data$celltype))
        ),
        levels = c("OA-BT", "RA-BT", "RA-AT")
    )) %>%
    mutate(new_group = ifelse(group == "OA-BT", "OA", "RA")) %>%
    mutate(patient = gsub("_[0-9]$", "", sampleid)) %>%
    mutate(treat = sapply(strsplit(as.character(group), "-", fixed = T), "[", 2)) %>% # nolint
    mutate(
        drug = factor(
            rep(
                c(rep("None", 3), rep(c("Adalimumab", "Tofacinitib"), each = 6)), # nolint
                length(unique(data$celltype))
            ),
            levels = c("None", "Adalimumab", "Tofacinitib")
        )
    )

percentage_graph <- function(x) {
    data_graph <- subset(data_matrix, celltype == x)
    p <- ggplot(data_graph, aes(group, value)) +
        stat_summary(fun.y = mean, geom = "bar", color = "black", fill = "white", width = .6) + # nolint
        stat_summary(fun.data = mean_se, geom = "errorbar", color = "black", width = .2) + # nolint
        geom_point(aes(color = group, shape = drug), size = 5) +
        geom_line(aes(group = patient), linewidth = 0.6, colour = "#9C9C9C") +
        scale_shape_manual(values = c(15, 1, 19)) +
        scale_color_manual(values = group_color3) +
        theme_classic() +
        stat_compare_means(comparisons = list(c("OA-BT", "RA-BT"), c("OA-BT", "RA-AT"), c("RA-BT", "RA-AT"))) + # nolint
        scale_y_continuous(expand = c(0, 0), limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0.0%", "25%", "50%", "75%", "100%")) + # nolint
        labs(y = "Proportion", title = "") +
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

final <- percentage_graph("SPP1+ Mac") + percentage_graph("S100A12+ Mac") +
    plot_layout(nrow = 1, guides = "collect") & theme(legend.position = "top")
ggsave("Fig3/b_percent.pdf", final, width = 4.8, height = 4)

## p2
data_matrix <- melt(table(data@meta.data[, c("sampleid", "celltype")]) / as.vector(table(data$sampleid))) # nolint
spp1_proportion <- subset(data_matrix, celltype == "SPP1+ Mac")
mat <- data.frame(
    SPP1 = spp1_proportion[4:15, "value"],
    DAS28 = das28_index[4:15],
    group = factor(rep(c("RA-BT", "RA-AT"), each = 6), levels = c("RA-BT", "RA-AT")), # nolint
    drug = rep(c("Ada", "Ada", "Ada", "Tof", "Tof", "Tof"), time = 2)
)
p <- ggplot(mat, aes(x = SPP1, y = DAS28)) +
    geom_point(aes(color = group, shape = drug), size = 4.5) +
    stat_cor() +
    scale_color_manual(values = c("#C46F68", "#C4AA69")) +
    scale_shape_manual(values = c(1, 19)) +
    geom_smooth(se = FALSE, method = "lm", color = "black") +
    labs(x = "Proportion of SPP1+ Macro", y = "DAS28") +
    theme_classic() +
    theme(legend.title = element_blank())
ggsave("Fig3/SPP1_proportion_DAS28_cor.pdf", p, width = 3.7, height = 3)
#############################################################


#################### Elisa analysis #########################
## p1
data_matrix <- read.table("input/elisa_input.txt", header = T, sep = "\t")
data_matrix <- data_matrix %>%
    mutate(GROUP = factor(GROUP, levels = c("OA-BT", "RA-BT", "RA-AT"))) %>%
    mutate(DRUG = factor(DRUG, levels = c("None", "TNF", "JAK"))) %>%
    mutate(group_drug = factor(paste(GROUP, DRUG, sep = "_"),
        levels = c("OA-BT_None", "RA-BT_TNF", "RA-AT_TNF", "RA-BT_JAK", "RA-AT_JAK") # nolint
    )) %>%
    mutate(
        group_acr20 = factor(paste(GROUP, ACR20, sep = "_"),
            levels = c("OA-BT_NA", "RA-BT_NA", "RA-BT_N", "RA-AT_N", "RA-BT_Y", "RA-AT_Y") # nolint
        ) # nolint
    )
data_matrix_select <- rbind(
    data_matrix[1:10, ],
    data_matrix[which(data_matrix$LINE %in% data_matrix[which(duplicated(data_matrix$LINE)), "LINE"]), ] # nolint
)

elisa_graph <- function(x) {
    p <- ggplot(data_matrix, aes(GROUP, get(x))) +
        stat_summary(fun.y = mean, geom = "bar", color = "black", fill = "white", width = .7) + # nolint
        stat_summary(fun.data = mean_se, geom = "errorbar", color = "black", width = .2) + # nolint
        geom_point(aes(color = GROUP, shape = DRUG), size = 4) +
        geom_line(aes(group = LINE), size = 0.6, colour = "#9C9C9C") +
        stat_compare_means(
            comparisons = list(c("OA-BT", "RA-BT"), c("RA-BT", "RA-AT")),
            label.y = ifelse(x == "SPP1", 45, 2700), size = 3, bracket.size = 0.2 # nolint
        ) +
        scale_shape_manual(values = c(15, 1, 19)) +
        scale_color_manual(values = group_color3) +
        theme_classic() +
        scale_y_continuous(expand = c(0, 0), limits = c(0, ifelse(x == "SPP1", 50, 3000))) + # nolint
        labs(y = x) +
        theme(
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 12, color = "black"),
            plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(size = 10, color = "black", angle = 45, hjust = 1), # nolint
            axis.text.y = element_text(size = 10, color = "black"),
            axis.ticks.length = unit(0.4, "lines"),
            legend.title = element_blank(),
            legend.text = element_text(size = 14, color = "black")
        )
    return(p)
}
final <- elisa_graph("SPP1") + elisa_graph("CCL2") +
    plot_layout(nrow = 1, guides = "collect")
ggsave("Fig3/c_elisa.pdf", final, width = 5, height = 3)

## p2
data_SDAI <- read.table("input/elisa_before_SDAI.txt", header = T, sep = "\t")
p <- ggplot(data_SDAI, aes(x = SPP1, y = SDAI)) +
    geom_point(size = 3, color = "#D26964") +
    geom_smooth(se = FALSE, method = "lm", color = "black") +
    stat_cor(method = "pearson") +
    labs(x = "Expression of SPP1") +
    theme_classic() +
    theme(
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 12, color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 8, color = "black")
    )

data_DAS28 <- read.table("input/elisa_before_DAS28.txt", header = T, sep = "\t")
p1 <- ggplot(data_DAS28, aes(x = SPP1, y = DAS28)) +
    geom_point(size = 3, color = "#D26964") +
    geom_smooth(se = FALSE, method = "lm", color = "black") +
    stat_cor(method = "pearson") +
    labs(x = "Expression of SPP1") +
    theme_classic() +
    theme(
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 12, color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 8, color = "black")
    )
p1

## p3
data_matrix <- read.table("input/elisa_input.txt", header = T, sep = "\t")
sub_matrix <- subset(data_matrix, GROUP == "RA-BT") %>%
    slice(c(1:6, 8, 10, 12:15, 18:23, 29:30)) %>%
    mutate(ACR20 = factor(gsub("^", "ACR20_", ACR20), levels = c("ACR20_N", "ACR20_Y"))) %>% # nolint
    mutate(DRUG = factor(DRUG, levels = c("TNF", "JAK")))

p1 <- ggplot(sub_matrix, aes(ACR20, SPP1)) +
    stat_summary(fun.y = mean, geom = "bar", color = "black", fill = "white", width = .7) + # nolint
    stat_summary(fun.data = mean_se, geom = "errorbar", color = "black", width = .2) + # nolint
    geom_point(aes(color = ACR20, shape = DRUG), size = 4) +
    stat_compare_means(
        comparisons = list(c("ACR20_N", "ACR20_Y")),
        label.y = 42, size = 4, bracket.size = 0.2
    ) +
    geom_line(aes(group = LINE), size = 0.6, colour = "#9C9C9C") +
    scale_color_manual(values = c("#87B1C8", "#A1568E")) +
    scale_shape_manual(values = c(1, 19)) +
    theme_classic() +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 45)) +
    labs(ylab = "SPP1") +
    theme(
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12, color = "black"),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 10, color = "black"), # nolint
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.length = unit(0.4, "lines"),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 14, color = "black")
    )
p1

p2 <- ggplot(sub_matrix, aes(ACR20, CCL2)) +
    stat_summary(fun.y = mean, geom = "bar", color = "black", fill = "white", width = .7) + # nolint
    stat_summary(fun.data = mean_se, geom = "errorbar", color = "black", width = .2) + # nolint
    geom_point(aes(color = ACR20, shape = DRUG), size = 4) +
    stat_compare_means(
        comparisons = list(c("ACR20_N", "ACR20_Y")),
        label.y = 2100, size = 4, bracket.size = 0.2
    ) +
    geom_line(aes(group = LINE), size = 0.6, colour = "#9C9C9C") +
    scale_color_manual(values = c("#87B1C8", "#A1568E")) +
    scale_shape_manual(values = c(1, 19)) +
    theme_classic() +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 2400)) +
    labs(ylab = "CCL2") +
    theme(
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12, color = "black"),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 10, color = "black"), # nolint
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.length = unit(0.4, "lines"),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 14, color = "black")
    )
p2
final <- p1 + p2 + plot_layout(guides = "collect", nrow = 1)
ggsave("Fig3/elisa_acr20.pdf", final, width = 4.5, height = 2.5)
#############################################################


################# Expression analysis #######################
## p1
sub <- subset(data, cells = rownames(data@meta.data)[which(data$celltype == "SPP1+ Mac")]) # nolint
mat <- data.frame(
    SPP1 = sub@assays$RNA@data["SPP1", ],
    CCL2 = sub@assays$RNA@data["CCL2", ]
)
p <- ggplot(mat, aes(x = SPP1, y = CCL2)) +
    geom_pointdensity(adjust = 4, show.legend = TRUE, size = 0.8) +
    geom_smooth(aes(x = SPP1, y = CCL2), se = FALSE, method = "lm", color = "black") + # nolint
    stat_cor(aes(x = SPP1, y = CCL2), method = "pearson") +
    theme_classic() +
    scale_color_distiller(palette = "Spectral", direction = -1) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0))
p
ggsave("Fig3/SPP1_CCL2.pdf", p, width = 4, height = 3.2)

## p2
mat <- mat %>%
    mutate(group = factor(sub$subtype, levels = unique(sub$subtype))) %>%
    mutate(group_drug = factor(paste(group, sub$drug, sep = "_"),
        levels = c(
            "OA_BT_None", "RA_BT_Adalimumab",
            "RA_AT_Adalimumab", "RA_BT_Tofacitinib", "RA_AT_Tofacitinib"
        )
    ))

mt <- melt(mat)
p <- ggplot(mt, aes(x = variable, y = value)) +
    geom_boxplot(aes(fill = group), outlier.colour = "white", outlier.size = 0) + # nolint
    theme_classic() +
    scale_fill_manual(values = group_color3) +
    scale_y_continuous(limits = c(0, 9), expand = c(0, 0))
ggsave("Fig3/exp.pdf", p, width = 4, height = 3.5)

## p3
p1 <- ggplot(mat, aes(x = group_drug, y = SPP1, fill = group)) +
    geom_boxplot(outlier.colour = "white", outlier.size = 0) +
    theme_classic() +
    scale_fill_manual(values = group_color5) +
    stat_compare_means(
        comparisons = list(
            c("OA_BT_None", "RA_BT_Adalimumab"),
            c("OA_BT_None", "RA_BT_Tofacitinib"),
            c("RA_BT_Adalimumab", "RA_AT_Adalimumab"),
            c("RA_BT_Tofacitinib", "RA_AT_Tofacitinib")
        ), label.y = c(8, 8.7, 8, 8), label = "p.signif"
    ) + # nolint
    scale_y_continuous(expand = c(0, 0), limits = c(0, 10)) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        axis.title.x = element_blank(),
        legend.position = "none"
    )

p2 <- ggplot(mat, aes(x = group_drug, y = CCL2, fill = group)) +
    geom_boxplot(outlier.colour = "white", outlier.size = 0) +
    theme_classic() +
    scale_fill_manual(values = group_color5) +
    stat_compare_means(
        comparisons = list(
            c("OA_BT_None", "RA_BT_Adalimumab"),
            c("OA_BT_None", "RA_BT_Tofacitinib"),
            c("RA_BT_Adalimumab", "RA_AT_Adalimumab"),
            c("RA_BT_Tofacitinib", "RA_AT_Tofacitinib")
        ), label.y = c(6.5, 7, 6.5, 6.5), label = "p.signif"
    ) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 8)) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        axis.title.x = element_blank(),
        legend.position = "none"
    )
final <- p1 + p2 + plot_layout(nrow = 1)
ggsave("Fig3/SPP1_CCL2_group.pdf", final, width = 5, height = 3.5)

## p4
vp_case1 <- function(cell, gene_signature, file_name, ncol, cc) {
    sub_data <- subset(data_input, cells = cell)
    sub_data <- subset(sub_data, cells = rownames(sub_data@meta.data)[which(sub_data$subcell == cc)]) # nolint
    sub_data@active.ident <- factor(sub_data$ACR20)

    plot_case1 <- function(signature) {
        VlnPlot(sub_data,
            features = signature,
            pt.size = 0,
            y.max = 1.2 * max(as.numeric(sub_data@assays$RNA@data[signature, ])), # nolint
            cols = c("#87B1C8", "#A1568E")
        ) +
            geom_boxplot(width = .1, col = "black", fill = "white", outlier.size = 0) + # nolint
            stat_compare_means(
                comparisons = list(levels(sub_data)),
                label = "p.signif",
                label.y = 1.1 * max(as.numeric(sub_data@assays$RNA@data[signature, ])) # nolint
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

genes <- c(
    "FCGR3A", "HLA-A", "MT2A", "FOS", "STAT1",
    "GBP1", "IFITM3", "IFI16", "IRF1", "IFNGR1"
)

vp_case1(
    cell = rownames(data_input@meta.data)[which(data_input$ACR20 %in% c("Y", "N"))], # nolint
    gene_signature = genes, file_name = "Fig3/SPP1_macrophage_ACR20",
    ncol = 5, cc = "SPP1+ Mac"
)

genes <- c("HIF1A", "IFI30", "IL1B", "FOS", "IRF1", "JUN")
vp_case1(
    cell = rownames(data_input@meta.data)[which(data_input$ACR20 %in% c("Y", "N"))], # nolint
    gene_signature = genes, file_name = "Fig3/S100A12_macrophage_ACR20",
    ncol = 3, cc = "S100A12+ Mac"
)
#############################################################


############# Flow cytometry marker select ##################
p1 <- plot_density(data, c("SPP1", "CD36"), joint = TRUE)
p2 <- plot_density(data, c("S100A12", "SELL"), joint = TRUE)

select <- subset(data, cells = rownames(data@meta.data)[which(data$disease == "RA")]) # nolint
mat_data <- as.data.frame(t(select@assays$RNA@data[c("SPP1", "S100A12", "CD36", "SELL"), ])) # nolint
mat_data <- mat_data %>% mutate(sample = select$sampleid) %>%
    mutate(spp1_pos = ifelse(SPP1 > 0, "pos", "neg")) %>%
    mutate(s100a12_pos = ifelse(S100A12 > 0, "pos", "neg")) %>%
    mutate(cd36_pos = ifelse(CD36 > 0 & SELL == 0, "pos", "neg")) %>%
    mutate(cd62l_pos = ifelse(SELL > 0, "pos", "neg"))

pos <- c()
neg <- c()
for (i in unique(mat_data$sample)) {
    mt <- subset(mat_data, sample == i)
    value <- table(mt[, c("spp1_pos", "cd36_pos")])
    print(value)
    if (nrow(value) > 1) {
        neg <- append(neg, value[1, 2] / rowSums(value)[1])
        pos <- append(pos, value[2, 2] / rowSums(value)[2])
    } else {
        pos <- append(pos, value[1, 2] / rowSums(value)[1])
        neg <- append(neg, NA)
    }
}
mk <- data.frame(pos = c(neg, pos), group = rep(c("neg", "pos"), each = 12))

p <- ggboxplot(mk,
    x = "group", y = "pos",
    color = "group", palette = "jama",
    add = "jitter"
) +
    stat_compare_means(comparisons = list(c("neg", "pos"))) +
    theme(axis.title.x = element_blank())
p

pos <- c()
neg <- c()
for (i in unique(mat_data$sample)) {
    mt <- subset(mat_data, sample == i)
    value <- table(mt[, c("s100a12_pos", "cd62l_pos")])
    print(value)
    if (nrow(value) > 1) {
        neg <- append(neg, value[1, 2] / rowSums(value)[1])
        pos <- append(pos, value[2, 2] / rowSums(value)[2])
    } else {
        pos <- append(pos, value[1, 2] / rowSums(value)[1])
        neg <- append(neg, NA)
    }
}
mk <- data.frame(pos = c(neg, pos), group = rep(c("neg", "pos"), each = 12))
p1 <- ggboxplot(mk,
    x = "group", y = "pos",
    color = "group", palette = "jama",
    add = "jitter"
) +
    stat_compare_means(comparisons = list(c("neg", "pos"))) +
    theme(axis.title.x = element_blank())

p1
final <- p + p1 + plot_layout(nrow = 1, guides = "collect")
ggsave("Fig3/FC_marker_select.pdf", final, width = 6.7, height = 3)
#############################################################
