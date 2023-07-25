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
setwd("/work/xiaxy/work/RA/NC")
source("/work/xiaxy/work/RA/NC/Dependent.R")
#############################################################


###################### Read input ###########################
data_input <- readRDS("data_rds/input.Rds")
#############################################################


################# Famous genes analysis #####################
## p1
genes <- c(
    "IL6", "TNF", "IL1A", "IL1B", "IL12A",
    "IL12B", "IL18", "JAK1", "JAK2", "JAK3"
)
data_input@meta.data <- data_input@meta.data %>% mutate(
    merge_name = factor(paste(celltype, subtype, drug, sep = "_"),
        levels = paste(rep(levels(data_input$celltype), each = 5),
            rep(c("OA_BT", "RA_BT", "RA_AT", "RA_BT", "RA_AT"), 8),
            rep(c(
                "None", "Adalimumab", "Adalimumab", "Tofacitinib", "Tofacitinib"
            ), 8),
            sep = "_"
        )[-34]
    )
)
data_input@active.ident <- data_input$merge_name

pdf("Fig2/famous_genes_dotplot.pdf", width = 8, height = 3)
DotPlot(data_input, features = rev(genes), cols = c("#BDA9A9", "#D4382C")) +
    theme_bw() +
    theme(
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    coord_flip()
dev.off()

## p2
famous_draw <- function(x, y) {
    sub_data <- subset(data_input, cells = rownames(data_input@meta.data)[which(data_input$celltype == x)]) # nolint
    sub_data@meta.data <- sub_data@meta.data %>%
        mutate(new_group = factor(paste(subtype, drug, sep = "_"), levels = c(
            "OA_BT_None", "RA_BT_Adalimumab",
            "RA_AT_Adalimumab", "RA_BT_Tofacitinib", "RA_AT_Tofacitinib"
        )))
    sub_data@active.ident <- factor(sub_data$new_group)

    p <- VlnPlot(sub_data,
        features = y, pt.size = 0, ncol = length(y),
        cols = group_color5
    ) &
        geom_boxplot(
            width = .2, col = "black",
            fill = "white", outlier.size = 0
        ) &
        stat_compare_means(comparisons = list(
            c("OA_BT_None", "RA_BT_Adalimumab"),
            c("OA_BT_None", "RA_BT_Tofacitinib"),
            c("RA_BT_Adalimumab", "RA_AT_Adalimumab"),
            c("RA_BT_Tofacitinib", "RA_AT_Tofacitinib")
        )) &
        theme(
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            legend.position = "none"
        )
    return(p)
}

gene <- c("JAK1", "JAK3", "TNF")
final <- famous_draw("Macrophage", gene) / famous_draw("T cell", gene) +
    plot_layout(nrow = 2)
ggsave("Fig2/famous_genes_violin.pdf", final, width = 5, height = 4)
#############################################################


################### DEG number analysis #####################
deg_calculate <- function(data_set) {
    number_matrix <- matrix(NA, nrow = length(unique(data_input$celltype)), ncol = 2) # nolint
    k <- 1
    for (i in unique(data_set$celltype)) {
        sub_data <- subset(data_set, cells = rownames(data_set@meta.data)[which(data_set$celltype == i)]) # nolint
        number_matrix[k, 1] <- i
        if (table(data_set@active.ident)[1] > 10 & table(data_set@active.ident)[2] > 10) { # nolint
            marker <- FindAllMarkers(sub_data, only.pos = T, min.pct = 0.1, logfc.threshold = 0.1) # nolint
            if (nrow(marker) > 0) {
                marker_select <- subset(marker, p_val_adj < 0.05)
                if (nrow(marker_select) > 0) {
                    number_matrix[k, 2] <- nrow(marker_select)
                } else {
                    number_matrix[k, 2] <- NA
                }
            }
        }
        k <- k + 1
    }
    colnames(number_matrix) <- c("celltype", "number")
    number_matrix <- as.data.frame(number_matrix)
    data_mat <- as.data.frame(cbind(data_input@meta.data, data_input@reductions$tsne@cell.embeddings)) # nolint
    data_merge <- merge(data_mat, number_matrix, by.x = "celltype", by.y = "celltype") # nolint
    return(data_merge)
}

deg_graph <- function(graph_input) {
    p1 <- ggplot(graph_input, aes(x = tSNE_1, y = tSNE_2, color = as.numeric(number))) + # nolint
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
            legend.position = "bottom"
        ) +
        guides(color = guide_colourbar(title.position = "top", title.hjust = 1)) # nolint
    return(p1)
}

data_subset <- function(dat_input, column, select) {
    return(subset(dat_input, cells = rownames(subset(dat_input@meta.data, get(column) %in% select)))) # nolint
}

data_select1 <- data_subset(data_input, "subtype", c("OA_BT", "RA_BT"))
data_select1@active.ident <- factor(data_select1$subtype)
data_select2 <- data_subset(data_input, "subtype", c("RA_BT", "RA_AT"))
data_select2@active.ident <- factor(data_select2$subtype)
deg_select1 <- deg_calculate(data_select1)
deg_select2 <- deg_calculate(data_select2)

final <- deg_graph(deg_select1) + deg_graph(deg_select2) + plot_layout(nrow = 1)
ggsave("Fig2/l_deg.pdf", final, width = 5.5, height = 3)

data_select3 <- data_subset(data_subset(data_input, "subtype", c("RA_BT", "RA_AT")), "sampleid", unique(data_input$sampleid)[4:9]) # nolint
data_select3@active.ident <- factor(data_select3$subtype)
data_select4 <- data_subset(data_subset(data_input, "subtype", c("RA_BT", "RA_AT")), "sampleid", unique(data_input$sampleid)[10:15]) # nolint
data_select4@active.ident <- factor(data_select4$subtype)
final <- deg_graph(deg_select3) + deg_graph(deg_select4) + plot_layout(nrow = 1)
ggsave("Fig2/l_deg2.pdf", final, width = 5.5, height = 3)
#############################################################


###################### DEG analysis #########################
## p1
gene_subset <- function(x, y) {
    return(x[which(x[, "cluster"] == y & x[, "p_val_adj"] < 10^-50), "gene"])
}

pos_calculate <- function(x) {
    return(100 * sum(x > 0) / length(x))
}

vol_graph <- function(x, y, z, k, m) {
    x[, "avg_log2FC"] <- ifelse(x[, "cluster"] == y, -1 * x[, "avg_log2FC"], x[, "avg_log2FC"]) # nolint
    x[, "p_val_adj"] <- ifelse(x[, "p_val_adj"] == 0,
        min(x[which(x[, "p_val_adj"] != 0), "p_val_adj"]), x[, "p_val_adj"]
    )
    x[, "label"] <- ifelse(x[, "gene"] %in% z, x[, "gene"], NA)
    x[, "group"] <- ifelse(x[, "p_val_adj"] >= 10^-10, "N.S.",
        ifelse(x[, "gene"] %in% k & x[, "avg_log2FC"] > 0, "Up",
            ifelse(x[, "gene"] %in% k & x[, "avg_log2FC"] < 0, "Down", "Sig")
        )
    )
    x[, "group"] <- factor(x[, "group"], levels = c("N.S.", "Sig", "Down", "Up")) # nolint
    x <- x[order(x$group), ]

    p <- ggplot(x, aes(
        x = avg_log2FC, y = -log(p_val_adj, 10),
        color = group
    )) +
        geom_point(size = 1) +
        scale_color_manual(
            values = c(
                "lightgrey", "gray",
                "#77A270", "#E3B66B"
            )
        ) +
        scale_y_continuous(limits = c(0, 310)) +
        geom_text_repel(aes(label = label), color = "#3750A1", max.overlaps = 30, size = 2.5) + # nolint
        theme_classic() +
        labs(title = m, y = "", x = "") + # nolint
        theme(
            plot.title = element_text(hjust = 0.5, size = 12),
            legend.position = "none",
            axis.text = element_text(size = 10, color = "black"),
            axis.title = element_text(size = 12, color = "black")
        )
    return(p)
}

plotenrich <- function(y) {
    p <- plotEnrichment(fgsea_sets[[y]], ranks) +
        labs(title = y, x = "Rank", y = "Enrichment score") +
        theme(
            plot.title = element_text(hjust = 0.5, color = "black", size = 8),
            axis.text = element_text(color = "black")
        )
    return(p)
}

marker_a <- read.table("input/macrophage_deg_OA_BT-RA_BT.txt", header = T, sep = "\t") # nolint
marker_b <- read.table("input/macrophage_deg_RA_BT-RA_AT.txt", header = T, sep = "\t") # nolint
marker_c <- read.table("input/macrophage_deg_RA_BT-RA_AT_ada.txt", header = T, sep = "\t") # nolint
marker_d <- read.table("input/macrophage_deg_RA_BT-RA_AT_tof.txt", header = T, sep = "\t") # nolint

inter_neg <- Reduce(intersect, list(gene_subset(marker_a, "OA_BT"), gene_subset(marker_b, "RA_AT"), gene_subset(marker_c, "RA_AT"), gene_subset(marker_d, "RA_AT"))) # nolint
inter_pos <- Reduce(intersect, list(gene_subset(marker_a, "RA_BT"), gene_subset(marker_b, "RA_BT"), gene_subset(marker_c, "RA_BT"), gene_subset(marker_d, "RA_BT"))) # nolint

label <- c(
    "S100A4", "SPP1", "CCL2", "MMP19", "CXCL2", "CXCL3", "GBP1", "GBP5",
    "STAT1", "HLA-A", "HLA-B", "CXCL8", "SLAMF9", "NFKBIA", "CCL7"
)
label1 <- c(label, "S100A9", "CD44", "MIF", "CCR1")
label2 <- c(label, "JUN", "VEGFA", "ICAM1")
label3 <- c(label, "JUN", "MIF", "CD44", "VEGFA")
label4 <- c(label, "JUN", "IRF1", "MT2A")

final <- vol_graph(marker_a, "OA_BT", label1, c(inter_pos, inter_neg), "") +
    vol_graph(marker_b, "RA_AT", label2, c(inter_pos, inter_neg), "") +
    vol_graph(marker_c, "RA_AT", label3, c(inter_pos, inter_neg), "Adalimumab") + # nolint
    vol_graph(marker_d, "RA_AT", label4, c(inter_pos, inter_neg), "Tofacitinib") + # nolint
    plot_layout(nrow = 1)
ggsave("Fig2/h_vol.pdf", final, width = 10, height = 2.5)

## p2
select_gene <- c(inter_neg, inter_pos)
mdb_c2 <- msigdbr(species = "Homo sapiens", category = "C2")
mdb_kegg <- mdb_c2[grep("^KEGG", mdb_c2$gs_name), ]
fgsea_sets <- mdb_kegg %>% split(x = .$gene_symbol, f = .$gs_name)

marker_a <- marker_a[which(marker_a$gene %in% select_gene), ]
marker_a$avg_log2FC <- ifelse(marker_a$cluster == "OA_BT",
    -1 * marker_a$avg_log2FC, marker_a$avg_log2FC
)

gene0 <- marker_a %>%
    arrange(desc(avg_log2FC)) %>%
    dplyr::select(gene, avg_log2FC)

ranks <- deframe(gene0)
fgseaRes <- fgsea(fgsea_sets, stats = ranks, nperm = 1000)

func <- c(
    "KEGG_CHEMOKINE_SIGNALING_PATHWAY",
    "KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION"
)

final <- plotenrich(func[1]) + plotenrich(func[2]) +
    plot_layout(nrow = 2)
final
ggsave("Fig2/j_gsea.pdf", final, width = 4, height = 5)

## p3
go_result <- read.table("input/macro_deg_go.txt", header = T, sep = "\t")
go_result$iterm <- factor(go_result$iterm, levels = rev(go_result$iterm))

p <- ggplot(go_result, aes(x = -log(qvalue, 10), y = iterm)) +
    geom_bar(stat = "identity", fill = "#AC5987") +
    theme_classic() +
    scale_x_continuous(expand = c(0, 0))
p
ggsave("Fig2/macrophage_go_result.pdf", p, width = 5, height = 3.5)

marker_a <- read.table("input/tcell_deg_OA_BT-RA_BT.txt", header = T, sep = "\t") # nolint
marker_b <- read.table("input/tcell_deg_RA_BT-RA_AT.txt", header = T, sep = "\t") # nolint
marker_c <- read.table("input/tcell_deg_RA_BT-RA_AT_ada.txt", header = T, sep = "\t") # nolint
marker_d <- read.table("input/tcell_deg_RA_BT-RA_AT_tof.txt", header = T, sep = "\t") # nolint

inter_neg <- Reduce(intersect, list(gene_subset(marker_a, "OA_BT"), gene_subset(marker_b, "RA_AT"), gene_subset(marker_c, "RA_AT"), gene_subset(marker_d, "RA_AT"))) # nolint
inter_pos <- Reduce(intersect, list(gene_subset(marker_a, "RA_BT"), gene_subset(marker_b, "RA_BT"), gene_subset(marker_c, "RA_BT"), gene_subset(marker_d, "RA_BT"))) # nolint

label <- c("S100A9", "S100A8", "SOCS3", "SLC2A3", "GBP1", "STAT1") # nolint
label1 <- c(label, "IRF1", "JUN", "CD44", "JUND", "JAK3")
label2 <- c(label, "IRF1", "ISG15", "GBP5", "GBP4", "JAK3", "GBP1", "GZMA", "IFI16") # nolint
label3 <- c(label, "CXCL8", "SPP1", "VEGFA", "CCL2", "CXCL2", "CD44", "CXCL3", "CCL7", "MIF") # nolint
label4 <- c(label, "SPP1", "CCL2", "MT2A", "CXCL3", "CXCL2", "GBP1", "SLAMF9", "IFI6", "CXCL8", "CCL7") # nolint

final <- vol_graph(marker_a, "OA_BT", label1, c(inter_pos, inter_neg), "") +
    vol_graph(marker_b, "RA_AT", label2, c(inter_pos, inter_neg), "") +
    vol_graph(marker_c, "RA_AT", label3, c(inter_pos, inter_neg), "Adalimumab") + # nolint
    vol_graph(marker_d, "RA_AT", label4, c(inter_pos, inter_neg), "Tofacitinib") + # nolint
    plot_layout(nrow = 1)
ggsave("Fig2/h_vol1.pdf", final, width = 10, height = 2.5)
#############################################################


###################### DEG markers ##########################
vp_case1 <- function(
    cell, gene_signature, file_name, ncol, cc) {
    sub_data <- subset(data_input, cells = cell)
    sub_data <- subset(sub_data, cells = rownames(sub_data@meta.data)[which(sub_data$celltype == cc)]) # nolint
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

vp_case1(
    cell = rownames(data_input@meta.data)[which(data_input$ACR20 %in% c("Y", "N"))], # nolint
    gene_signature = c(
        "STAT1", "FOS", "IRF1", "GBP5"
    ), file_name = "Fig2/ACR20_Y_macrophage", ncol = 4, cc = "Macrophage"
)

vp_case1(
    cell = rownames(data_input@meta.data)[which(data_input$ACR20 %in% c("Y", "N"))], # nolint
    gene_signature = c(
        "GZMA", "GZMH", "GZMK", "CCL5"
    ), file_name = "Fig2/ACR20_Y_tcell", ncol = 4, cc = "T cell"
)
#############################################################
