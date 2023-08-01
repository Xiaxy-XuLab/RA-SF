#!/usr/bin/env R
# coding utf-8

###################### Color Setting ###########################
group_color3 <- c(
    "#7FBCE8", "#E26260", "#DDAF54"
)
celltype_color <- c(
    "#DCA263CC", "#B87475CC", "#8491B4E5", "#3C5488E5",
    "#91D1C2E5", "#00A087E5", "#F39B7FE5", "#7E6148FF"
)
group_color5 <- c(
    "#509AAF", "#D1473D", "#F0C05F", "#A04246", "#F7EF60"
)
dotplot_color <- c(
    "#FFFFFF", "#F1C4AF", "#DE785C", "#C83A2F", "#641C19"
)
deg_color <- c(
    "#8F2A47", "#B83A4D", "#D25C52", "#E27E56",
    "#ECA86B", "#F4CB85", "#F8E8A2", "#FAF8C7", "#EBF0AF",
    "#CEE2A2", "#ABD3A6", "#82C3A5", "#609EB0", "#4C78B1",
    "#5C519B"
)
macrophage_color <- c(
    "#F39B7FE5", "#C6B4D3", "#F1B26E",
    "#4F89BB", "#609EB0", "#C66654"
)
tcell_color <- c(
    "#4DBBD5E5", "#F39B7FE5", "#E64B35E5", "#951C55",
    "#7E6148E5", "#8491B4E5", "#91D1C2E5", "#00A087E5",
    "#3C5488E5", "#EC762F"
)
################################################################


################### clinical indicators ########################
das28_index <- c(
    NA, NA, NA,
    5.29, 3.71, 4.43, 3.45, 4.15, 3.28,
    5.16, 3.7, 4.22, 3.63, 6.93, 4.42
)
crp_index <- c(
    6.94, 6.10, 6.06,
    17.5, 3.89, 31.6, 15.2, 19.5, 5.12,
    76.1, 13, 5.89, 2.55, 74.3, 10.2
)
sdai_index <- c(
    NA, NA, NA,
    28.75, 14.39, 19.16, 10.52, 16.95, 10.51,
    27.61, 12.3, 19.58, 15.25, 51.43, 19.02
)
cdai_index <- c(
    NA, NA, NA,
    27, 14, 16, 9, 15, 10,
    20, 11, 19, 15, 44, 18
)
samples <- c(
    "OA_1_1", "OA_2_1", "OA_3_1", "RA_1_1",
    "RA_1_2", "RA_2_1", "RA_2_2", "RA_3_1",
    "RA_3_2", "RA_4_1", "RA_4_2", "RA_5_1",
    "RA_5_2", "RA_6_1", "RA_6_2"
)

clinical_index_single <- data.frame(
    row.names = samples,
    das28_index = das28_index, crp_index = crp_index,
    sdai_index = sdai_index, cdai_index = cdai_index,
    subtype = factor(
        c(rep("OA-BT", 3), rep(c("RA-BT", "RA-AT"), 6)),
        levels = c("OA-BT", "RA-BT", "RA-AT")
    ),
    patient = c(paste("OA", 1:3), rep(paste("RA", 1:6), each = 2)),
    lab = c(NA, NA, NA, 1, NA, 2, NA, 3, NA, 4, NA, 5, NA, 6, NA),
    drug = factor(
        c(rep("None", 3), rep(c("Adalimumab", "Tofacitinib"), each = 6)),
        levels = c("None", "Adalimumab", "Tofacitinib")
    )
)
################################################################


######################### signature ############################
M1 <- c(
    "IL12", "IL23", "TNF", "IL6", "CD86", "MHCII", "IL1B",
    "SPP1", "CCL2", "STAT1", "MARCO", "iNOS", "CD64", "CD80",
    "CXCR10", "CXCL9", "CXCL10", "CXCL11", "IL1A", "TNF", "CCL5",
    "IRF5", "IRF1", "CD40", "IDO1", "KYNU", "CCR7", "CD45",
    "CD68", "CD115", "HLA-DRA", "CD205", "CD14"
)
tmarkerlist <- list(
    "Naive T" = c("IL7R", "CCR7"),
    "ZNF683+ CD8+ Trm" = c("ZNF683", "CCL5", "ITGA1", "CD8A", "GZMH", "TRGC2"), # nolint
    "GZMK+ CD8+ Tem1" = c("GZMK", "CD8A", "EOMES", "GZMA", "GZMH", "NKG7", "CCL5", "KLRG1", "GZMB"), # nolint
    "GZMK+ CD8+ Tem2" = c("MYADM", "LMNA", "ANXA1", "TUBA1A", "JUND", "KLF6", "GZMK", "CCL5", "RGCC", "NKG7"), # nolint
    "GZMK+ CD8+ NKTem" = c("KLRD1", "TRGC2", "TRDC", "CCL4", "CCL5", "NKG7", "KLRB1", "KLRG1"), # nolint
    "GNLY+ NK" = c("GNLY", "TYROBP", "KLRC1", "FCER1G", "KLRD1", "TRDC", "XCL2", "XCL1", "KLRF1", "KLRB1", "CEBPD"), # nolint
    Treg = c("FOXP3", "IL2RA", "RTKN2", "TBC1D4", "CTLA4", "TIGIT"), # nolint
    "CD4+ Tex" = c("CXCL13", "CD4", "MAF", "TNFRSF18", "CTLA4", "TIGIT"), # nolint
    "Proliferating CD8+ T" = c("STMN1", "TYMS", "MKI67", "TOP2A", "UBE2C", "CD8A"), # nolint
    "Proliferating NK" = c("UBE2C", "BIRC5", "KIF20A", "TYMS", "SPC25", "TOP2A", "MKI67", "KLRC1", "KLRD1", "NKG7") # nolint
)
diagnostic_signature <- unique(c(
    "SPP1", "CD44", "PTGER4", "LAMP1", "FAM3C", "CD74", "COPA",
    "ICAM1", "SPN", "ITGAL", "TNFRSF1A", "FASLG", "NRP1", "VEGFB",
    "FLT1", "CLEC2D", "TNF", "ICOS", "FAS", "CD40", "CD40LG", "LGALS9",
    "HAVCR2", "SORL1", "PLXNB2", "SEMA4D", "MIF", "TNFRSF14", "PDCD1",
    "TNFRSF1B", "TFRC", "TNFSF13B", "HLA-F", "LILRB2", "LILRB1", "CCL4",
    "CCR5", "SLC7A1", "IL15", "BTLA", "TNFRSF14", "CXCR6", "CXCL16",
    "CCL7", "IFNG", "CD2", "CD58", "PTPRC", "MRC1", "CD28", "CD86", "SELL",
    "SELPLG", "CD99", "PILRA", "SIGLEC1", "CD6", "ALCAM", "IGFL2", "IGFLR1",
    "CCL2", "CCL3", "CCL4", "CCL7", "CXCL16", "CCL3L1", "CCR1", "CCR5", "SDC4",
    "CCR2", "CXCR3", "CXCL13", "CXCR6", "IL1B", "IL1RAP", "ICAM3", "HLA-DRA",
    "TNFSF12", "IL1RN", "MMP9", "ITGB1", "HBEGF", "HLA-A", "A2M", "ITGB2",
    "SEMA4D", "SELL", "PECAM1", "TNFSF13B", "HLA-E", "ADAM17", "ITGAM",
    "ITGB7", "IL18", "IL15", "PIK3CB", "OSM", "RAB27A", "NCOR1", "YY1",
    "BTG2", "B4GALT1", "GNAS", "LRRFIP1", "IKZF3", "PBXIP1", "MAF", "CLEC2B",
    "MLLT6", "JUND", "TPM4", "SOD2", "S100A9", "PTPN1", "ITGA4", "GZMB", "VIM",
    "STAT5B", "STAT3", "SOCS3", "SBNO2", "S100A8", "NEAT1", "MCL1", "KLF6",
    "JUNB", "JUN", "IRF1", "IL6ST", "FAS", "CDKN1B", "BRD2", "BHLHE40", "AHNAK",
    "LTB", "TIMP1", "LGALS3", "CD320", "AREG", "SEMA4D", "FN1", "ADAM12",
    "ADAM17", "HMGB1", "NAMPT", "UTS2", "ANXA1", "EMP1", "ACTN1", "FTH1",
    "LDHA", "RTN4", "PIM3", "ADAM9", "NEDD9", "JMJD1C", "ZMIZ1", "SKIL",
    "FOSL2", "TAGLN2", "CD55", "ADAM17", "VCAN", "TAPBP", "CXCL3", "TNS3",
    "NDRG1", "S100A10", "RBPJ", "IL6R", "NEAT1", "KDM6B", "BHLHE40", "TPM4",
    "TIMP1", "SERPINE1", "ITGAX", "ITGA5", "GSTO1", "VEGFA", "TLR2", "STAT1",
    "SOCS3", "PTGER2", "PLAUR", "PIM1", "NFKBIA", "IRF1", "ICAM1", "GBP1",
    "FCGR3A", "CXCL2", "FASLG", "TNFRSF10B", "TNFRSF10A", "TNFSF10", "CSF1R",
    "CSF1", "CD55", "ADGRE5", "ALOX5", "ALOX5AP", "SELPLG", "C5AR1", "RPS19",
    "LTBR", "TIGIT", "NECTIN2", "CD226", "PDCD1", "SLC7A1", "TNFRSF10D",
    "CSF1", "SIRPA", "TNFSF12", "TNFRSF12A", "CTLA4"
))
disease_severity_signature <- c("SPP1", "STAT1", "ETV7", "IRF1", "EOMES")
ACR20_Y_signature <- c(
    "STAT1", "IRF1", "GBP5", "GZMA", "GZMH",
    "GZMK", "CCL5", "FCGR3A", "HLA-A", "MT2A", "GBP1",
    "IFITM3", "IFI16", "IFNGR1", "LAG3", "PDCD1", "NFKBIA"
)
ACR20_N_signature <- c(
    "HIF1A", "IFI30", "IL1B", "JUN", "CXCL13", "TIGIT",
    "CTLA4", "S100A9", "MIF", "IFITM1", "CXCR4", "JUND", "ISG15"
)
#################################################################
