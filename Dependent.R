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
