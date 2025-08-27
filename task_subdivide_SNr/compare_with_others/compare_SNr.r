# project mouse whole cortex
# subdivision of SNr
# compare with previous studies

SNr_gao <- nat::read.nrrd(
  file = "/run/media/admin/LeGAO2/project_mouse_cortex/analysis/task_subdivide_SNr/K6_coronal.nrrd",
  ReadByteAsRaw = FALSE)

dim(SNr_gao)


SNr_foster_label <-
  snp::snp_read_itksnap_label("/run/media/admin/LeGAO2/project_mouse_cortex/analysis/task_subdivide_SNr/compare_with_others/Foster_SNr_labels.label")


SNr_foster_label$ColorHex <-
  rgb(red = SNr_foster_label$R,
      green = SNr_foster_label$G,
      blue = SNr_foster_label$B,
      maxColorValue = 255)



##
SNr_level_81 <-
  nat::read.nrrd(file = "/run/media/admin/LeGAO2/project_mouse_cortex/analysis/task_subdivide_SNr/compare_with_others/Allen_level_81_50um_seg.nrrd",
                 ReadByteAsRaw = F)

SNr_level_83 <-
  nat::read.nrrd(file = "/run/media/admin/LeGAO2/project_mouse_cortex/analysis/task_subdivide_SNr/compare_with_others/Allen_level_83_50um_seg.nrrd",
                 ReadByteAsRaw = F)

SNr_level_85 <-
  nat::read.nrrd(file = "/run/media/admin/LeGAO2/project_mouse_cortex/analysis/task_subdivide_SNr/compare_with_others/Allen_level_85_50um_seg.nrrd",
                 ReadByteAsRaw = F)

SNr_level_87 <-
  nat::read.nrrd(file = "/run/media/admin/LeGAO2/project_mouse_cortex/analysis/task_subdivide_SNr/compare_with_others/Allen_level_87_50um_seg.nrrd",
                 ReadByteAsRaw = F)

SNr_level_89 <-
  nat::read.nrrd(file = "/run/media/admin/LeGAO2/project_mouse_cortex/analysis/task_subdivide_SNr/compare_with_others/Allen_level_89_50um_seg.nrrd",
                 ReadByteAsRaw = F)

SNr_level_91 <-
  nat::read.nrrd(file = "/run/media/admin/LeGAO2/project_mouse_cortex/analysis/task_subdivide_SNr/compare_with_others/Allen_level_91_50um_seg.nrrd",
                 ReadByteAsRaw = F)


dim(SNr_level_81)


SNr_foster <- array(dim = c(228, 160, 6))
SNr_foster[,,1] <- SNr_level_81
SNr_foster[,,2] <- SNr_level_83
SNr_foster[,,3] <- SNr_level_85
SNr_foster[,,4] <- SNr_level_87
SNr_foster[,,5] <- SNr_level_89
SNr_foster[,,6] <- SNr_level_91

nat::write.nrrd(x = SNr_foster, file = "/tmp/SNr_foster.nrrd", dtype = "byte")





#part_SNr_gao <- SNr_gao[,,c(160, 164, 168, 172, 176, 180, 184)]

part_SNr_gao <- SNr_gao[,,c(162, 166, 170, 174, 178, 182)]


matrix_foster_gao <- matrix(nrow = 6, ncol = 14, data = 0)
rownames(matrix_foster_gao) <- paste0("Gao-", 1:6)
colnames(matrix_foster_gao) <- SNr_foster_label$LABEL

for (iLabel in 1:6) {
  index_voxel <- which(part_SNr_gao == iLabel)
  for (jLabel in 1:14) {
    matrix_foster_gao[iLabel, jLabel] <-
      length(which(SNr_foster[index_voxel] == jLabel))
  }
}




library(ComplexHeatmap)

Heatmap(matrix = matrix_foster_gao,
        cluster_rows = F,
        cluster_columns = F,
        row_names_side = "left",
        column_names_side = "top")





# enrichment analysis ----

SNr_matrix_enrich_p <- matrix(nrow = 6, ncol = 14)
rownames(SNr_matrix_enrich_p) <- paste0("Gao-", 1:6)
colnames(SNr_matrix_enrich_p) <- SNr_foster_label$LABEL

SNr_matrix_enrich_e <- matrix(nrow = 6, ncol = 14)
rownames(SNr_matrix_enrich_e) <- paste0("Gao-", 1:6)
colnames(SNr_matrix_enrich_e) <- SNr_foster_label$LABEL


iLabel <- 1
for (iLabel in 1:6) {
  jLabel <- 1
  for (jLabel in 1:14) {
    tmp1 <- matrix_foster_gao[iLabel, jLabel]
    tmp2 <- sum(matrix_foster_gao[iLabel,]) - tmp1
    tmp3 <- sum(matrix_foster_gao[,jLabel]) - tmp1
    tmp4 <- sum(matrix_foster_gao) - tmp1 - tmp2 - tmp3

    tmpTest <-
      fisher.test(x = matrix(data = c(tmp1, tmp2, tmp3, tmp4),
                             byrow = TRUE, nrow = 2))

    SNr_matrix_enrich_p[iLabel, jLabel] <- tmpTest$p.value
    SNr_matrix_enrich_e[iLabel, jLabel] <- tmpTest$estimate
  }
}



#matrix_enrich_e[is.infinite(matrix_enrich_e)] <- 0

Heatmap(matrix = -log10(SNr_matrix_enrich_p),
        name = "-log(P)",
        cell_fun = function(j, i, x, y, width, height, fill) {
          if (SNr_matrix_enrich_p[i,j] < 0.001 &
              (SNr_matrix_enrich_e[i,j] > 1 | is.infinite(SNr_matrix_enrich_e[i,j]))) {
            #grid.text(sprintf("%.1f", matrix_enrich_p[i, j]), x, y, gp = gpar(fontsize = 10))
            grid.circle(x = x, y = y,
                        r = 0.04,
                        gp = gpar(fill = NA, color = "black"))
          }
        },
        cluster_rows = F,
        cluster_columns = F,
        row_names_side = "left",
        #row_labels = gsub("\\.", "\\\\.", rownames(matrix_enrich_e)),
        column_names_side = "top")


