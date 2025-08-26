
# test. read LUT file generated in FIJI ----
lut_raw <-
  readBin(con = "/run/media/admin/LeGAO2/project_mouse_cortex/analysis/task_subdivide_GPe/compare_with_others/Untitled.lut",
        what = "raw", n = 768)

lut_vals <- as.integer(lut_raw)

red <- lut_vals[1:256]
green <- lut_vals[257:512]
blue <- lut_vals[513:768]

lut_df <- data.frame(R = red, G = green, B = blue)
head(lut_df)


# 1. write LUT for GPe subdomains from Foster et al., 2021, Nature ----
itksnap_label <-
  snp::snp_read_itksnap_label("/run/media/admin/LeGAO2/project_mouse_cortex/analysis/task_subdivide_GPe/compare_with_others/Untitled.label")

R_0_255 <- c(0, itksnap_label$R, rep(0, (256 - 1 - nrow(itksnap_label))))
G_0_255 <- c(0, itksnap_label$G, rep(0, (256 - 1 - nrow(itksnap_label))))
B_0_255 <- c(0, itksnap_label$B, rep(0, (256 - 1 - nrow(itksnap_label))))

R_bytes <- as.raw(as.integer(R_0_255))
G_bytes <- as.raw(as.integer(G_0_255))
B_bytes <- as.raw(as.integer(B_0_255))

lut_raw <- c(R_bytes, G_bytes, B_bytes)

writeBin(object = lut_raw,
         con = "/run/media/admin/LeGAO2/project_mouse_cortex/analysis/task_subdivide_GPe/compare_with_others/Untitled.lut")


# 2. read subdivisions in current strudy and Foster et al. ----
gao_GPe <-
  nat::read.nrrd(file = "/run/media/admin/LeGAO2/project_mouse_cortex/analysis/task_subdivide_GPe/coronal_substack_116_120_124_128_132_136.nrrd",
                 ReadByteAsRaw = F)

table(gao_GPe)

foster_GPe <-
  nat::read.nrrd(file = "/run/media/admin/LeGAO2/project_mouse_cortex/analysis/task_subdivide_GPe/compare_with_others/Substack_116_120_124_128_132_136_seg.nrrd",
                 ReadByteAsRaw = F)

table(foster_GPe)

# 3. calculate number of voxels assigned as different subdivisions ----
GPe_matrix_correspondence <- matrix(nrow = 6, ncol = 36)
rownames(GPe_matrix_correspondence) <- paste0("Gao-", 1:6)
colnames(GPe_matrix_correspondence) <- paste0("Foster-", 1:36)


GPe_matrix_Dice <- matrix(nrow = 6, ncol = 36)
rownames(GPe_matrix_Dice) <- paste0("Gao-", 1:6)
colnames(GPe_matrix_Dice) <- paste0("Foster-", 1:36)


for (iSubdomain in 1:6) {
  index_voxel <- which(gao_GPe == iSubdomain)

  for (jSubdomain in 1:36) {
    GPe_matrix_correspondence[iSubdomain, jSubdomain] <-
      length(which(foster_GPe[index_voxel] == jSubdomain))
    GPe_matrix_Dice[iSubdomain, jSubdomain] <-
      2 * length(which(foster_GPe[index_voxel] == jSubdomain)) / (length(index_voxel) + length(which(foster_GPe == jSubdomain)))
  }
}

library(ComplexHeatmap)

Heatmap(matrix = GPe_matrix_Dice,
        name = "Number of voxels",
        cluster_rows = F,
        cluster_columns = F,
        row_names_side = "left",
        column_names_side = "top")

# 4. enrichment analysis ----

GPe_matrix_enrich_p <- matrix(nrow = 6, ncol = 36)
rownames(GPe_matrix_enrich_p) <- paste0("Gao-", 1:6)
colnames(GPe_matrix_enrich_p) <- paste0("Foster-", 1:36)

GPe_matrix_enrich_e <- matrix(nrow = 6, ncol = 36)
rownames(GPe_matrix_enrich_e) <- paste0("Gao-", 1:6)
colnames(GPe_matrix_enrich_e) <- paste0("Foster-", 1:36)


iLabel <- 1
for (iLabel in 1:6) {
  jLabel <- 1
  for (jLabel in 1:36) {
    tmp1 <- GPe_matrix_correspondence[iLabel, jLabel]
    tmp2 <- sum(GPe_matrix_correspondence[iLabel,]) - tmp1
    tmp3 <- sum(GPe_matrix_correspondence[,jLabel]) - tmp1
    tmp4 <- sum(GPe_matrix_correspondence) - tmp1 - tmp2 - tmp3

    tmpTest <-
      fisher.test(x = matrix(data = c(tmp1, tmp2, tmp3, tmp4),
                             byrow = TRUE, nrow = 2))

    GPe_matrix_enrich_p[iLabel, jLabel] <- tmpTest$p.value
    GPe_matrix_enrich_e[iLabel, jLabel] <- tmpTest$estimate
  }
}


GPe_matrix_enrich_e[is.infinite(GPe_matrix_enrich_e)] <- 0


Heatmap(matrix = -log10(GPe_matrix_enrich_p),
  name = "-log10(P)",
        #matrix = GPe_matrix_correspondence,
        #name = "Number of voxels",
        cell_fun = function(j, i, x, y, width, height, fill) {
          if (GPe_matrix_enrich_p[i,j] < 0.001 &
              (GPe_matrix_enrich_e[i,j] > 1 | is.infinite(GPe_matrix_enrich_e[i,j]))) {
            #grid.text(sprintf("%.1f", matrix_enrich_p[i, j]), x, y, gp = gpar(fontsize = 10))
            grid.circle(x = x, y = y,
                        r = 0.04,
                        gp = gpar(fill = NA, color = "black"))
          }
        },
        cluster_rows = F,
        cluster_columns = F,
        row_names_side = "left",
        column_names_side = "top")


