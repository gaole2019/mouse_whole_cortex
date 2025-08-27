# project mouse whole cortex
# compare GPe subdivisions with others

# 1. write LUT for GPe subdomains from Foster et al., 2021, Nature ----
itksnap_label <-
  snp::snp_read_itksnap_label("/run/media/admin/LeGAO2/project_mouse_cortex/analysis/task_subdivide_striatum/compare_with_others/hintiryan.label")

R_0_255 <- c(0, itksnap_label$R, rep(0, (256 - 1 - nrow(itksnap_label))))
G_0_255 <- c(0, itksnap_label$G, rep(0, (256 - 1 - nrow(itksnap_label))))
B_0_255 <- c(0, itksnap_label$B, rep(0, (256 - 1 - nrow(itksnap_label))))

R_bytes <- as.raw(as.integer(R_0_255))
G_bytes <- as.raw(as.integer(G_0_255))
B_bytes <- as.raw(as.integer(B_0_255))

lut_raw <- c(R_bytes, G_bytes, B_bytes)

writeBin(object = lut_raw,
         con = "/run/media/admin/LeGAO2/project_mouse_cortex/analysis/task_subdivide_striatum/compare_with_others/hintiryan.lut")



# 2. read subdivisions in current strudy and Hintiryan et al. ----
gao_CP <-
  nat::read.nrrd(file = "/run/media/admin/LeGAO2/project_mouse_cortex/analysis/task_subdivide_striatum/gao_Substack_41_53_61_67_74.nrrd",
                 ReadByteAsRaw = F)

table(gao_CP)

hintiryan_CP <-
  nat::read.nrrd(file = "/run/media/admin/LeGAO2/project_mouse_cortex/analysis/task_subdivide_striatum/compare_with_others/Substack_41_53_61_67_74_seg_and_ACB.nrrd",
                 ReadByteAsRaw = F)

table(hintiryan_CP)
#add ACB as the 31th subdivision


# 3. calculate number of voxels assigned as different subdivisions ----

CP_matrix_correspondence <- matrix(nrow = 14, ncol = 31)
rownames(CP_matrix_correspondence) <- paste0("Gao-", 1:14)
colnames(CP_matrix_correspondence) <- c(itksnap_label$LABEL, "ACB")#paste0("Hintiryan-", 1:31)


CP_matrix_Dice <- matrix(nrow = 14, ncol = 31)
rownames(CP_matrix_Dice) <- paste0("Gao-", 1:14)
colnames(CP_matrix_Dice) <- c(itksnap_label$LABEL, "ACB")#paste0("Hintiryan-", 1:31)


for (iSubdomain in 1:14) {
  index_voxel <- which(gao_CP == iSubdomain)

  for (jSubdomain in 1:31) {
    CP_matrix_correspondence[iSubdomain, jSubdomain] <-
      length(which(hintiryan_CP[index_voxel] == jSubdomain))
    CP_matrix_Dice[iSubdomain, jSubdomain] <-
      2 * length(which(hintiryan_CP[index_voxel] == jSubdomain)) / (length(index_voxel) + length(which(hintiryan_CP == jSubdomain)))
  }
}

library(ComplexHeatmap)

Heatmap(matrix = CP_matrix_Dice,
        name = "Number of voxels",
        cluster_rows = F,
        cluster_columns = F,
        row_names_side = "left",
        column_names_side = "top")





# 4. enrichment analysis ----

CP_matrix_enrich_p <- matrix(nrow = 14, ncol = 31)
rownames(CP_matrix_enrich_p) <- paste0("Gao-", 1:14)
colnames(CP_matrix_enrich_p) <- c(itksnap_label$LABEL, "ACB")

CP_matrix_enrich_e <- matrix(nrow = 14, ncol = 31)
rownames(CP_matrix_enrich_e) <- paste0("Gao-", 1:14)
colnames(CP_matrix_enrich_e) <- c(itksnap_label$LABEL, "ACB")


iLabel <- 1
for (iLabel in 1:14) {
  jLabel <- 1
  for (jLabel in 1:31) {
    tmp1 <- CP_matrix_correspondence[iLabel, jLabel]
    tmp2 <- sum(CP_matrix_correspondence[iLabel,]) - tmp1
    tmp3 <- sum(CP_matrix_correspondence[,jLabel]) - tmp1
    tmp4 <- sum(CP_matrix_correspondence) - tmp1 - tmp2 - tmp3

    tmpTest <-
      fisher.test(x = matrix(data = c(tmp1, tmp2, tmp3, tmp4),
                             byrow = TRUE, nrow = 2))

    CP_matrix_enrich_p[iLabel, jLabel] <- tmpTest$p.value
    CP_matrix_enrich_e[iLabel, jLabel] <- tmpTest$estimate
  }
}





Heatmap(matrix = -log10(CP_matrix_enrich_p),
        name = "-log10(P)",
        #matrix = GPe_matrix_correspondence,
        #name = "Number of voxels",
        cell_fun = function(j, i, x, y, width, height, fill) {
          if (CP_matrix_enrich_p[i,j] < 0.001 &
              (CP_matrix_enrich_e[i,j] > 1 | is.infinite(CP_matrix_enrich_e[i,j]))) {
            #grid.text(sprintf("%.1f", matrix_enrich_p[i, j]), x, y, gp = gpar(fontsize = 10))
            grid.circle(x = x, y = y,
                        r = 0.02,
                        gp = gpar(fill = NA, color = "black"))
          }
        },
        cluster_rows = F,
        cluster_columns = F,
        row_names_side = "left",
        column_names_side = "top")

