

chon_ontology <-
  read.csv("/home/admin/shared/25750983/UnifiedAtlas_Label_ontology_v2.csv")

index_CP <- which(chon_ontology$name == "Caudate Putamen ")

id_CP <- chon_ontology$id[index_CP]

chon_ontology[which(chon_ontology$structure_id_path == id_CP),]

flag_CP <- stringr::str_detect(string = chon_ontology$name, pattern = "Caudate*")

table(flag_CP)

which(flag_CP == TRUE)

chon_ontology[559,]


chon_anno <-
  nat::read.nrrd(file = "~/shared/25750983/UnifiedAtlas_Label_v2_20um-isotropic_flipZ.nrrd", ReadByteAsRaw = F)

chon_anno_100um <-
  nat::read.nrrd(file = "~/shared/25750983/UnifiedAtlas_Label_v2_100um_coronal.nrrd",
                 ReadByteAsRaw = FALSE)


gao_anno_100um <-
  nat::read.nrrd(file = "/run/media/admin/LeGAO2/project_mouse_cortex/analysis/task_subdivide_striatum/striatum_subdomains_K14.nrrd",
                 ReadByteAsRaw = FALSE)

chon_ontology_CP <-
  read.csv("/run/media/admin/LeGAO2/project_mouse_cortex/analysis/task_subdivide_striatum/chon_CP_ontology.csv")



region_size <- vector(length = nrow(chon_ontology_CP))
names(region_size) <- chon_ontology_CP$name
for (id_region in 1:nrow(chon_ontology_CP)) {
  tmp_index_voxel <- which(chon_anno == chon_ontology_CP$id[id_region])
  region_size[id_region] <- length(tmp_index_voxel)
}



region_size2 <- vector(length = nrow(chon_ontology_CP))
names(region_size2) <- chon_ontology_CP$name
for (id_region in 1:nrow(chon_ontology_CP)) {
  tmp_index_voxel <- which(chon_anno_100um == chon_ontology_CP$id[id_region])
  region_size2[id_region] <- length(tmp_index_voxel)
}

names(region_size2[region_size2 > 0])

chon_ontology_CP_valid <- chon_ontology_CP[region_size2 > 0, ]

matrix_gao_chon <- matrix(nrow = 14, ncol = 29, data = 0)
rownames(matrix_gao_chon) <- paste0("Subdomain", 1:14)
colnames(matrix_gao_chon) <- chon_ontology_CP_valid$acronym

for (iLabel in 1:14) {
  index_voxel <- which(gao_anno_100um == iLabel)
  for (jLabel in 1:29) {
    matrix_gao_chon[iLabel, jLabel] <-
      length(which(chon_anno_100um[index_voxel] == chon_ontology_CP_valid$id[jLabel]))
  }
}

matrix_gao_chon_normalized <- apply(matrix_gao_chon, 2, function(x) x/sum(x))

library(ComplexHeatmap)

Heatmap(matrix = t(matrix_gao_chon_normalized),
        cluster_rows = F,
        cluster_columns = F,
        row_names_side = "left",
        column_names_side = "top")


# enrichment analysis ----

matrix_enrich_p <- matrix(nrow = 14, ncol = 29)
rownames(matrix_enrich_p) <- paste0("Subdomain", 1:14)
colnames(matrix_enrich_p) <- chon_ontology_CP_valid$acronym

matrix_enrich_e <- matrix(nrow = 14, ncol = 29)
rownames(matrix_enrich_e) <- paste0("Subdomain", 1:14)
colnames(matrix_enrich_e) <- chon_ontology_CP_valid$acronym


iLabel <- 1
for (iLabel in 1:14) {
  jLabel <- 1
  for (jLabel in 1:29) {
    tmp1 <- matrix_gao_chon[iLabel, jLabel]
    tmp2 <- sum(matrix_gao_chon[iLabel,]) - tmp1
    tmp3 <- sum(matrix_gao_chon[,jLabel]) - tmp1
    tmp4 <- sum(matrix_gao_chon) - tmp1 - tmp2 - tmp3

    tmpTest <-
      fisher.test(x = matrix(data = c(tmp1, tmp2, tmp3, tmp4), byrow = TRUE, nrow = 2))

    matrix_enrich_p[iLabel, jLabel] <- tmpTest$p.value
    matrix_enrich_e[iLabel, jLabel] <- tmpTest$estimate
  }
}



Heatmap(matrix = log2(matrix_enrich_e + 1),
        name = "Log2(Odds ratio + 1)",
        cell_fun = function(j, i, x, y, width, height, fill) {
          if (matrix_enrich_p[i,j] < 0.0001 & matrix_enrich_e[i,j] > 3) {
            #grid.text(sprintf("%.1f", matrix_enrich_p[i, j]), x, y, gp = gpar(fontsize = 10))
            grid.circle(x = x, y = y, r = 0.01, gp = gpar(fill = NA, color = "black"))
          }
        },
        cluster_rows = F,
        cluster_columns = F,
        row_names_side = "left",
        column_names_side = "top")




