# project mouse whole cortex
# subdivide GPi
#

## 1. encode voxels ----
library(nat)

mask_GPi <-
  nat::read.nrrd(
    file = "/run/media/admin/LeGAO2/project_mouse_cortex/analysis/task_subdivide_GPi/structure_1031_GPi_50um_ccf2017.nrrd",
    ReadByteAsRaw = FALSE)

dim(mask_GPi)

mask_GPi[,,115:228] <- 0

index_voxel <- which(mask_GPi == 1)

label_voxel <- 1:length(index_voxel)

anno_GPi <- mask_GPi

anno_GPi[index_voxel] <- label_voxel

nat::write.nrrd(
  x = anno_GPi,
  file = "/run/media/admin/LeGAO2/project_mouse_cortex/analysis/task_subdivide_GPi/data_reproduce/anno_GPi.nrrd",
  enc = "gzip",
  dtype = "ushort")



## 2. obtain neurons ----
index_PT_GPi <-
  intersect(which(shared_proj_Combined[,snp::getStructureFromAcronym("l_GPi")] > 100),
            df_All_index_AxonType$Index[df_All_index_AxonType$AxonType == "PT"])


saveRDS(
  object = index_PT_GPi,
  file = "/run/media/admin/LeGAO2/project_mouse_cortex/analysis/task_subdivide_GPi/data_reproduce/index_PT_GPi.rds",
  version = 3,
  compress = "gzip")

## 3. calculate axon lengths in the cubes ----
proj_cube_GPi <- snp::snp_calculateProjections(
  PathOfNeuron = shared_SWCPath_axon_allen[index_PT_GPi],
  Hemisphere = info_tmp$Hemisphere[index_PT_GPi],
  StructureIDOfTargets = 1:length(index_voxel),
  Annotation = anno_GPi,
  Resolution = 50,
  NumOfCores = 3)


saveRDS(
  object = proj_cube_GPi,
  file = "/run/media/admin/LeGAO2/project_mouse_cortex/analysis/task_subdivide_GPi/data_reproduce/proj_cube_GPi.rds",
  version = 3,
  compress = "gzip")


## 4. exclude cubes with no axon inputs ----

index_cube_valid <- which(colSums(proj_cube_GPi) > 0)

proj_cube_GPi_valid <- proj_cube_GPi[,index_cube_valid]
colnames(proj_cube_GPi_valid) <- paste0("CubeID-", index_cube_valid)


# 4. clustering ----
ff_similarity <-
  snp::snp_bigcorPar(
    x = proj_cube_GPi_valid,
    size = 400,
    verbose = TRUE,
    ncore = 4,
    method = "spearman"
  )

matrix_similarity <- as.matrix(as.data.frame(ff_similarity[, ]))



saveRDS(
  object = matrix_similarity,
  file = "/run/media/admin/LeGAO2/project_mouse_cortex/analysis/task_subdivide_GPi/data_reproduce/matrix_similarity.rds",
  version = 3,
  compress = "gzip")



library(fastcluster)

h <-
  fastcluster::hclust(d = as.dist(1 - matrix_similarity), method = "ward.D2")

plot(h, hang = -1, labels = FALSE)
rect.hclust(tree = h, k = 5)


saveRDS(
  object = h,
  file = "/run/media/admin/LeGAO2/project_mouse_cortex/analysis/task_subdivide_GPi/data_reproduce/h.rds",
  version = 3,
  compress = "gzip")



NumCluster <- 5
label_cluster_cube <- cutree(tree = h, k = NumCluster)

labelOnDend <- unique(label_cluster_cube[h$order])

label_cluster_real <- label_cluster_cube
for (iLabel in 1:NumCluster) {
  label_cluster_real[which(label_cluster_cube == labelOnDend[iLabel])] <-
    iLabel
}


saveRDS(
  object = label_cluster_real,
  file = "/run/media/admin/LeGAO2/project_mouse_cortex/analysis/task_subdivide_GPi/data_reproduce/label_cluster_real.rds",
  version = 3,
  compress = "gzip")


##
cluster_GPi <- mask_GPi
cluster_GPi[index_voxel[index_cube_valid]] <- label_cluster_real

nat::write.nrrd(
  x = cluster_GPi,
  file = "/run/media/admin/LeGAO2/project_mouse_cortex/analysis/task_subdivide_GPi/cluster_GPi_K5.nrrd",
  enc = "gzip",
  dtype = "byte")



# plot heatmap of sperman's r ----
matrix_similarity_ordered <- matrix_similarity[h$order, h$order]

col_fun <- circlize::colorRamp2(breaks = c(-0.05, 0, 0.1),
                                colors = c("blue", "white", "red"))

library(ComplexHeatmap)
ht <-
  Heatmap(
    matrix = matrix_similarity_ordered[seq(1, nrow(matrix_similarity), 5),
                                       seq(1, nrow(matrix_similarity), 5)],
    name = "Spearman's r",
    #col = col_fun,
    clustering_distance_columns = "spearman",
    clustering_distance_rows = "spearman",
    clustering_method_columns = "ward.D2",
    clustering_method_rows = "ward.D2",

    cluster_rows = FALSE,
    #as.dendrogram(h),
    show_row_dend = TRUE,
    row_dend_side = "left",
    show_row_names = FALSE,


    cluster_columns = FALSE,
    #as.dendrogram(h),
    show_column_dend = TRUE,
    show_column_names = FALSE,

    column_dend_side = "top",
    #left_annotation = row_ha,


    width = unit(10, "cm"),
    height = unit(10, "cm"),
    column_names_rot = 45,
    column_names_gp = gpar(fontsize = 10),
    use_raster = TRUE,
    border = TRUE
  )
draw(ht)


# calculate single-neuron projection preference ----
NumNeuron <- length(index_PT_GPi)

proj_domain <- array(data = 0, dim = c(NumNeuron, NumCluster))
volume_domain <- array(data = 0, dim = c(NumCluster, 1))
for (iDomain in 1:NumCluster) {
  volume_domain[iDomain, 1] <-
    length(which(label_cluster_real == iDomain))
  proj_domain[, iDomain] =
    rowSums(proj_cube_GPi_valid[, which(label_cluster_real == iDomain)]) / volume_domain[iDomain, 1]
}


label_domain_neuron <- array(data = 0, dim = c(NumNeuron, 1))
for (iNeuron in 1:NumNeuron) {
  label_domain_neuron[iNeuron, 1] <-
    which(max(proj_domain[iNeuron, ]) == proj_domain[iNeuron, ])[1]
}

# flatmap whole cortex ----
ColorSubdomain <-
  c(
    "1" = "#ff0000",
    "2" = "#00ff00",
    "3" = "#0000ff",
    "4" = "#ff7f00",
    "5" = "#00ffff",
    "6" = "#ff00ff",
    "7" = "#007f7f"
  )

filename <- paste0("/tmp/GPi_K5.svg")
#
U <- shared_soma_XYZUV_CC$U[index_PT_GPi]
V <- shared_soma_XYZUV_CC$V[index_PT_GPi]
Color <- "red"

Radius <- 6

conn <- file("/home/admin/shared/flatmap_Isocortex_template.svg")
svg_PFC <- readLines(conn)
close(conn)
fileConn <- file(filename, "w")
NumElement <- length(U)
svg_element <- c()
if (length(Color) == 1) {
  Color <- rep(Color, NumElement)
}
if (length(Radius) == 1) {
  Radius <- rep(Radius, NumElement)
}
for (iElement in 1:NumElement) {
  svg_element <-
    c(
      svg_element,
      sprintf(
        "<circle cx=\"%f\" cy=\"%f\" r=\"%f\" style=\"stroke: none; fill: %s\"/>",
        U[iElement],
        V[iElement],
        Radius[iElement],
        ColorSubdomain[label_domain_neuron[iElement]]
      )
    )
}

out <- c(svg_PFC[1:50], svg_element, svg_PFC[51:length(svg_PFC)])
writeLines(out, fileConn, sep = "\n")
close(fileConn)





# plot heatmap of projection density at the region level ----
colnames(proj_domain) <- paste0("Subdomain-", 1:5)

df_neuron_domain <-
  data.frame(SourceRegion = info_tmp$Region[index_PT_GPi], proj_domain)

library(dplyr)
library(tidyr)

region_num <-
  df_neuron_domain %>% group_by(SourceRegion) %>% summarize(Count = n())

source_region_valid <- region_num$SourceRegion[region_num$Count > 2]

df_region_domain <-
  df_neuron_domain %>%
  filter(SourceRegion %in% source_region_valid) %>%
  group_by(SourceRegion) %>%
  summarize_all(mean)

matrix_region_domain <- df_region_domain[, -1]
rownames(matrix_region_domain) <- df_region_domain$SourceRegion


#write.csv(x = matrix_region_domain, file = "./result_heatmap_IT.csv")


library(ComplexHeatmap)

ht_region <-
  Heatmap(
    matrix_region_domain,
    cluster_rows = T,
    clustering_distance_rows = "spearman",
    clustering_method_rows = "ward.D2",
    cluster_columns = F,
    row_names_side = "left",
    column_names_side = "top"
  )

draw(ht_region)
