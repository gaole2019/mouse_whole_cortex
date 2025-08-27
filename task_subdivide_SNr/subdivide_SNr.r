# project mouse cortex
# figure 4
# PT
# subdivide SNr


# 1. ----
Hemisphere_Path_PT_SNr <-
  data.frame(Hemisphere = info_tmp$Hemisphere[index_axon_PT_SNr],
             Path = shared_SWCPath_axon_allen[index_axon_PT_SNr])

write.table(
  x = Hemisphere_Path_PT_SNr,
  file = "/run/media/admin/LeGAO/Project_Cortex/task_PT/subdivide_SNr/Hemisphere_Path_PT_SNr.csv",
  row.names = FALSE,
  col.names = FALSE)

Hemisphere_Path_PT_SNr <-
  read.table(
    file = "/run/media/admin/LeGAO/Project_Cortex/task_PT/subdivide_SNr/Hemisphere_Path_PT_SNr.csv")


matrix_proj_SNr <-
  read.csv(file = "/run/media/admin/LeGAO/Project_Cortex/task_PT/subdivide_SNr/axonLength.csv", row.names = 1)
matrix_proj_SNr <- as.matrix(matrix_proj_SNr)

# some cubes are empty
index_cube_valid <- which(colSums(matrix_proj_SNr) > 0)


matrix_proj_SNr_valid <- matrix_proj_SNr[,index_cube_valid]

## Calculate similarity
ff_similarity <-
  snp::snp_bigcorPar(x = matrix_proj_SNr_valid,
                     size = 400, verbose = TRUE,
                     ncore = 2, method = "spearman")

matrix_similarity <- as.matrix(as.data.frame(ff_similarity[,]))


#matrix_similarity = read.table("./matrix_similarity_800.txt")


library(fastcluster)

h <- fastcluster::hclust(d = as.dist(1 - matrix_similarity), method = "ward.D2")
plot(h, hang = -1, labels = FALSE)
rect.hclust(tree = h, k = 6)

#
NumCluster <- 6
label_cluster_cube <- cutree(tree = h, k = NumCluster)

labelOnDend <- unique(label_cluster_cube[h$order])

label_cluster_real <- label_cluster_cube
for (iLabel in 1:NumCluster) {
  label_cluster_real[which(label_cluster_cube == labelOnDend[iLabel])] <- iLabel
}

# 2. tsne ----



ColorStellate <-
  c("L4 spiny stellate"="red","Other"="green")

set.seed(1)
tsne_out <- Rtsne::Rtsne(X = 1 - matrix_similarity, dims = 2,
                         is_distance = T)

plot(
  tsne_out$Y,
  #col = as.factor(dendrite_info$W),
  #col = ColorSourceRegion[match(dendrite_info$BrainRegion, names(ColorSourceRegion))],
  #col = ColorPyN[match(dendrite_info_revised$IsApicalDendriteVisible,names(ColorPyN))],
  col = ColorSubdomain[match(label_cluster_real, names(ColorSubdomain))],
  #col = as.factor(label_cluster_real),
  #pch = 17,
  cex = 0.6,
  asp = 1)



# 3. heatmap ----

matrix_similarity_ordered = matrix_similarity[h$order, h$order]

library(circlize)

col_fun <-
  circlize::colorRamp2(breaks = c(-0.1,0.1,0.3), colors = c("blue","white","red"))

library(ComplexHeatmap)
ht <-
  Heatmap(
    matrix = matrix_similarity_ordered[seq(1,nrow(matrix_similarity),5),
                                       seq(1,nrow(matrix_similarity),5)],
    name = "Spearman's r",
             col = col_fun,
             clustering_distance_columns = "spearman",
             clustering_distance_rows = "spearman",
             clustering_method_columns = "ward.D2",
             clustering_method_rows = "ward.D2",

             cluster_rows = FALSE, #as.dendrogram(h),
             show_row_dend = TRUE,
             row_dend_side = "left",
             show_row_names = FALSE,


             cluster_columns = FALSE, #as.dendrogram(h),
             show_column_dend = TRUE,
             show_column_names = FALSE,

             column_dend_side = "top",
             #left_annotation = row_ha,


             width = unit(10, "cm"),
             height = unit(10, "cm"),
             column_names_rot = 45,
             column_names_gp = gpar(fontsize = 10),
             use_raster = TRUE,
             border = TRUE)
draw(ht)



# 4. projections to SNr subdomains ----
label_domain_full <- 1:ncol(matrix_proj_SNr)
label_domain_full[-index_cube_valid] <- 0
for (iLabel in 1:NumCluster) {
  label_domain_full[index_cube_valid[which(label_cluster_real == iLabel)]] <- iLabel
}


write.table(
  x = label_domain_full,
  file = "/run/media/admin/LeGAO/Project_Cortex/task_PT/subdivide_SNr/label_domain_full_K6.csv",
  sep = ",",
  row.names = FALSE,
  col.names = FALSE)


# Flatmap. preference of single-neuron projection ----

NumNeuron <- length(index_axon_PT_SNr)

proj_domain <- array(data = 0, dim = c(NumNeuron, NumCluster))
volume_domain <- array(data = 0, dim = c(NumCluster, 1))
for(iDomain in 1:NumCluster) {
  volume_domain[iDomain, 1] = length(which(label_cluster_real == iDomain))
  proj_domain[,iDomain] =
    rowSums(matrix_proj_SNr_valid[,which(label_cluster_real == iDomain)]) / volume_domain[iDomain, 1]
}

label_domain_neuron = array(data = 0, dim = c(NumNeuron, 1))
for (iNeuron in 1:NumNeuron) {
  label_domain_neuron[iNeuron, 1] <-
    which(max(proj_domain[iNeuron,]) == proj_domain[iNeuron,])[1]
}


ColorSubdomain <-
  c("1" = "#ff0000",
    "2" = "#00ff00",
    "3" = "#0000ff",
    "4" = "#ff7f00",
    "5" = "#00ffff",
    "6" = "#ff00ff",
    "7" = "#007f7f")

## flatmap whole cortex test
filename <- paste0("/tmp/PT_SNr_K4.svg")
#
U <- soma_UV$U[index_axon_PT_SNr]
V <- soma_UV$V[index_axon_PT_SNr]
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
    c(svg_element, sprintf("<circle cx=\"%f\" cy=\"%f\" r=\"%f\" style=\"stroke: none; fill: %s\"/>",
                           U[iElement], V[iElement], Radius[iElement], ColorSubdomain[label_domain_neuron[iElement]]))
}

out <- c(svg_PFC[1:50], svg_element, svg_PFC[51:length(svg_PFC)])
writeLines(out, fileConn, sep = "\n")
close(fileConn)




# region level ----

colnames(proj_domain) <- paste0("Subdomain-", 1:NumCluster)

df_neuron_domain <-
  data.frame(SourceRegion = info_tmp$Region[index_axon_PT_SNr], proj_domain)

library(dplyr)
library(tidyr)

df_region_domain <-
  df_neuron_domain %>%
  group_by(SourceRegion) %>%
  summarize_all(mean)

matrix_region_domain <- df_region_domain[,-1]
rownames(matrix_region_domain) <- df_region_domain$SourceRegion


#write.csv(x = matrix_region_domain, file = "./result_heatmap_IT.csv")


library(ComplexHeatmap)

ht_IT <-
  Heatmap(matrix_region_domain,
          cluster_rows = T,
          clustering_distance_rows = "spearman",
          clustering_method_rows = "ward.D2",
          cluster_columns = F,
          row_names_side = "left",
          column_names_side = "top")

draw(ht_IT)



# visualization. locatino of SNr ----

open3d()
shade3d(mesh_root, color = "white", alpha = 0.1)
shade3d(mesh_l_SNr, color = "red", alpha = 0.5)
rgl.viewpoint(userMatrix = userMatrix_coronal, fov = 0)


