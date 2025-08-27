# project mouse whole cortex
# figure 5. cortico-thalamic projections
# PT and CT
# subdivide RT


# 1. select neurons projecting to ipsilateral RT ----
index_RT_PT <-
  which(shared_proj_Combined[, snp::getStructureFromAcronym("l_RT")] > 100 &
          1:21582 %in% df_fnt_dist$Index[which(df_fnt_dist$AxonType == "PT")])

index_RT_CT <-
  which(shared_proj_Combined[, snp::getStructureFromAcronym("l_RT")] > 100 &
          1:21582 %in% df_fnt_dist$Index[which(df_fnt_dist$AxonType == "CT")])

index_RT_PT_CT <- c(index_RT_PT, index_RT_CT)

Hemisphere_Path_PT_CT_RT <-
  data.frame(Hemisphere = info_tmp$Hemisphere[index_RT_PT_CT],
             Path = shared_SWCPath_axon_allen[index_RT_PT_CT])

write.table(
  x = Hemisphere_Path_PT_CT_RT,
  file = "/run/media/admin/LeGAO2/project_mouse_cortex/analysis/task_PT_CT_subdivide_RT/Hemisphere_Path_PT_CT_RT2.csv",
  row.names = FALSE,
  col.names = FALSE
)


# 2. prepare anno_RT ----

# 3. calculate the axon length in each cube ----


# 4. clustering ----
matrix_proj_RT <-
  read.csv(file = "/run/media/admin/LeGAO/Project_Cortex/task_PT_CT_subdivide_RT/axonLength.csv",
           row.names = 1)
matrix_proj_RT <- as.matrix(matrix_proj_RT)


ff_similarity <-
  snp::snp_bigcorPar(
    x = matrix_proj_RT,
    size = 400,
    verbose = TRUE,
    ncore = 4,
    method = "spearman"
  )

matrix_similarity = as.matrix(as.data.frame(ff_similarity[, ]))


#matrix_similarity = read.table("./matrix_similarity_800.txt")


library(fastcluster)

h <-
  fastcluster::hclust(d = as.dist(1 - matrix_similarity), method = "ward.D2")

plot(h, hang = -1, labels = FALSE)
rect.hclust(tree = h, k = 7)



NumCluster <- 7
label_cluster_cube <- cutree(tree = h, k = NumCluster)

labelOnDend <- unique(label_cluster_cube[h$order])

label_cluster_real <- label_cluster_cube
for (iLabel in 1:NumCluster) {
  label_cluster_real[which(label_cluster_cube == labelOnDend[iLabel])] <-
    iLabel
}

write.table(
  x = label_cluster_real,
  file = "/run/media/admin/LeGAO2/project_mouse_cortex/analysis/task_PT_CT_subdivide_RT/K7/label_cluster_real_K7.csv",
  sep = ",",
  row.names = FALSE,
  col.names = FALSE
)



# plot heatmap of sperman's r ----
matrix_similarity_ordered <- matrix_similarity[h$order, h$order]

col_fun <- circlize::colorRamp2(breaks = c(-0.05, 0, 0.1),
                                colors = c("blue", "white", "red"))

library(ComplexHeatmap)
ht <-
  Heatmap(
    matrix = matrix_similarity_ordered[seq(1, 5793, 5), seq(1, 5793, 5)],
    name = "Spearman's r",
    col = col_fun,
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
NumNeuron <- length(index_RT_PT_CT)

proj_domain <- array(data = 0, dim = c(NumNeuron, NumCluster))
volume_domain <- array(data = 0, dim = c(NumCluster, 1))
for (iDomain in 1:NumCluster) {
  volume_domain[iDomain, 1] <-
    length(which(label_cluster_real == iDomain))
  proj_domain[, iDomain] =
    rowSums(matrix_proj_RT[, which(label_cluster_real == iDomain)]) / volume_domain[iDomain, 1]
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

filename <- paste0("/tmp/PT_CT_RT_K7.svg")
#
U <- soma_UV$U[index_axon_PT_CT_RT]
V <- soma_UV$V[index_axon_PT_CT_RT]
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
colnames(proj_domain) <- paste0("Subdomain-", 1:7)

df_neuron_domain <-
  data.frame(SourceRegion = info_tmp$Region[index_RT_PT_CT], proj_domain)

library(dplyr)
library(tidyr)

df_region_domain <-
  df_neuron_domain %>%
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
