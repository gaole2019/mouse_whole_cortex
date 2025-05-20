# project mouse whole cortex
# figure 5. cortico-thalamic projections
# CT neurons in VISp


# 1. select neurons with extensive projections to selected thalamic nuclei ----
index_CT_VISp <- which(
  rowSums(shared_proj_Combined[,snp::getStructureFromAcronym(c("l_LGd", "l_LP"))]) > 1000 &
  1:21582 %in% df_fnt_dist$Index[df_fnt_dist$AxonType == "CT"] &
  info_tmp$Region == "VISp")

Target_acronym <- paste0("l_", c("LGd", "LP", "RT"))

proj_CT_VISp <-
  shared_proj_Combined[index_CT_VISp,
                       snp::getStructureFromAcronym(Target_acronym)]

rownames(proj_CT_VISp) <-
  paste0(substring(info_tmp$SampleID,2),
         "_",
         substring(info_tmp$NeuronID,2))[index_CT_VISp]

colnames(proj_CT_VISp) <- Target_acronym


# 2. calculate projection preference ----
CT_ratio_LGd <-
  proj_CT_VISp[,"l_LGd"] / (proj_CT_VISp[,"l_LGd"] + proj_CT_VISp[,"l_LP"])
order_ratio_LGd <- order(ratio_LGd, decreasing = T)

CT_VISp_label_pattern_proj <- vector(length = length(index_CT_VISp))
for (iNeuron in 1:length(index_CT_VISp)) {
  CT_VISp_label_pattern_proj[iNeuron] <- "LGd and LP"
  if (CT_ratio_LGd[iNeuron] > 0.736)
    CT_VISp_label_pattern_proj[iNeuron, 1] <- "LGd prefer"
  if (CT_ratio_LGd[iNeuron] < 0.219)
    CT_VISp_label_pattern_proj[iNeuron, 1] <- "LP prefer"
}


# 3. plot heatmap ----
library(ComplexHeatmap)
ht <- Heatmap(
  matrix = t(proj_CT_VISp),
  top_annotation =
    HeatmapAnnotation(
      AxonProj = CT_VISp_label_pattern_proj,
      Ratio = anno_points(x = CT_ratio_LGd),
      col = list(AxonProj = c("LGd prefer" = "red",
                              "LGd and LP" = "green",
                              "LP prefer" = "blue"))
    ),
  name = "Axon length (um)",
  heatmap_legend_param	=
    list(title_gp = gpar(fontsize = 8),
         labels_gp = gpar(fontsize = 8),
         title_position = "leftcenter-rot"),

  cluster_rows = F,
  row_names_side = "left",
  row_names_rot = 0,
  row_names_gp = gpar(fontsize = 8),

  column_order = CT_order_ratio_LGd,

  show_column_names = F,
  cluster_columns = F,
  use_raster = F)
draw(ht)


# 4. load mesh ----
mesh_LGd <- readOBJ("/home/admin/Documents/YanLab/database_neuron/ABA_mesh/170_LGd.obj")

mesh_LP <- readOBJ("/home/admin/Documents/YanLab/database_neuron/ABA_mesh/218_LP.obj")

mesh_l_LGd <- getMeshLeft(mesh_LGd)
mesh_l_LP <- getMeshLeft(mesh_LP)

# 5. example neurons ----
## demo: LGd-projecting neurons ----
index_full_CT_VISp_LGd <-
  index_full_CT_VISp[which(VISp_label_pattern_proj == "LGd prefer")][11]


nl_axon_CT_VISp_LGd <-
  nat::read.neurons(
    paths = shared_SWCPath_axon_allen[index_full_CT_VISp_LGd],
    neuronnames = shared_NeuronName[index_full_CT_VISp_LGd])

nl_dendrite_CT_VISp_LGd <-
  nat::read.neurons(
    paths = shared_SWCPath_dendrite_allen[index_full_CT_VISp_LGd],
    neuronnames = shared_NeuronName[index_full_CT_VISp_LGd])

open3d()
snp::plot_morphology(
  neuronList = nl_axon_CT_VISp_LGd,
  hemisphere = info_tmp$Hemisphere[index_full_CT_VISp_LGd],
  flip = TRUE, soma_size = 100, color = "blue")

snp::plot_morphology(
  neuronList = nl_dendrite_CT_VISp_LGd,
  hemisphere = info_tmp$Hemisphere[index_full_CT_VISp_LGd],
  flip = TRUE, soma_size = 100, color = "red")



shade3d(mesh_l_LGd, color = "blue", alpha = 0.1)

shade3d(mesh_l_LP, color = "red", alpha = 0.1)
shade3d(mesh_root, color = "white", alpha = 0.1)

rgl.viewpoint(fov = 0, userMatrix = userMatrix_coronal)




## demo: LGd- and LP-projecting neurons ----

index_full_CT_VISp_LGd_LP <-
  index_full_CT_VISp[which(VISp_label_pattern_proj == "LGd and LP")][1]


nl_axon_CT_VISp_LGd_LP <-
  nat::read.neurons(
    paths = shared_SWCPath_axon_allen[index_full_CT_VISp_LGd_LP],
    neuronnames = shared_NeuronName[index_full_CT_VISp_LGd_LP])

nl_dendrite_CT_VISp_LGd_LP <-
  nat::read.neurons(
    paths = shared_SWCPath_dendrite_allen[index_full_CT_VISp_LGd_LP],
    neuronnames = shared_NeuronName[index_full_CT_VISp_LGd_LP])

open3d()
snp::plot_morphology(
  neuronList = nl_axon_CT_VISp_LGd_LP,
  hemisphere = info_tmp$Hemisphere[index_full_CT_VISp_LGd_LP],
  flip = TRUE, soma_size = 100, color = "blue")

snp::plot_morphology(
  neuronList = nl_dendrite_CT_VISp_LGd_LP,
  hemisphere = info_tmp$Hemisphere[index_full_CT_VISp_LGd_LP],
  flip = TRUE, soma_size = 100, color = "red")



shade3d(mesh_l_LGd, color = "blue", alpha = 0.1)

shade3d(mesh_l_LP, color = "red", alpha = 0.1)
#shade3d(mesh_root, color = "white", alpha = 0.1)

rgl.viewpoint(fov = 0, userMatrix = userMatrix_coronal)






## demo: LP-projecting neurons ----

index_full_CT_VISp_LP <-
  index_full_CT_VISp[which(VISp_label_pattern_proj == "LP prefer")][6]


nl_axon_CT_VISp_LP <-
  nat::read.neurons(
    paths = shared_SWCPath_axon_allen[index_full_CT_VISp_LP],
    neuronnames = shared_NeuronName[index_full_CT_VISp_LP])

nl_dendrite_CT_VISp_LP <-
  nat::read.neurons(
    paths = shared_SWCPath_dendrite_allen[index_full_CT_VISp_LP],
    neuronnames = shared_NeuronName[index_full_CT_VISp_LP])

open3d()
snp::plot_morphology(
  neuronList = nl_axon_CT_VISp_LP,
  hemisphere = info_tmp$Hemisphere[index_full_CT_VISp_LP],
  flip = TRUE, soma_size = 100, color = "blue")

snp::plot_morphology(
  neuronList = nl_dendrite_CT_VISp_LP,
  hemisphere = info_tmp$Hemisphere[index_full_CT_VISp_LP],
  flip = TRUE, soma_size = 100, color = "red")



shade3d(mesh_l_LGd, color = "blue", alpha = 0.1)

shade3d(mesh_l_LP, color = "red", alpha = 0.1)
#shade3d(mesh_root, color = "white", alpha = 0.1)

rgl.viewpoint(fov = 0, userMatrix = userMatrix_coronal)

