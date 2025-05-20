# project mouse whole cortex
# figure 5. cortico-thalamic projections
# CT neurons in AUD


# 1. select neurons with extensive projections to selected thalamic nuclei ----
index_CT_AUD <- which(
  rowSums(shared_proj_Combined[,snp::getStructureFromAcronym(c("l_MGv", "l_MGd"))]) > 500 &
  1:21582 %in% df_fnt_dist$Index[df_fnt_dist$AxonType == "CT"] &
  info_tmp$Region %in% c("AUDd", "AUDp", "AUDpo", "AUDv"))


Target_acronym <- paste0("l_", c("MGv", "MGd", "RT"))

proj_CT_AUD <-
  shared_proj_Combined[index_CT_AUD,
                       snp::getStructureFromAcronym(Target_acronym)]

rownames(proj_CT_AUD) <-
  paste0(substring(info_tmp$SampleID,2),
         "_",
         substring(info_tmp$NeuronID,2))[index_CT_AUD]

colnames(proj_CT_AUD) <- Target_acronym


# 2. calculate their projection preference ----
CT_ratio_MGv <-
  proj_CT_AUD[,"l_MGv"] / (proj_CT_AUD[,"l_MGv"] + proj_CT_AUD[,"l_MGd"])

CT_order_ratio_MGv <- order(CT_ratio_MGv, decreasing = T)


CT_AUD_label_pattern_proj <-
  vector(length = length(index_CT_AUD))
for (iNeuron in 1:length(index_CT_AUD)) {
  CT_AUD_label_pattern_proj[iNeuron] <- "MGv and MGd"
  if (CT_ratio_MGv[iNeuron] > 0.8)
    CT_AUD_label_pattern_proj[iNeuron] <- "MGv prefer"
  if (CT_ratio_MGv[iNeuron] < 0.2)
    CT_AUD_label_pattern_proj[iNeuron] <- "MGd prefer"
}



# 3. heatmap ----
library(ComplexHeatmap)
ht <- Heatmap(
  matrix = t(proj_CT_AUD),
  top_annotation =
    HeatmapAnnotation(
      AxonProj = CT_AUD_label_pattern_proj,

      Ratio = anno_points(x = CT_ratio_MGv),


      col = list(AxonProj = c("MGv prefer" = "red",
                              "MGv and MGd" = "green",
                              "MGd prefer" = "blue"))
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

  column_order = CT_order_ratio_MGv,

  show_column_names = F,
  cluster_columns = F,
  use_raster = F)
draw(ht)




# 4. load meshes ----
mesh_lr_MGd <- readOBJ("/home/admin/Documents/YanLab/database_neuron/ABA_mesh/1072_MGd.obj")

mesh_lr_MGv <- readOBJ("/home/admin/Documents/YanLab/database_neuron/ABA_mesh/1079_MGv.obj")

mesh_l_MGd <- getMeshLeft(mesh_lr_MGd)
mesh_l_MGv <- getMeshLeft(mesh_lr_MGv)


# 5. demo neurons ----
## MGd-projecting neurons ----

index_full_CT_AUD_MGd <-
  index_full_CT_AUD[which(AUD_label_pattern_proj == "MGd prefer")][3]


nl_axon_CT_AUD_MGd <-
  nat::read.neurons(
    paths = shared_SWCPath_axon_allen[index_full_CT_AUD_MGd],
    neuronnames = shared_NeuronName[index_full_CT_AUD_MGd])

nl_dendrite_CT_AUD_MGd <-
  nat::read.neurons(
    paths = shared_SWCPath_dendrite_allen[index_full_CT_AUD_MGd],
    neuronnames = shared_NeuronName[index_full_CT_AUD_MGd])

open3d()
snp::plot_morphology(
  neuronList = nl_axon_CT_AUD_MGd,
  hemisphere = info_tmp$Hemisphere[index_full_CT_AUD_MGd],
  flip = TRUE, soma_size = 100, color = "blue")

snp::plot_morphology(
  neuronList = nl_dendrite_CT_AUD_MGd,
  hemisphere = info_tmp$Hemisphere[index_full_CT_AUD_MGd],
  flip = TRUE, soma_size = 100, color = "red")



shade3d(mesh_l_MGd, color = "blue", alpha = 0.1)
shade3d(mesh_l_MGv, color = "red", alpha = 0.1)
#shade3d(mesh_root, color = "white", alpha = 0.1)

rgl.viewpoint(fov = 0, userMatrix = userMatrix_sagittal)



## demo: MGd- and MGv-projecting neurons ----

index_full_CT_AUD_MGd_MGv <-
  index_full_CT_AUD[which(AUD_label_pattern_proj == "MGd and MGv")][1]


nl_axon_CT_AUD_MGd_MGv <-
  nat::read.neurons(
    paths = shared_SWCPath_axon_allen[index_full_CT_AUD_MGd_MGv],
    neuronnames = shared_NeuronName[index_full_CT_AUD_MGd_MGv])

nl_dendrite_CT_AUD_MGd_MGv <-
  nat::read.neurons(
    paths = shared_SWCPath_dendrite_allen[index_full_CT_AUD_MGd_MGv],
    neuronnames = shared_NeuronName[index_full_CT_AUD_MGd_MGv])

open3d()
snp::plot_morphology(
  neuronList = nl_axon_CT_AUD_MGd_MGv,
  hemisphere = info_tmp$Hemisphere[index_full_CT_AUD_MGd_MGv],
  flip = TRUE, soma_size = 100, color = "blue")

snp::plot_morphology(
  neuronList = nl_dendrite_CT_AUD_MGd_MGv,
  hemisphere = info_tmp$Hemisphere[index_full_CT_AUD_MGd_MGv],
  flip = TRUE, soma_size = 100, color = "red")



shade3d(mesh_l_MGd, color = "blue", alpha = 0.1)
shade3d(mesh_l_MGv, color = "red", alpha = 0.1)
shade3d(mesh_root, color = "white", alpha = 0.1)

rgl.viewpoint(fov = 0, userMatrix = userMatrix_sagittal)



## demo: MGv-projecting neurons ----

index_full_CT_AUD_MGv <-
  index_full_CT_AUD[which(AUD_label_pattern_proj == "MGv prefer")][3]


nl_axon_CT_AUD_MGv <-
  nat::read.neurons(
    paths = shared_SWCPath_axon_allen[index_full_CT_AUD_MGv],
    neuronnames = shared_NeuronName[index_full_CT_AUD_MGv])

nl_dendrite_CT_AUD_MGv <-
  nat::read.neurons(
    paths = shared_SWCPath_dendrite_allen[index_full_CT_AUD_MGv],
    neuronnames = shared_NeuronName[index_full_CT_AUD_MGv])

open3d()
snp::plot_morphology(
  neuronList = nl_axon_CT_AUD_MGv,
  hemisphere = info_tmp$Hemisphere[index_full_CT_AUD_MGv],
  flip = TRUE, soma_size = 100, color = "blue")

snp::plot_morphology(
  neuronList = nl_dendrite_CT_AUD_MGv,
  hemisphere = info_tmp$Hemisphere[index_full_CT_AUD_MGv],
  flip = TRUE, soma_size = 100, color = "red")



shade3d(mesh_l_MGd, color = "blue", alpha = 0.1)
shade3d(mesh_l_MGv, color = "red", alpha = 0.1)
shade3d(mesh_root, color = "white", alpha = 0.1)

rgl.viewpoint(fov = 0, userMatrix = userMatrix_sagittal)

