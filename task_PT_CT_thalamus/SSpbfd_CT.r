# project mouse whole cortex
# figure 5. cortico-thalamic projections
# CT neurons in SSp-bfd


# 1. select neurons with extensive projections to selected thalamic nuclei ----
index_CT_SSpbfd <-
  which(
    rowSums(shared_proj_Combined[,snp::getStructureFromAcronym(c("l_VPM", "l_PO"))]) > 500 &
      1:21582 %in% df_fnt_dist$Index[df_fnt_dist$AxonType == "CT"] &
      info_tmp$Region %in% c("SSp-bfd"))

Target_acronym <- paste0("l_", c("VPM", "PO", "RT"))


proj_CT_SSpbfd <-
  shared_proj_Combined[index_CT_SSpbfd,
                       snp::getStructureFromAcronym(Target_acronym)]

rownames(proj_CT_SSpbfd) <-
  paste0(substring(info_tmp$SampleID,2),
         "_",
         substring(info_tmp$NeuronID,2))[index_CT_SSpbfd]

colnames(proj_CT_SSpbfd) <- Target_acronym


# 2. calculate their projection preference ----
CT_ratio_VPM <-
  proj_CT_SSpbfd[,"l_VPM"] / (proj_CT_SSpbfd[,"l_VPM"] + proj_CT_SSpbfd[,"l_PO"])
CT_order_ratio_VPM <- order(CT_ratio_VPM, decreasing = T)


CT_SSpbfd_label_pattern_proj <-
  vector(length = length(index_CT_SSpbfd))
for (iNeuron in 1:length(index_CT_SSpbfd)) {
  CT_SSpbfd_label_pattern_proj[iNeuron] <- "VPM and PO"
  if (CT_ratio_VPM[iNeuron] > 0.8)
    CT_SSpbfd_label_pattern_proj[iNeuron] <- "VPM prefer"
  if (CT_ratio_VPM[iNeuron] < 0.2)
    CT_SSpbfd_label_pattern_proj[iNeuron] <- "PO prefer"
}


# 3. plot heatmap ----
library(ComplexHeatmap)
ht <- Heatmap(
  matrix = t(proj_CT_SSpbfd),
  top_annotation =
    HeatmapAnnotation(
      AxonProj = CT_SSpbfd_label_pattern_proj,
      Ratio = anno_points(x = CT_ratio_VPM),
      col = list(AxonProj = c("VPM prefer" = "red",
                              "VPM and PO" = "green",
                              "PO prefer" = "blue"))
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

  column_order = CT_order_ratio_VPM,

  show_column_names = F,
  cluster_columns = F,
  use_raster = F)
draw(ht)



# 4. load meshes ----

mesh_VPM <-
  readOBJ("/home/admin/Documents/YanLab/database_neuron/ABA_mesh/733_VPM.obj")
mesh_l_VPM <- getMeshLeft(mesh_VPM)

mesh_RT <-
  readOBJ("/home/admin/Documents/YanLab/database_neuron/ABA_mesh/262_RT.obj")
mesh_l_RT <- getMeshLeft(mesh_RT)

mesh_PO <-
  readOBJ("/home/admin/Documents/YanLab/database_neuron/ABA_mesh/1020_PO.obj")
mesh_l_PO <- getMeshLeft(mesh_PO)


mesh_root <- readOBJ("/home/admin/Documents/YanLab/database_neuron/ABA_mesh/997_root.obj")


# 5. plot example neurons ----
## demo: VPM-projecting neurons ----

index_full_CT_SSpbfd_VPM <-
  index_full_CT_SSpbfd[which(label_pattern_proj == "VPM prefer")][13]


nl_axon_CT_SSbfd_VPM <-
  nat::read.neurons(paths = shared_SWCPath_axon_allen[index_full_CT_SSpbfd_VPM],
                    neuronnames = shared_NeuronName[index_full_CT_SSpbfd_VPM])

nl_dendrite_CT_SSbfd_VPM <-
  nat::read.neurons(paths = shared_SWCPath_dendrite_allen[index_full_CT_SSpbfd_VPM],
                    neuronnames = shared_NeuronName[index_full_CT_SSpbfd_VPM])

open3d()
snp::plot_morphology(
  neuronList = nl_axon_CT_SSbfd_VPM,
  hemisphere = info_tmp$Hemisphere[index_full_CT_SSpbfd_VPM],
  flip = TRUE, soma_size = 100, color = "blue")

snp::plot_morphology(
  neuronList = nl_dendrite_CT_SSbfd_VPM,
  hemisphere = info_tmp$Hemisphere[index_full_CT_SSpbfd_VPM],
  flip = TRUE, soma_size = 100, color = "red")


shade3d(mesh_root, color = "white", alpha = 0.1)

shade3d(mesh_l_VPM, color = "blue", alpha = 0.1)
#shade3d(mesh_lr_RT, color = "blue", alpha = 0.1)

shade3d(mesh_l_PO, color = "red", alpha = 0.1)

rgl.viewpoint(fov = 0, userMatrix = userMatrix_coronal)



## demo: VPM- and PO- projecting neurons ----
index_full_CT_SSpbfd_VPM_PO <-
  index_full_CT_SSpbfd[which(label_pattern_proj == "VPM and PO")][5]


nl_axon_CT_SSbfd_VPM_PO <-
  nat::read.neurons(
    paths = shared_SWCPath_axon_allen[index_full_CT_SSpbfd_VPM_PO],
    neuronnames = shared_NeuronName[index_full_CT_SSpbfd_VPM_PO])

nl_dendrite_CT_SSbfd_VPM_PO <-
  nat::read.neurons(
    paths = shared_SWCPath_dendrite_allen[index_full_CT_SSpbfd_VPM_PO],
    neuronnames = shared_NeuronName[index_full_CT_SSpbfd_VPM_PO])

open3d()
snp::plot_morphology(
  neuronList = nl_axon_CT_SSbfd_VPM_PO,
  hemisphere = info_tmp$Hemisphere[index_full_CT_SSpbfd_VPM_PO],
  flip = TRUE, soma_size = 100, color = "blue")

snp::plot_morphology(
  neuronList = nl_dendrite_CT_SSbfd_VPM_PO,
  hemisphere = info_tmp$Hemisphere[index_full_CT_SSpbfd_VPM_PO],
  flip = TRUE, soma_size = 100, color = "red")


#shade3d(mesh_root, color = "white", alpha = 0.1)

shade3d(mesh_l_VPM, color = "blue", alpha = 0.1)
#shade3d(mesh_lr_RT, color = "blue", alpha = 0.1)

shade3d(mesh_l_PO, color = "red", alpha = 0.1)

rgl.viewpoint(fov = 0, userMatrix = userMatrix_coronal)


##  demo: PO-projecting neurons ----
index_full_CT_SSpbfd_PO <-
  index_full_CT_SSpbfd[which(label_pattern_proj == "PO prefer")][11]


nl_axon_CT_SSbfd_PO <-
  nat::read.neurons(
    paths = shared_SWCPath_axon_allen[index_full_CT_SSpbfd_PO],
    neuronnames = shared_NeuronName[index_full_CT_SSpbfd_PO])

nl_dendrite_CT_SSbfd_PO <-
  nat::read.neurons(
    paths = shared_SWCPath_dendrite_allen[index_full_CT_SSpbfd_PO],
    neuronnames = shared_NeuronName[index_full_CT_SSpbfd_PO])

open3d()
snp::plot_morphology(
  neuronList = nl_axon_CT_SSbfd_PO,
  hemisphere = info_tmp$Hemisphere[index_full_CT_SSpbfd_PO],
  flip = TRUE, soma_size = 100, color = "blue")

snp::plot_morphology(
  neuronList = nl_dendrite_CT_SSbfd_PO,
  hemisphere = info_tmp$Hemisphere[index_full_CT_SSpbfd_PO],
  flip = TRUE, soma_size = 100, color = "red")



shade3d(mesh_l_VPM, color = "blue", alpha = 0.1)

shade3d(mesh_l_PO, color = "red", alpha = 0.1)
shade3d(mesh_root, color = "white", alpha = 0.1)

rgl.viewpoint(fov = 0, userMatrix = userMatrix_coronal)

