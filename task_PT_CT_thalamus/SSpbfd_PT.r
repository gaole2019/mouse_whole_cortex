# project mouse whole cortex
# figure 5. cortico-thalamic projections
# PT neurons in SSp-bfd
# VPM, first-order; PO, high-order


# 1. select neurons with extensive projections to selected thalamic nuclei ----
index_PT_SSpbfd <-
  which(
    rowSums(shared_proj_Combined[,snp::getStructureFromAcronym(c("l_VPM", "l_PO"))]) > 500 &
      1:21582 %in% df_fnt_dist$Index[df_fnt_dist$AxonType == "PT"] &
      info_tmp$Region %in% c("SSp-bfd"))

Target_acronym <- paste0("l_", c("VPM", "PO", "RT"))


proj_PT_SSpbfd <-
  shared_proj_Combined[index_PT_SSpbfd,
                       snp::getStructureFromAcronym(Target_acronym)]

rownames(proj_PT_SSpbfd) <-
  paste0(substring(info_tmp$SampleID,2),
         "_",
         substring(info_tmp$NeuronID,2))[index_PT_SSpbfd]

colnames(proj_PT_SSpbfd) <- Target_acronym


# 2. calculate projection preference ----
PT_ratio_VPM <-
  proj_PT_SSpbfd[,"l_VPM"] / (proj_PT_SSpbfd[,"l_VPM"] + proj_PT_SSpbfd[,"l_PO"])
PT_order_ratio_VPM <- order(PT_ratio_VPM, decreasing = T)

PT_SSpbfd_label_pattern_proj <-
  vector(length = length(index_PT_SSpbfd))
for (iNeuron in 1:length(index_PT_SSpbfd)) {
  PT_SSpbfd_label_pattern_proj[iNeuron] <- "VPM and PO"
  if (PT_ratio_VPM[iNeuron] > 0.8)
    PT_SSpbfd_label_pattern_proj[iNeuron] <- "VPM prefer"
  if (PT_ratio_VPM[iNeuron] < 0.2)
    PT_SSpbfd_label_pattern_proj[iNeuron] <- "PO prefer"
}


# 3. plot heatmap ----
library(ComplexHeatmap)
ht <- Heatmap(
  matrix = t(proj_PT_SSpbfd),
  top_annotation =
    HeatmapAnnotation(
      AxonProj = PT_SSpbfd_label_pattern_proj,
      Ratio = anno_points(x = PT_ratio_VPM),
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

  column_order = PT_order_ratio_VPM,

  show_column_names = F,
  cluster_columns = F,
  use_raster = F)
draw(ht)


# 4. load meshes ----
mesh_VPM <-
  readOBJ("/home/admin/Documents/YanLab/database_neuron/ABA_mesh/733_VPM.obj")

mesh_PO <-
  readOBJ("/home/admin/Documents/YanLab/database_neuron/ABA_mesh/1020_PO.obj")

mesh_l_VPM <- getMeshLeft(mesh_VPM)
mesh_l_PO <- getMeshLeft(mesh_PO)

# 5. example neurons ----
## VPM ----
index_demo_PT_VPM <-
  index_PT_SSpbfd[which(PT_SSpbfd_label_pattern_proj == "VPM prefer" &
                       proj_PT_SSpbfd$l_VPM > 1000)]

demo_PT_VPM <-
  nat::read.neuron("/run/media/admin/LeGAO2/project_mouse_cortex/swc_files/220237/swc_allen_axon/001.swc")


open3d()
plot3d(demo_PT_VPM, WithNodes = F, color = "blue", soma = 200)

shade3d(mesh_l_VPM, color = "blue", alpha = 0.1)
shade3d(mesh_l_PO, color = "red", alpha = 0.1)
shade3d(mesh_root, color = "white", alpha = 0.1)


## VPM and PO ----

index_demo_PT_VPM_PO <-
  index_PT_SSpbfd[which(PT_SSpbfd_label_pattern_proj == "VPM and PO" &
                          proj_PT_SSpbfd$l_VPM > 3000 &
                          proj_PT_SSpbfd$l_PO > 3000)]

demo_PT_VPM_PO <-
  nat::read.neuron("/run/media/admin/LeGAO2/project_mouse_cortex/swc_files/220235/swc_allen_axon/012.swc")
demo_PT_VPM_PO$d$Z <- 11400 - demo_PT_VPM_PO$d$Z

open3d()
plot3d(demo_PT_VPM_PO, WithNodes = F, color = "blue", soma = 200)

shade3d(mesh_l_VPM, color = "blue", alpha = 0.1)
shade3d(mesh_l_PO, color = "red", alpha = 0.1)
shade3d(mesh_root, color = "white", alpha = 0.1)


## PO ----

index_demo_PT_PO <-
  index_PT_SSpbfd[which(PT_SSpbfd_label_pattern_proj == "PO prefer" &
                          proj_PT_SSpbfd$l_PO > 4000 &
                          proj_PT_SSpbfd$l_VPM == 0)]

demo_PT_PO <-
  nat::read.neuron("/run/media/admin/LeGAO2/project_mouse_cortex/swc_files/221070/swc_allen_axon/015.swc")
demo_PT_PO$d$Z <- 11400 - demo_PT_PO$d$Z

open3d()
plot3d(demo_PT_PO, WithNodes = F, color = "blue", soma = 200)

shade3d(mesh_l_VPM, color = "blue", alpha = 0.1)
shade3d(mesh_l_PO, color = "red", alpha = 0.1)
shade3d(mesh_root, color = "white", alpha = 0.1)



PT_neurite_in_PO <-
  nat::prune_in_volume(
    x = demo_PT_PO,
    surf = nat::as.hxsurf(mesh_l_PO),
    invert = T, OmitFailures = T)

open3d()
plot3d(PT_neurite_in_PO, WithNodes = F, color = "blue")

shade3d(mesh_l_VPM, color = "blue", alpha = 0.1)
shade3d(mesh_l_PO, color = "red", alpha = 0.1)


