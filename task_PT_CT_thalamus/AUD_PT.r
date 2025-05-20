# project mouse whole cortex
# figure 5. cortico-thalamic projections
# PT neurons in AUD
# MGv, first-order; MGd, high-order

# 1. select neurons with extensive projections to selected thalamic nuclei ----
index_PT_AUD <-
  which(
    rowSums(shared_proj_Combined[,snp::getStructureFromAcronym(c("l_MGv", "l_MGd"))]) > 500 &
    1:21582 %in% df_fnt_dist$Index[df_fnt_dist$AxonType == "PT"] &
    info_tmp$Region %in% c("AUDd", "AUDp", "AUDpo", "AUDv"))

Target_acronym <- paste0("l_", c("MGv", "MGd", "RT"))


proj_PT_AUD <-
  shared_proj_Combined[index_PT_AUD,
                       snp::getStructureFromAcronym(Target_acronym)]

rownames(proj_PT_AUD) <-
  paste0(substring(info_tmp$SampleID,2),
         "_",
         substring(info_tmp$NeuronID,2))[index_PT_AUD]

colnames(proj_PT_AUD) <- Target_acronym


# 2. calculate their projection preference ----
PT_ratio_MGv <-
  proj_PT_AUD[,"l_MGv"] / (proj_PT_AUD[,"l_MGv"] + proj_PT_AUD[,"l_MGd"])
PT_order_ratio_MGv <- order(PT_ratio_MGv, decreasing = T)

PT_AUD_label_pattern_proj <- vector(length = length(index_PT_AUD))
for (iNeuron in 1:length(index_PT_AUD)) {
  PT_AUD_label_pattern_proj[iNeuron] <- "MGv and MGd"
  if (PT_ratio_MGv[iNeuron] > 0.8)
    PT_AUD_label_pattern_proj[iNeuron] <- "MGv prefer"
  if (PT_ratio_MGv[iNeuron] < 0.2)
    PT_AUD_label_pattern_proj[iNeuron] <- "MGd prefer"
}


# 3. plot heatmap ----
library(ComplexHeatmap)
ht <- Heatmap(
  matrix = t(proj_PT_AUD),
  top_annotation =
    HeatmapAnnotation(
      AxonProj = PT_AUD_label_pattern_proj,
      Ratio = anno_points(x = PT_ratio_MGv),
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

  column_order = PT_order_ratio_MGv,

  show_column_names = F,
  cluster_columns = F,
  use_raster = F)
draw(ht)


# 4. load meshes ----

mesh_MGd <-
  readOBJ("/home/admin/Documents/YanLab/database_neuron/ABA_mesh/1072_MGd.obj")

mesh_MGv <-
  readOBJ("/home/admin/Documents/YanLab/database_neuron/ABA_mesh/1079_MGv.obj")

mesh_l_MGd <- getMeshLeft(mesh_MGd)
mesh_l_MGv <- getMeshLeft(mesh_MGv)

# 5. example PT neurons ----
## MGv ----
index_demo_PT_MGv <-
  index_PT_AUD[which(PT_AUD_label_pattern_proj == "MGv prefer" &
                     proj_PT_AUD$l_MGv > 1000)]

demo_PT_MGv <-
  nat::read.neuron("/run/media/admin/LeGAO2/project_mouse_cortex/swc_files/221530/swc_allen_axon/061.swc")

demo_PT_MGv$d$Z <- 11400 - demo_PT_MGv$d$Z


open3d()
plot3d(demo_PT_MGv, WithNodes = F, color = "blue")

shade3d(mesh_l_MGd, color = "blue", alpha = 0.1)
shade3d(mesh_l_MGv, color = "red", alpha = 0.1)
shade3d(mesh_root, color = "white", alpha = 0.1)



## MGv and MGd ----
index_demo_PT_MGv_MGd <-
  index_PT_AUD[which(PT_AUD_label_pattern_proj == "MGv and MGd" &
                       proj_PT_AUD$l_MGv > 1000 &
                       proj_PT_AUD$l_MGd > 1000)]

demo_PT_MGv_MGd <-
  nat::read.neuron("/run/media/admin/LeGAO2/project_mouse_cortex/swc_files/221474/swc_allen_axon/155.swc")

demo_PT_MGv_MGd$d$Z <- 11400 - demo_PT_MGv_MGd$d$Z


open3d()
plot3d(demo_PT_MGv_MGd, WithNodes = F, color = "blue")

shade3d(mesh_l_MGd, color = "blue", alpha = 0.1)
shade3d(mesh_l_MGv, color = "red", alpha = 0.1)
shade3d(mesh_root, color = "white", alpha = 0.1)


## MGd ----



