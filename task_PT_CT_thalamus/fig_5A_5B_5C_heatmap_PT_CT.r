# project mouse whole cortex
# figure 5. cortico-thalamic projection
# PT and CT


# 1. 44 thalamic targets ----
target_acronym <-
  c("VAL", "VM", "VPL", "VPLpc", "VPM", "VPMpc", "SPFm", "SPFp", "SPA", "PP", "MG", "LGd", # DORsm
    "LP", "PO", "POL", "SGN", "Eth", # LAT
    "AV", "AM", "AD", "IAM", "IAD", "LD", # ATN
    "IMD", "MD", "SMT", "PR", #MED
    "PVT", "PT", "RE", "Xi", # MTN
    "RH", "CM", "PCN", "CL", "PF", "PIL",# ILM
    "RT",
    "IGL", "IntG", "LGv", "SubG", # GENv
    "MH", "LH" #EPI
  )

target_acronym <- paste0("l_", target_acronym)



index_axon_PT_CT <- which(
  1:21582 %in% df_All_index_AxonType$Index[df_All_index_AxonType$AxonType %in% c("PT", "CT")] &
  rowSums(shared_proj_Combined[,snp::getStructureFromAcronym(target_acronym)]) > 10000)


proj_PT_CT_TH <-
  shared_proj_Combined[index_axon_PT_CT,
                       snp::getStructureFromAcronym(target_acronym)]

colnames(proj_PT_CT_TH) <- target_acronym



map_level1_level2 <-
  read.csv(file = "/run/media/admin/LeGAO2/project_mouse_cortex/analysis/shared/map_level1_level2.csv",
           header = T)


label_region_L2 <-
  map_level1_level2$Region2[match(info_tmp$Region[index_axon_PT_CT],
                                  map_level1_level2$Region)]



# 3. color map of region level 2 ----
color_K17 <- read.csv("~/shared/colors_17.txt")
color_K17 <- rgb(color_K17$R, color_K17$G, color_K17$B, maxColorValue = 255)

names(color_K17) <-
  c(unique(map_level1_level2$Region2), "NA", "NA2", "NA3" ,"NA4")

col <- list()
col[["Region"]] <- color_K17
col[["AxonType"]] <- c("PT" = "blue", "CT" = "red")

library(ComplexHeatmap)

# 4. clustering: inter-region ----
ht <-
  Heatmap(
    matrix = t(proj_PT_CT_TH),
    top_annotation =
      HeatmapAnnotation(
        Region = label_region_L2,
        AxonType = df_All_index_AxonType$AxonType[match(index_axon_PT_CT, df_All_index_AxonType$Index)],
        col = col),
    column_split = label_region_L2,
    clustering_distance_rows = "spearman",
    clustering_method_rows = "ward.D2",

    cluster_column_slices = T,

    clustering_distance_columns = "spearman",
    clustering_method_columns = "ward.D2",
    show_column_names = F,

    show_row_names = T,
    row_names_side = "left",
    use_raster = T)

draw(ht)



region_order <- names(column_order(ht))



# 5. clustering: intra-region ----
ht <-
  Heatmap(
    matrix = t(proj_PT_CT_TH),
    top_annotation =
      HeatmapAnnotation(
        Region = label_region_L2,
        AxonType = df_All_index_AxonType$AxonType[match(index_axon_PT_CT, df_All_index_AxonType$Index)],
        col = col),
    column_split = factor(label_region_L2, level = region_order),
    clustering_distance_rows = "spearman",
    clustering_method_rows = "ward.D2",

    cluster_column_slices = F,

    clustering_distance_columns = "spearman",
    clustering_method_columns = "ward.D2",
    show_column_names = F,

    show_row_names = T,
    row_names_side = "left",
    use_raster = T)

draw(ht)


# 6. PT vs. CT neurons in RSP ----
df_PT_CT_RSP <-
  data.frame(
    proj_PT_CT_TH[label_region_L2 == "RSP", c("RT", "LGd", "LGv")],
    AxonType = info_tmp$AxonType[index_axon_PT_CT][label_region_L2 == "RSP"]) %>%
  pivot_longer(cols = 1:3, names_to = "Region", values_to = "Length")

df_PT_CT_RSP$Region <- factor(x = df_PT_CT_RSP$Region, levels = c("RT", "LGd", "LGv"))



library(ggpubr)
ggplot(data = df_PT_CT_RSP, mapping = aes(x = Region, y = Length, fill = AxonType)) +
  geom_boxplot() +
  #geom_violin() +
  #geom_jitter(width = 0.1) +
  theme_classic() +stat_compare_means() +
  ggtitle("PT vs. CT. RSP")


## read RSP neurons ----

library(nat)
nl_neurite_PT_CT_RSP <-
  nat::read.neurons(
    paths = SWCPath_neurite[index_axon_PT_CT[label_region_L2 == "RSP"]],
    neuronnames = paste0(info_tmp$SampleID, "-", info_tmp$NeuronID)[index_axon_PT_CT[label_region_L2 == "RSP"]])

hemisphere <- info_tmp$Hemisphere[index_axon_PT_CT[label_region_L2 == "RSP"]]


## get index of PT ----
which(proj_PT_CT_TH$RT[label_region_L2 == "RSP"] == 0 &
        proj_PT_CT_TH$LGd[label_region_L2 == "RSP"] > 800 &
        proj_PT_CT_TH$LGv[label_region_L2 == "RSP"] > 1000 &
        info_tmp$Region[index_axon_PT_CT][label_region_L2 == "RSP"] == "RSPd" &
        info_tmp$AxonType[index_axon_PT_CT][label_region_L2 == "RSP"] == "PT")

index_tmp_RSP_PT <- c(40,45)



## get index of IT ----
which(proj_PT_CT_TH$RT[label_region_L2 == "RSP"] > 2000 &
        proj_PT_CT_TH$LGd[label_region_L2 == "RSP"] == 0 &
        proj_PT_CT_TH$LGv[label_region_L2 == "RSP"] == 0 &
        info_tmp$Region[index_axon_PT_CT][label_region_L2 == "RSP"] == "RSPd" &
        info_tmp$AxonType[index_axon_PT_CT][label_region_L2 == "RSP"] == "CT")

index_tmp_RSP_CT <- c(128, 129)

open3d()
snp::plot_morphology(
  neuronList = nl_neurite_PT_CT_RSP[index_tmp_RSP_PT[1]],
  hemisphere = hemisphere[index_tmp_RSP_PT[1]],
  color = "blue")

snp::plot_morphology(
  neuronList = nl_neurite_PT_CT_RSP[index_tmp_RSP_CT[1]],
  hemisphere = hemisphere[index_tmp_RSP_CT[1]],
  color = "red")

## read mesh of RT, LGd, LGv ----

mesh_RT <- readOBJ("~/Documents/YanLab/database_neuron/ABA_mesh/262_RT.obj")
mesh_l_RT <- getMeshLeft(mesh_RT)

mesh_LGd <- readOBJ("~/Documents/YanLab/database_neuron/ABA_mesh/170_LGd.obj")
mesh_l_LGd <- getMeshLeft(mesh_LGd)

mesh_LGv <- readOBJ("~/Documents/YanLab/database_neuron/ABA_mesh/178_LGv.obj")
mesh_l_LGv <- getMeshLeft(mesh_LGv)


shade3d(mesh_l_LGd, color = "red", alpha = 0.1)
shade3d(mesh_l_LGv, color = "red", alpha = 0.1)

shade3d(mesh_l_RT, color = "green", alpha = 0.1)
shade3d(mesh_root, color = "white", alpha = 0.1)

rgl::view3d(userMatrix = userMatrix_horizontal, fov = 0)




# 7. PT vs. CT. PL ----
df_PT_CT_PL <-
  data.frame(
    proj_PT_CT_TH[label_region_L2 == "PL", c("MD", "RT", "PF", "SPA")],
    AxonType = info_tmp$AxonType[index_axon_PT_CT][label_region_L2 == "PL"]) %>%
  pivot_longer(cols = c("MD", "RT", "PF", "SPA"),
               names_to = "Region",
               values_to = "Length")

df_PT_CT_PL$Region <- factor(x = df_PT_CT_PL$Region,
                             levels = c("MD", "RT", "PF", "SPA"))



library(ggpubr)
ggplot(data = df_PT_CT_PL, mapping = aes(x = AxonType, y = Length, fill = AxonType)) +
  facet_wrap(facets = "Region", nrow = 1, scales = "free_y") +
  geom_boxplot() +
  #geom_violin() +
  #geom_jitter(width = 0.1) +
  theme_classic() +stat_compare_means() +
  ggtitle("PT vs. CT. PL")


## read PL neurons ----

library(nat)
nl_neurite_PT_CT_PL <-
  nat::read.neurons(
    paths = SWCPath_neurite[index_axon_PT_CT[label_region_L2 == "PL"]],
    neuronnames = paste0(info_tmp$SampleID, "-", info_tmp$NeuronID)[index_axon_PT_CT[label_region_L2 == "PL"]])

hemisphere <- info_tmp$Hemisphere[index_axon_PT_CT[label_region_L2 == "PL"]]

## get index of PT ----
which(proj_PT_CT_TH$RT[label_region_L2 == "PL"] < 1000 &
        proj_PT_CT_TH$PF[label_region_L2 == "PL"] > 500 &
        proj_PT_CT_TH$SPA[label_region_L2 == "PL"] > 500 &
        info_tmp$Region[index_axon_PT_CT][label_region_L2 == "PL"] == "PL" &
        info_tmp$ProjectID[index_axon_PT_CT][label_region_L2 == "PL"] == "ISO" &
        info_tmp$AxonType[index_axon_PT_CT][label_region_L2 == "PL"] == "PT")

index_tmp_PL_PT <- c(125, 127)


## get index of CT ----
which(proj_PT_CT_TH$RT[label_region_L2 == "PL"] > 2000 &
        proj_PT_CT_TH$MD[label_region_L2 == "PL"] > 10000 &
        proj_PT_CT_TH$PF[label_region_L2 == "PL"] == 0 &
        proj_PT_CT_TH$SPA[label_region_L2 == "PL"] == 0 &
        info_tmp$Region[index_axon_PT_CT][label_region_L2 == "PL"] == "PL" &
        info_tmp$ProjectID[index_axon_PT_CT][label_region_L2 == "PL"] == "ISO" &
        info_tmp$AxonType[index_axon_PT_CT][label_region_L2 == "PL"] == "CT")

index_tmp_PL_CT <- c(246, 251)

open3d()
snp::plot_morphology(
  neuronList = nl_neurite_PT_CT_PL[index_tmp_PL_PT[1]],
  hemisphere = hemisphere[index_tmp_PL_PT[1]],
  color = "blue")

snp::plot_morphology(
  neuronList = nl_neurite_PT_CT_PL[index_tmp_PL_CT[1]],
  hemisphere = hemisphere[index_tmp_PL_CT[1]],
  color = "red")


## read mesh of RT, MD, PF, SPA ----

mesh_MD <- readOBJ("~/Documents/YanLab/database_neuron/ABA_mesh/362_MD.obj")
mesh_l_MD <- getMeshLeft(mesh_MD)

mesh_RT <- readOBJ("~/Documents/YanLab/database_neuron/ABA_mesh/262_RT.obj")
mesh_l_RT <- getMeshLeft(mesh_RT)

mesh_PF <- readOBJ("~/Documents/YanLab/database_neuron/ABA_mesh/930_PF.obj")
mesh_l_PF <- getMeshLeft(mesh_PF)

mesh_SPA <- readOBJ("~/Documents/YanLab/database_neuron/ABA_mesh/609_SPA.obj")
mesh_l_SPA <- getMeshLeft(mesh_SPA)


shade3d(mesh_l_MD, color = "green", alpha = 0.1)
shade3d(mesh_l_RT, color = "green", alpha = 0.1)

shade3d(mesh_l_PF, color = "green", alpha = 0.1)
shade3d(mesh_SPA, color = "green", alpha = 0.1)

shade3d(mesh_root, color = "white", alpha = 0.1)

rgl::view3d(userMatrix = userMatrix_horizontal, fov = 0)

