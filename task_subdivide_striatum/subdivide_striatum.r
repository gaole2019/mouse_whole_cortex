# project mouse whole cortex
# figure 4. cortico-basal ganglia projections
# subdivide striatum



# 1. fig4K. example IT and PT neuron projecting to striatum ----

index_example_IT <-
  which(shared_proj_Combined[,snp::getStructureFromAcronym("l_CP")] > 30000 &
        1:21582 %in% df_fnt_dist$Index[df_fnt_dist$AxonType %in% c("IT")] &
          info_tmp$Region == "SSp-bfd")

index_example_PT <-
  which(shared_proj_Combined[,snp::getStructureFromAcronym("l_CP")] > 30000 &
          1:21582 %in% df_fnt_dist$Index[df_fnt_dist$AxonType %in% c("PT")] &
          info_tmp$Region == "SSp-bfd")


nl_neuron_IT <- nat::read.neuron(f = SWCPath_neurite[index_example_IT[2]])
if (nl_neuron_IT$d$Z[1] > 5700) {
  nl_neuron_IT$d$Z <- 11400 - nl_neuron_IT$d$Z
}

nl_neuron_PT <- nat::read.neuron(f = SWCPath_neurite[index_example_PT[2]])
if (nl_neuron_PT$d$Z[1] > 5700) {
  nl_neuron_PT$d$Z <- 11400 - nl_neuron_PT$d$Z
}

mesh_root <- rgl::readOBJ("~/Documents/YanLab/database_neuron/ABA_mesh/997_root.obj")
mesh_CP <- rgl::readOBJ("~/Documents/YanLab/database_neuron/ABA_mesh/672_CP.obj")
mesh_ACB <- rgl::readOBJ("~/Documents/YanLab/database_neuron/ABA_mesh/56_ACB.obj")

mesh_l_CP <- getMeshLeft(mesh = mesh_CP)
mesh_l_ACB <- getMeshLeft(mesh = mesh_ACB)



open3d()
plot3d(x = nl_neuron_IT, color = "red", WithNodes = F, soma = 200)
# shade3d(mesh_CP, color = "red", alpha = 0.1)
# shade3d(mesh_root, color = "white", alpha = 0.1)

#open3d()
plot3d(x = nl_neuron_PT, color = "blue", WithNodes = F, soma = 200)


shade3d(mesh_l_CP, color = "white", alpha = 0.1)
shade3d(mesh_l_ACB, color = "white", alpha = 0.1)

shade3d(mesh_root, color = "white", alpha = 0.1)


# axon segments in striatum ----
IT_nl_neurite_in_CP <-
  nat::prune_in_volume(
    x = nl_neuron_IT,
    surf = nat::as.hxsurf(mesh_l_CP),
    invert = T, OmitFailures = T)

PT_nl_neurite_in_CP <-
  nat::prune_in_volume(
    x = nl_neuron_PT,
    surf = nat::as.hxsurf(mesh_l_CP),
    invert = T, OmitFailures = T)


open3d()
plot3d(x = IT_nl_neurite_in_CP, color = "red", WithNodes = F)
plot3d(x = PT_nl_neurite_in_CP, color = "blue", WithNodes = F)

shade3d(mesh_l_CP, color = "white", alpha = 0.1)
shade3d(mesh_l_ACB, color = "white", alpha = 0.1)




# 2. fig 4M. IT and PT neurons projecting to striatum ----
index_striatum_IT <-
  which(rowSums(shared_proj_Combined[,snp::getStructureFromAcronym(c("l_CP", "l_ACB"))]) > 1000 &
        1:21582 %in% df_fnt_dist$Index[df_fnt_dist$AxonType %in% c("IT")])

index_striatum_PT <-
  which(rowSums(shared_proj_Combined[,snp::getStructureFromAcronym(c("l_CP", "l_ACB"))]) > 1000 &
        1:21582 %in% df_fnt_dist$Index[df_fnt_dist$AxonType %in% c("PT")])



table(info_tmp$Hemisphere[info_tmp$ProjectID == "PFC"])


df_subdivide_striatum <-
data.frame(
  Hemisphere = info_tmp$Hemisphere[c(index_striatum_IT, index_striatum_PT)],
  ProjectID = info_tmp$ProjectID[c(index_striatum_IT, index_striatum_PT)],
  AxonType = c(rep("IT", length(index_striatum_IT)),
               rep("PT", length(index_striatum_PT))),
  SWCPath = shared_SWCPath_axon_allen[c(index_striatum_IT, index_striatum_PT)])

df_subdivide_striatum$Hemisphere[df_subdivide_striatum$ProjectID == "PFC"] <- "Left"


write.csv(x = df_subdivide_striatum,
          file = "/run/media/admin/LeGAO2/project_mouse_cortex/analysis/task_subdivide_striatum/df_subdivide_striatum.csv")



axon_length_subdomain <-
  read.csv(
    file = "/run/media/admin/LeGAO2/project_mouse_cortex/analysis/task_subdivide_striatum/axon_length_in_strital_subdomains.csv",
    row.names = 1)


identical(df_subdivide_striatum$SWCPath, rownames(axon_length_subdomain))


# 4. read annotation file of striatal subdomains ----
anno_subdomain <-
  nat::read.nrrd(
    file = "/run/media/admin/LeGAO/Project_Cortex/task_subdivide_striatum/label_cluster_real_K9_800_add_nonPFC.nrrd",
    ReadByteAsRaw = F)

# 5. calculate volumes of each subdomain ----
volume_subdomain <- sapply(1:14, FUN = function(x) length(which(anno_subdomain == x)))

label_major_subdomain <-
  apply(X = axon_length_subdomain, MARGIN = 1, FUN = function(x) which.max(x/volume_subdomain))

label_major_subdomain_ratio <-
  apply(X = axon_length_subdomain, MARGIN = 1,
        FUN = function(x) {
          proj_density <- x / volume_subdomain
          index_max <- which.max(proj_density)
          proj_density[index_max] / mean(proj_density)
          })



# 6. flatmap of neurons preferentially projecting to each striatal subdomain ----
ColorSubdomain <- read.csv("~/shared/colors_17.txt")
ColorSubdomain <- rgb(red = ColorSubdomain$R,
                      green = ColorSubdomain$G,
                      blue = ColorSubdomain$B,
                      maxColorValue = 255)
names(ColorSubdomain) <- 1:17



#
U <- shared_soma_XYZUV$U[c(index_striatum_IT, index_striatum_PT)][df_subdivide_striatum$AxonType == "PT"]
V <- shared_soma_XYZUV$V[c(index_striatum_IT, index_striatum_PT)][df_subdivide_striatum$AxonType == "PT"]
Color <- ColorSubdomain[match(label_major_subdomain[df_subdivide_striatum$AxonType == "PT"],
                             names(ColorSubdomain))]
Radius = 6
filename = "/tmp/PT_strital_subdomains.svg"

conn <- file("/home/admin/shared/image_flatmap_color_v3.svg")
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
                           U[iElement], V[iElement], Radius[iElement], Color[iElement]))
}

out <- c(svg_PFC[1:18], svg_element, svg_PFC[20:length(svg_PFC)])
writeLines(out, fileConn, sep = "\n")
close(fileConn)


# 7. figure 4N. example neurons in ORBm, AIv, ECT projecting to strital subdomain 11 ----

## ORBm ----
index_example_IT_subdomain_11_ORBm <-
  intersect(which(info_tmp$Region == "ORBm" &
        1:21582 %in% df_fnt_dist$Index[df_fnt_dist$AxonType %in% c("IT")]),
        c(index_striatum_IT, index_striatum_PT)[label_major_subdomain == 11])


table(info_tmp$Region[index_example_IT_subdomain_11_ORBm])


nl_neuron_example_IT_subdomain_11_ORBm <-
  nat::read.neuron(f = shared_SWCPath_axon_allen[index_example_IT_subdomain_11_ORBm[1]])



## AIv ----
index_example_IT_subdomain_11_AIv <-
  intersect(which(info_tmp$Region == "AIv" &
                    1:21582 %in% df_fnt_dist$Index[df_fnt_dist$AxonType %in% c("IT")]),
            c(index_striatum_IT, index_striatum_PT)[label_major_subdomain == 11])
table(info_tmp$Region[index_example_IT_subdomain_11_AIv])


nl_neuron_example_IT_subdomain_11_AIv <-
  nat::read.neuron(f = shared_SWCPath_axon_allen[index_example_IT_subdomain_11_AIv[1]])

nl_neuron_example_IT_subdomain_11_AIv$d$Z <-
  11400 - nl_neuron_example_IT_subdomain_11_AIv$d$Z



## ECT ----
index_example_IT_subdomain_11_ECT <-
  intersect(which(info_tmp$Region == "ECT" &
                    1:21582 %in% df_fnt_dist$Index[df_fnt_dist$AxonType %in% c("IT")]),
            c(index_striatum_IT, index_striatum_PT)[which(label_major_subdomain == 11 &
                                                          label_major_subdomain_ratio > 9)])
table(info_tmp$Region[index_example_IT_subdomain_11_ECT])


nl_neuron_example_IT_subdomain_11_ECT <-
  nat::read.neuron(f = shared_SWCPath_axon_allen[index_example_IT_subdomain_11_ECT[1]])

nl_neuron_example_IT_subdomain_11_ECT$d$Z <-
  11400 - nl_neuron_example_IT_subdomain_11_ECT$d$Z


mesh_subdomain_11 <-
  rgl::readOBJ(con = "/run/media/admin/LeGAO/Project_Cortex/task_subdivide_striatum/rendering/structures_mesh/11.tiff_smoothed.obj")


mesh_CP <- rgl::readOBJ("~/Documents/YanLab/database_neuron/ABA_mesh/672_CP.obj")
mesh_ACB <- rgl::readOBJ("~/Documents/YanLab/database_neuron/ABA_mesh/56_ACB.obj")

mesh_l_CP <- getMeshLeft(mesh = mesh_CP)
mesh_l_ACB <- getMeshLeft(mesh = mesh_ACB)

library(rgl)
open3d()
plot3d(nl_neuron_example_IT_subdomain_11_ORBm, color = "red", WithNodes = F, soma = 200)
plot3d(nl_neuron_example_IT_subdomain_11_AIv, color = "green", WithNodes = F, soma = 200)
plot3d(nl_neuron_example_IT_subdomain_11_ECT, color = "blue", WithNodes = F, soma = 200)

shade3d(mesh_subdomain_11, color = "orange", alpha = 0.1)
shade3d(mesh_root, color = "white", alpha = 0.1)


mesh_subdomain_11 <-  Morpho::invertFaces(mesh = mesh_subdomain_11)
Rvcg::checkFaceOrientation(mesh_subdomain_11)


## axons in subdomain 11. zoomin ----
IT_nl_neurite_in_subdomain_11_ORBm <-
  nat::prune_in_volume(
    x = nl_neuron_example_IT_subdomain_11_ORBm,
    surf = nat::as.hxsurf(mesh_subdomain_11),
    invert = T, OmitFailures = T)

IT_nl_neurite_in_subdomain_11_AIv <-
  nat::prune_in_volume(
    x = nl_neuron_example_IT_subdomain_11_AIv,
    surf = nat::as.hxsurf(mesh_subdomain_11),
    invert = T, OmitFailures = T)

IT_nl_neurite_in_subdomain_11_ECT <-
  nat::prune_in_volume(
    x = nl_neuron_example_IT_subdomain_11_ECT,
    surf = nat::as.hxsurf(mesh_subdomain_11),
    invert = T, OmitFailures = T)


library(rgl)
open3d()
plot3d(IT_nl_neurite_in_subdomain_11_ORBm, color = "red", WithNodes = F, soma = 0)
plot3d(IT_nl_neurite_in_subdomain_11_AIv, color = "green", WithNodes = F, soma = 0)
plot3d(IT_nl_neurite_in_subdomain_11_ECT, color = "blue", WithNodes = F, soma = 00)

shade3d(mesh_subdomain_11, color = "white", alpha = 0.1)




# 8. figure 4O. example neurons in ACAv, VISp projecting to strital subdomain 11 ----

## ACA ----
index_example_IT_subdomain_14_ACAd <-
  intersect(which(info_tmp$Region == "ACAd" &
                    1:21582 %in% df_fnt_dist$Index[df_fnt_dist$AxonType %in% c("IT")]),
            c(index_striatum_IT, index_striatum_PT)[which(
              label_major_subdomain == 14 &
              label_major_subdomain_ratio > 9 &
              rowSums(axon_length_subdomain) > 20000)])

table(info_tmp$Region[index_example_IT_subdomain_14_ACAd])


nl_neuron_example_IT_subdomain_14_ACAd <-
  nat::read.neuron(f = shared_SWCPath_axon_allen[index_example_IT_subdomain_14_ACAd[1]])



## VISp ----
index_example_IT_subdomain_14_VISp <-
  intersect(which(info_tmp$Region == "VISp" &
                    1:21582 %in% df_fnt_dist$Index[df_fnt_dist$AxonType %in% c("IT")]),
            c(index_striatum_IT, index_striatum_PT)[which(
              label_major_subdomain == 14 &
              label_major_subdomain_ratio > 9 &
              rowSums(axon_length_subdomain) > 20000)])

table(info_tmp$Region[index_example_IT_subdomain_14_VISp])


nl_neuron_example_IT_subdomain_14_VISp <-
  nat::read.neuron(f = shared_SWCPath_axon_allen[index_example_IT_subdomain_14_VISp[1]])



mesh_subdomain_14 <-
  rgl::readOBJ(con = "/run/media/admin/LeGAO/Project_Cortex/task_subdivide_striatum/rendering/structures_mesh/14.tiff_smoothed.obj")


library(rgl)
open3d()
plot3d(nl_neuron_example_IT_subdomain_14_ACAd, color = "red", WithNodes = F, soma = 200)
plot3d(nl_neuron_example_IT_subdomain_14_VISp, color = "green", WithNodes = F, soma = 200)

shade3d(mesh_subdomain_14, color = "orange", alpha = 0.1)
shade3d(mesh_root, color = "white", alpha = 0.1)




## axons in subdomain 14. zoomin ----


mesh_subdomain_14 <-  Morpho::invertFaces(mesh = mesh_subdomain_14)
Rvcg::checkFaceOrientation(mesh_subdomain_14)


IT_nl_neurite_in_subdomain_14_ACAd <-
  nat::prune_in_volume(
    x = nl_neuron_example_IT_subdomain_14_ACAd,
    surf = nat::as.hxsurf(mesh_subdomain_14),
    invert = T, OmitFailures = T)

IT_nl_neurite_in_subdomain_14_VISp <-
  nat::prune_in_volume(
    x = nl_neuron_example_IT_subdomain_14_VISp,
    surf = nat::as.hxsurf(mesh_subdomain_14),
    invert = T, OmitFailures = T)


library(rgl)
open3d()
plot3d(IT_nl_neurite_in_subdomain_14_ACAd, color = "red", WithNodes = F, soma = 0)
plot3d(IT_nl_neurite_in_subdomain_14_VISp, color = "green", WithNodes = F, soma = 0)

shade3d(mesh_subdomain_14, color = "white", alpha = 0.1)

