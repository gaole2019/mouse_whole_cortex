# project mouse whole cortex
# cortical columns
# volume of CC255 and ABA43

statistics_CC <- read.csv(
  file = "/run/media/admin/LeGAO/CodeMatlab/Isocortex_ROI/result/statistics_volumes.csv",
  header = T)

mean(statistics_CC$Number.Of.Voxels[2:256] * 25*25*25 / 1000/1000/1000)

sd(statistics_CC$Number.Of.Voxels[2:256] * 25*25*25 / 1000/1000/1000) / sqrt(255)

summary(statistics_CC$Number.Of.Voxels[2:256] * 25*25*25 / 1000/1000/1000)

# 1. obtain the basic regions of each of 43 cortical region, since anno_ABA.nrrd only contains the basic regions.
region_StructureID <-
  substring(snp::getStructureFromAcronym(paste0("l_", shared_Isocortex_subregion)), 3)

label_group <- c()
subregion_StructureID <- c()
subregion_Acronym <- c()

NumRegion <- length(region_StructureID)
iRegion <- 1
for (iRegion in 1:NumRegion) {
  tmp_subregion_StructureID <-
    snp::shared_allen_anno$CurrentID[snp::shared_allen_anno$ParentID %in% region_StructureID[iRegion]]
  tmp_subregion_Acronym <-
    snp::shared_allen_anno$Acronym[snp::shared_allen_anno$ParentID %in% region_StructureID[iRegion]]

  label_group <- c(label_group, rep(iRegion, length(tmp_subregion_StructureID)))
  subregion_StructureID <- c(subregion_StructureID, tmp_subregion_StructureID)
  subregion_Acronym <- c(subregion_Acronym, tmp_subregion_Acronym)
}

df_RegionID_StructureID <- data.frame(RegionID = label_group,
                                      StuctureID = subregion_StructureID,
                                      Acronym = subregion_Acronym)

write.csv(
  x = df_RegionID_StructureID,
  file = "/run/media/admin/LeGAO2/project_mouse_cortex/analysis/task_volume/ABA_43.csv",
  row.names = F)





