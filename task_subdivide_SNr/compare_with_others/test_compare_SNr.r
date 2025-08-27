# project mouse whole cortex
# compare with SNr domains identified in Foster et al. Nature, 2021
# test. use the raw data deposited in the website

SNr_level_81 <-
  read.csv(file = "/home/admin/Downloads/SNr81_boxgrid_data.csv",
           row.names = 1)
SNr_level_81[is.na(SNr_level_81)] <- 0

tmp_split <-
sapply(X = colnames(SNr_level_81), FUN = function(x) {
  tmp_split <- strsplit(x = x, split = "r.")[[1]][2]
  strsplit(x = tmp_split, split = "[.]")[[1]]
})

coord_X <- as.numeric(tmp_split[1,])
coord_Y <- as.numeric(tmp_split[2,])


summary(coord_X)
summary(coord_Y)

intensity_sum <- colSums(SNr_level_81)

matrix_intensity <- matrix(data = 0, nrow = max(coord_X), ncol = max(coord_Y))

index_grid <-
sapply(X = 1:length(coord_X), FUN = function(x) {
  nat::sub2ind(dims = dim(matrix_intensity), indices = c(coord_X[x], coord_Y[x]))
})



matrix_intensity[index_grid] <- intensity_sum

Heatmap(matrix = t(matrix_intensity[75:90,55:69]),
        cluster_rows = F, cluster_columns = F, height = 5, width = 5)

