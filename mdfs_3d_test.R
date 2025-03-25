
load("mean_data.RData")

library(MDFS)
set.seed(23423)
  MDFS_3D_result <- MDFS(
    data =  meta.mean.final,
    decision = decision.mean.final,
    dimensions = 3,
    divisions = 1,
    discretizations = 1,
    range = 0,
    seed = current_seed,
    level = 0.05
  )
saveRDS(MDFS_3D_result.rds,"MDFS_3D_result.rds")
