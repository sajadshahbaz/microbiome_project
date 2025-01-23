library(MDFS)

# Initialize lists to store results
MDFS_group <- list()
rel_group <- list()
MDFS_1D_rel <- list()
  MDFS_2D_rel <- list()
  MDFS_1D2D_rel <- list()
  MDFS_1D <- list()
  MDFS_2D <- list()
  MDFS_seed_no<-vector()
for (i in 1:30) {
  # Initialize lists for each iteration
  current_seed <- i  # Use 'i' as the seed value
  MDFS_seed_no[length(MDFS_seed_no)+1]<-current_seed
  set.seed(current_seed)
  # Perform MDFS with 1D
  MDFS_1D_result <- MDFS(
    data = randomdataset[[i]][["binary_taxa_sampled"]],
    decision = randomdataset[[i]][["y_sampled"]],
    dimensions = 1,
    divisions = 1,
    discretizations = 1,
    range = 0,
    seed= current_seed,  # Set the seed for reproducibility
    level = 0.05
  )

  MDFS_1D[[length(MDFS_1D) + 1]] <- MDFS_1D_result
  MDFS_1D_rel[[length(MDFS_1D_rel) + 1]] <- colnames(randomdataset[[i]][["binary_taxa_sampled"]])[MDFS_1D_result$relevant.variables]

  # Perform MDFS with 2D
  MDFS_2D_result <- MDFS(
    data = randomdataset[[i]][["binary_taxa_sampled"]],
    decision = randomdataset[[i]][["y_sampled"]],
    dimensions = 2,
    divisions = 1,
    discretizations = 1,
    range = 0,
    seed = current_seed,
    level = 0.05
  )

  MDFS_2D[[length(MDFS_2D) + 1]] <- MDFS_2D_result
  MDFS_2D_rel[[length(MDFS_2D_rel) + 1]] <- colnames(randomdataset[[i]][["binary_taxa_sampled"]])[MDFS_2D_result$relevant.variables]

  # Combine 1D and 2D relevant variables
  MDFS_1D2D_rel[[length(MDFS_1D2D_rel) + 1]] <- union(MDFS_1D_rel[[length(MDFS_1D_rel)]], MDFS_2D_rel[[length(MDFS_2D_rel)]])
}

# Combine results into groups
MDFS_group <- list(MDFS_1D = MDFS_1D, MDFS_2D = MDFS_2D)
rel_group <- list(MDFS_1D_rel = MDFS_1D_rel, MDFS_2D_rel = MDFS_2D_rel, MDFS_1D2D_rel = MDFS_1D2D_rel)

# Save results
save(rel_group, MDFS_group, file = "rel_MDFS_group.RData")
#
