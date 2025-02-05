library(ggrasp)
set.seed(1234)

args <- commandArgs(trailingOnly = TRUE)
species <- args[1]

distance_matrix_path <- paste0("./genomes_per_species/", species, "/converted_matrix.mat")
output_path <- paste0("./selections/", species, "/ggrasp")

gg.1 <- ggrasp.load(file=distance_matrix_path, file.format="matrix", tree.method="single")
gg.2 <- ggrasp.cluster(gg.1)

writeLines(gg.2@medoids, output_path)