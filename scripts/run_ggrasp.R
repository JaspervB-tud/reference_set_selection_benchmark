library(ggrasp)
set.seed(1234)

args <- commandArgs(trailingOnly = TRUE)
distance_matrix_path <- args[1]
output_path <- args[2]

gg.1 <- ggrasp.load(file=distance_matrix_path, file.format="matrix", tree.method="single")
gg.2 <- ggrasp.cluster(gg.1)

if (!dir.exists(output_path)) {
    dir.create(output_path, recursive = TRUE)
}
output_file <- file.path(output_path, "ggrasp")
writeLines(gg.2@medoids, output_file)