# Run from the repository root.
#Creates a summary csv for each variant in an input genotype file
#Outputs into results folder

source("R/variant_functions.R")

#Settings
input_path  <- "data/genotypes.csv"
out_path  <- "results"
min_maf     <- 0.25
max_missing <- 0.25

genotypes <- read.csv(input_path, colClasses = "numeric")
out <- data.frame()

#Summary generation
for (v in names(genotypes)) {
    summ <- summarize_variant(genotypes[[v]],
                              min_maf = min_maf,
                              max_missing = max_missing)
    out <- rbind(out, cbind(variant_id = v, summ))
}

if (!dir.exists(out_path)) {
    dir.create(out_path)
}

write.csv(out, file = file.path(out_path, "summary.csv"), quote=FALSE, row.names = FALSE)
