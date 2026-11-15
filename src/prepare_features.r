args <- commandArgs(trailingOnly = TRUE)
suppressMessages(library(jsonlite))
config <- jsonlite::fromJSON(args[1])
settings <- config$settings
directories <- config$directories
cores <- as.numeric(args[2])


# --- FUNCTIONS ---

process_mirna <- function(i) {
    if (as.logical(settings$use_caching) && file.exists(file.path(directories$features_prepared, feature_files[i]))) {
        message("Feature preparation ", i, "/", n, " - loaded from cache.")
    } else {
        features <- read.table(file.path(directories$features_cons_shape, feature_files[i]), sep = "\t", header = TRUE, stringsAsFactors = TRUE)
        features$rel_utr_pos = features$binding_site_pos / features$X3_utr_length
        features$X3_utr_length = log10(features$X3_utr_length)
        features$cds_length = log10(features$cds_length)
        features$dist_closest_utr_end = log10(features$dist_closest_utr_end)
        features$binding_site_pos = log10(features$binding_site_pos)
        features$gu_1 <- NULL
        features$gu_8 <- NULL
        features$gu_9 <- NULL
        features$gu_10 <- NULL
        features$rnafold_direction <- NULL
        features$perfect_pair_9 <- NULL
        features$perfect_pair_10 <- NULL
        features$perfect_pair_1 <- NULL
        features$mrna_8 <- NULL
        write.table(features, file.path(directories$features_prepared, feature_files[i]), quote = FALSE, sep = "\t", row.names = FALSE)

        message("Feature preparation ", i, "/", n, " - done.")
    }
}



# --- ENTRY POINT ---

feature_files <- dir(directories$features_cons_shape, pattern = ".tsv")
if (settings$mirna_id_filter != "") {
    feature_files <- feature_files[feature_files %in% paste0(strsplit(settings$mirna_id_filter, ",")[[1]], ".tsv")]
}
n <- length(feature_files)

result <- mclapply(1:n, process_mirna, mc.cores = cores)
