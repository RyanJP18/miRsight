args <- commandArgs(trailingOnly = TRUE)
suppressMessages(library(jsonlite))
config <- jsonlite::fromJSON(args[1])
settings <- config$settings
directories <- config$directories

suppressMessages(library(mice))
suppressMessages(library(parallel))


# --- FUNCTIONS ---

process_mirna <- function(i) {
    if (as.logical(settings$use_caching) && file.exists(file.path(directories$features_full_imputed, feature_files[i]))) {
        message("Imputation ", i, "/", n, " - loaded from cache.")
    } else {
        features <- read.table(file.path(directories$features_cons_shape, feature_files[i]), sep = "\t", header = TRUE, stringsAsFactors=TRUE)

        # for (j in seq_len(nrow(features))) {
        #     if (is.na(features$au_content_5[j])) {
        #         features$au_content_5_weighted[j] <- NA # TODO: currently there is a bug in the way weightings are computed when the source feature has an NA, rather than fix it here need to fix it in the extraction logic ideally
        #     }
        # }

        features$rel_utr_pos = features$binding_site_pos / features$X3_utr_length

        features$X3_utr_length = log10(features$X3_utr_length)
        features$cds_length = log10(features$cds_length)
        features$dist_closest_utr_end = log10(features$dist_closest_utr_end)
        features$binding_site_pos = log10(features$binding_site_pos)

        features_no_rn_l2fc <- features[c('gu_1', 'gu_8', 'perfect_pair_9', 'gu_9', 'perfect_pair_10', 'gu_10', 'perfect_pair_1', 'longest_any_sequence_12_17', 
            'longest_any_sequence_start_12_17', 'perfect_pair_count_09_20', 'any_pair_avg_dist_09_20', 'gu_count_09_20', 'mrna_binding_spread_full', 'site_abundance_6mer', 
            'site_abundance_6off', 'site_abundance_7mer', 'site_abundance_8mer', 'site_abundance_6cds', 'site_abundance_7cds', 'site_abundance_8cds', 'au_content_3', 
            'au_content_sup', 'au_content_5_weighted', 'rnafold_mfe', 'rnafold_direction', 'rnacofold_full_mfe', 'rnacofold_seed_mfe', 'rnaplfold_seed', 'rnaplfold_sup', 
            'X3_utr_length', 'cds_length', 'dist_closest_utr_end', 'best_abundance', 'phylo100_seed', 'phylo100_sup', 'phylo100_3', 'phylo100_5', 'shape_seed', 'shape_sup', 
            'rel_utr_pos', 'seed_binding_type', 'mirna_1', 'mirna_8', 'mrna_8')]

        # mice imputation
        m <- complete(mice(features_no_rn_l2fc, pred = quickpred(features_no_rn_l2fc), method = "cart", print = FALSE))
        features$shape_seed <- m$shape_seed
        features$shape_sup <- m$shape_sup
        features$phylo100_sup <- m$phylo100_sup
        features$phylo100_seed <- m$phylo100_seed
        features$rnaplfold_sup <- m$rnaplfold_sup
        features$rnaplfold_seed <- m$rnaplfold_seed
        features$au_content_sup <- m$au_content_sup
        features$au_content_3 <- m$au_content_3
        features$au_content_5_weighted <- m$au_content_5_weighted
        features$phylo100_5 <- m$phylo100_5
        features$phylo100_3 <- m$phylo100_3
        write.table(features, file.path(directories$features_full_imputed, feature_files[i]), quote = FALSE, sep = "\t", row.names = FALSE)

        message("Imputation ", i, "/", n, " - done.")
    }
}



# --- ENTRY POINT ---

feature_files <- dir(directories$features_cons_shape, pattern = ".tsv")

n <- length(feature_files)
max_cores <- as.numeric(settings$max_cores)
cores <- ifelse(max_cores == -1, detectCores() - 1, cores <- max_cores)

result <- mclapply(1:n, process_mirna, mc.cores = cores)