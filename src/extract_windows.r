args <- commandArgs(trailingOnly = TRUE)
suppressMessages(library(jsonlite))
config <- jsonlite::fromJSON(args[1])
settings <- config$settings
directories <- config$directories

suppressMessages(library(dplyr))
suppressMessages(library(parallel))



# --- FUNCTIONS ---

create_new_frame <- dget("src/functions/create_new_frame.r")
reverse_complement <- dget("src/functions/reverse_complement.r")

# process and write out fold information to file
store_folding_windows <- function(folding_windows, experiment_name) {
  write.table(folding_windows$fold_window_lr, file.path(directories$windows_rnafold_lr, paste0(experiment_name, ".txt")), sep = "", col.names = FALSE, row.names = FALSE, quote = FALSE)
  write.table(folding_windows$fold_window_rl, file.path(directories$windows_rnafold_rl, paste0(experiment_name, ".txt")), sep = "", col.names = FALSE, row.names = FALSE, quote = FALSE)
  write.table(folding_windows$fold_window_ctr, file.path(directories$windows_rnafold_ctr, paste0(experiment_name, ".txt")), sep = "", col.names = FALSE, row.names = FALSE, quote = FALSE)

  # convert into a single column dataframe to better suit RNAcofold and then output to file
  cofold_full_interleaved <- c(rbind(folding_windows$cofold_window_full, folding_windows$cofold_constraint_full))
  write.table(cofold_full_interleaved, file.path(directories$windows_rnacofold_full, paste0(experiment_name, ".txt")), sep = "", col.names = FALSE, row.names = FALSE, quote = FALSE)

  # same again but with seed only
  cofold_seed_interleaved <- c(rbind(folding_windows$cofold_window_seed, folding_windows$cofold_constraint_seed))
  write.table(cofold_seed_interleaved, file.path(directories$windows_rnacofold_seed, paste0(experiment_name, ".txt")), sep = "", col.names = FALSE, row.names = FALSE, quote = FALSE)

  write.table(folding_windows$rnaplfold_window, file.path(directories$windows_rnaplfold, paste0(experiment_name, ".txt")), sep = "", col.names = FALSE, row.names = FALSE, quote = FALSE)
}

# find any 6mer or 6mer offset binding sites and count them up to determine their abundance values, also check the cds
extract_folding_windows <- function(expanded_binding_sites, utrs, mirna_sequence) {
    folding_windows <- create_new_frame(
        c("ensembl_transcript_id_version",
        "fold_window_lr", "fold_window_rl", "fold_window_ctr",
        "cofold_window_full", "cofold_constraint_full",
        "cofold_window_seed", "cofold_constraint_seed",
        "rnaplfold_window", "rnaplfold_6mer_pos"),
        NULL, nrow(utrs))
    folding_windows$ensembl_transcript_id_version <- utrs$ensembl_transcript_id_version

    # find the current miRNA sequence, search for its 6mer site and RC it to get the 6mer target site
    site_6mer <- substr(mirna_sequence, start = 2, stop = 7) # miRNA is fetched from 5` -> 3`
    target_6mer <- reverse_complement(site_6mer)

    for (i in seq_len(nrow(expanded_binding_sites))) {
        utr <- utrs$X3utr[i]
        binding_site_pos <- expanded_binding_sites$binding_site_pos[i]

        # rnafold windows
        fold_window_lr <- substr(utr, start = binding_site_pos - folding_window_size - 1, stop = binding_site_pos + 6)
        folding_windows$fold_window_lr[i] <- fold_window_lr # Get a window from left to right, ending with the seed portion
        folding_windows$fold_window_rl[i] <-
            substr(utr, start = binding_site_pos - 1, stop = binding_site_pos - 1 + 6 + folding_window_size + 1) # Get a window from right to left, beginning with the seed portion
        folding_windows$fold_window_ctr[i] <-
            substr(utr, start = binding_site_pos - 1 - (folding_window_size * 0.5), stop = binding_site_pos + 6 + (folding_window_size * 0.5)) # Half window, seed, half window

        # rnacofold full window
        full_mirna_constraint <- gsub(site_6mer, "||||||", mirna_sequence)
        full_mirna_constraint <- gsub("[[:alpha:]]", ".", full_mirna_constraint)

        full_mrna_constraint <- substr(fold_window_lr, start = 8, stop = nchar(fold_window_lr)) # multiple seed targets could exist in the window and we only want the last
        full_mrna_constraint <- paste0(gsub("[[:alpha:]]", ".", full_mrna_constraint), "||||||.")

        folding_windows$cofold_window_full[i] <- paste0(mirna_sequence, "&", fold_window_lr)
        folding_windows$cofold_constraint_full[i] <- paste0(full_mirna_constraint, "&", full_mrna_constraint)

        # rnacofold seed windows
        mirna_seed_only <- substr(mirna_sequence, start = 1, stop = 8)
        mrna_seed_only <- substr(utr, start = binding_site_pos - 1, stop = binding_site_pos + 6) # Get just the seed portion of the window

        seed_only_mrna_constraint <- gsub(target_6mer, "||||||", mrna_seed_only)
        seed_only_mrna_constraint <- gsub("[[:alpha:]]", ".", seed_only_mrna_constraint)

        seed_only_mirna_constraint <- gsub(site_6mer, "||||||", mirna_seed_only)
        seed_only_mirna_constraint <- gsub("[[:alpha:]]", ".", seed_only_mirna_constraint)

        folding_windows$cofold_window_seed[i] <- paste0(mirna_seed_only, "&", mrna_seed_only)
        folding_windows$cofold_constraint_seed[i] <- paste0(seed_only_mirna_constraint, "&", seed_only_mrna_constraint)

        # rnaplfold window
        substr_start <- binding_site_pos - 1 - half_rnapl_window_size
        substr_stop <- binding_site_pos + 6 + half_rnapl_window_size

        rnaplfold_start_6mer <- 1
        if (substr_start >= 1) {
            rnaplfold_start_6mer <- half_rnapl_window_size + 2
        } else if (substr_start == 0) {
            rnaplfold_start_6mer <- half_rnapl_window_size + 1
        } else {
            rnaplfold_start_6mer <- half_rnapl_window_size + substr_start + 1
        }

        folding_windows$rnaplfold_window[i] <- substr(utr, start = substr_start, stop = substr_stop)
        folding_windows$rnaplfold_6mer_pos[i] <- rnaplfold_start_6mer
    }

    return(folding_windows)
}

process_mirna <- function(i) {

    if (as.logical(settings$use_caching) && file.exists(file.path(directories$windows, binding_site_files[i]))) {
        message(paste0("Window extraction ", i, "/", n, " - loaded from cache"))
    } else {
        start_time <- Sys.time()

        binding_site_filename <- binding_site_files[i]
        mirna_id <- mirna_ids[i]
        window_path <- file.path(directories$windows, binding_site_files[i])
        
        raw_binding_sites <- read.table(file.path(directories$bindings_raw, binding_site_filename), sep = "\t", header = TRUE)
        expanded_binding_sites <- read.table(file.path(directories$bindings, binding_site_filename), sep = "\t", header = TRUE)

        utrs <- utrs[rep(1:nrow(raw_binding_sites), raw_binding_sites$site_abundance_6mer), ]

        mirna_sequence <- mirna_sequences[mirna_sequences$mirna_id == mirna_id, ]$mirna_sequence

        folding_windows <- extract_folding_windows(expanded_binding_sites, utrs, mirna_sequence)

        write.table(folding_windows, window_path, quote = FALSE, sep = "\t", row.names = FALSE)

        store_folding_windows(folding_windows, gsub("*.tsv", "", binding_site_filename))

        message(paste0("Window extraction ", i, "/", n, " - done.", "-- Elapsed-- ", Sys.time() - start_time))
    }
}


# --- ENTRY POINT ---

mirna_sequences <- read.table(file.path(directories$annotations, "mirna_sequences.tsv"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
if (settings$mirna_id_filter != "") {
    mirna_sequences <- mirna_sequences[mirna_sequences$mirna_id %in% strsplit(settings$mirna_id_filter, ",")[[1]], ]
}

annotations <- read.table(file.path(directories$annotations, "annotations.tsv"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
mane <- read.table(file.path(directories$annotations, "mane.tsv"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
utrs <- read.table(file.path(directories$annotations, "utr_sequences.tsv"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# apply filters
if (settings$chromosome_filter != "") {
    annotations <- annotations[annotations$chromosome_name %in% strsplit(settings$chromosome_filter, ",")[[1]], ]
}
# if (settings$ensembl_transcript_id_filter != "") {
#     annotations <- annotations[sub("\\..*", "", annotations$ensembl_transcript_id_version) %in% strsplit(settings$ensembl_transcript_id_filter, ",")[[1]], ]
# }
# if (settings$ensembl_gene_id_filter != "") {
#     annotations <- annotations[annotations$ensembl_gene_id %in% strsplit(settings$ensembl_gene_id_filter, ",")[[1]], ]
# }
# if (settings$external_gene_id_filter != "") {
#     annotations <- annotations[annotations$external_gene_id %in% strsplit(settings$external_gene_id_filter, ",")[[1]], ]
# }

# filter to only transcripts where we have utr sequence annotations
annotations <- annotations[annotations$ensembl_transcript_id_version %in% mane$ensembl_transcript_id_version, ]
utrs <- utrs[utrs$ensembl_transcript_id_version %in% annotations$ensembl_transcript_id_version, ]
utrs <- utrs[utrs$X3utr != "Sequence unavailable", ]
utrs_all <- utrs

folding_window_size <- as.numeric(settings$folding_window_size)
rnaplfold_window_size <- as.numeric(settings$rnaplfold_window_size)
half_rnapl_window_size <- rnaplfold_window_size * 0.5

binding_site_files <- dir(directories$bindings, pattern = ".tsv")
binding_site_files <- head(binding_site_files, length(binding_site_files) - 1)
mirna_ids <- gsub("*.tsv", "", binding_site_files)

n <- length(binding_site_files)
max_cores <- as.numeric(settings$max_cores)
cores <- ifelse(max_cores == -1, detectCores() - 1, cores <- max_cores)

result <- mclapply(1:n, process_mirna, mc.cores = cores)