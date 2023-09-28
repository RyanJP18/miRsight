args <- commandArgs(trailingOnly = TRUE)
suppressMessages(library(jsonlite))
config <- jsonlite::fromJSON(args[1])
settings <- config$settings
directories <- config$directories



# --- FUNCTIONS ---

create_new_frame <- dget("src/functions/create_new_frame.r")
reverse_complement <- dget("src/functions/reverse_complement.r")

# compare all rnafolds to find the one with the strongst binding, return the data of that one only
determine_most_likely_folds <- function(rnafolds, rnafolds_lr, rnafolds_rl, rnafolds_ctr) {
    for (i in seq_len(nrow(rnafolds))) {
        mfe_lr <- as.numeric(as.character(paste0("-", gsub("^\\D+", "", gsub("[^0-9.-]", "", rnafolds_lr$rnafold_struct[i])))))
        mfe_rl <- as.numeric(as.character(paste0("-", gsub("^\\D+", "", gsub("[^0-9.-]", "", rnafolds_rl$rnafold_struct[i])))))
        mfe_ctr <- as.numeric(as.character(paste0("-", gsub("^\\D+", "", gsub("[^0-9.-]", "", rnafolds_ctr$rnafold_struct[i])))))
        if (mfe_lr <= mfe_rl) {
            if (mfe_lr <= mfe_ctr) {
                rnafolds[i, 2:4] <- c(substr(as.character(rnafolds_lr$rnafold_struct[i]), 1, nchar(as.character(rnafolds_lr$rnafold_struct[i])) - 9), mfe_lr, -1)
            } else {
                rnafolds[i, 2:4] <- c(substr(as.character(rnafolds_ctr$rnafold_struct[i]), 1, nchar(as.character(rnafolds_ctr$rnafold_struct[i])) - 9), mfe_ctr, 0)
            }
        } else {
            if (mfe_rl <= mfe_ctr) {
                rnafolds[i, 2:4] <- c(substr(as.character(rnafolds_rl$rnafold_struct[i]), 1, nchar(as.character(rnafolds_rl$rnafold_struct[i])) - 9), mfe_rl, 1)
            } else {
                rnafolds[i, 2:4] <- c(substr(as.character(rnafolds_ctr$rnafold_struct[i]), 1, nchar(as.character(rnafolds_ctr$rnafold_struct[i])) - 9), mfe_ctr, 0)
            }
        }
    }

    return(rnafolds)
}

# get the rnafold data for the mirna and put it into a single dataframe
load_rnafolds <- function(experiment_name, expanded_binding_sites, load_rnafolds_lr, load_rnafolds_rl) {
    rnafolds <- create_new_frame(c("ensembl_transcript_id_version", "rnafold_struct", "rnafold_mfe", "rnafold_direction"), NULL, nrow(expanded_binding_sites))
    rnafolds$ensembl_transcript_id_version <- expanded_binding_sites$ensembl_transcript_id_version

    rnafolds_ctr_raw <- read.table(file.path(directories$folds_rnafold_ctr, paste0(experiment_name, ".csv")), sep = ",", header = FALSE)
    rnafolds_ctr <- do.call("cbind", split(rnafolds_ctr_raw, rep(c(1, 2), length.out = nrow(rnafolds_ctr_raw))))
    colnames(rnafolds_ctr) <- c("mirna_sequence", "rnafold_struct")

    return(determine_most_likely_folds(rnafolds, rnafolds_lr, rnafolds_rl, rnafolds_ctr))
}

# get the rnacofold data for the mirna and put it into a single dataframe
load_rnacofolds <- function(experiment_name, expanded_binding_sites) {
    rnacofolds <- create_new_frame(c("ensembl_transcript_id_version"), NULL, nrow(expanded_binding_sites))
    rnacofolds$ensembl_transcript_id_version <- expanded_binding_sites$ensembl_transcript_id_version

    # get the rna cofold data and bind it to the "targets" data frame so we have everything in one place (transcript IDs, rnafold changs and cofolds)
    rnacofolds_full <- read.table(file.path(directories$folds_rnacofold_full, paste0(experiment_name, ".csv")), sep = ",", header = TRUE)
    rnacofolds_full <- rnacofolds_full[, !(names(rnacofolds_full) %in% c("seq_num", "seq_id"))] # drop unneeded columns
    colnames(rnacofolds_full) <- c("rnacofold_full_sequence", "rnacofold_full_struct", "rnacofold_full_mfe")
    rnacofolds <- cbind(rnacofolds, rnacofolds_full) # bind the targets and cofolds together

    # get the seed cofold data and bind it to the "targets" data frame so we have everything in one place (transcript IDs, rnafold changs and cofolds)
    rnacofolds_seed <- read.table(file.path(directories$folds_rnacofold_seed, paste0(experiment_name, ".csv")), sep = ",", header = TRUE)
    rnacofolds_seed <- rnacofolds_seed[, !(names(rnacofolds_seed) %in% c("seq_num", "seq_id"))] # drop unneeded columns
    colnames(rnacofolds_seed) <- c("rnacofold_seed_sequence", "rnacofold_seed_struct", "rnacofold_seed_mfe")
    rnacofolds <- cbind(rnacofolds, rnacofolds_seed) # bind the targets and cofolds together

    return(rnacofolds)
}

load_rnaplfolds <- function(experiment_name, folding_windows, seed_features) {
    rnaplfolds <- create_new_frame(c("ensembl_transcript_id_version", "rnaplfold_seed", "rnaplfold_sup"), NULL, nrow(folding_windows))
    rnaplfolds$ensembl_transcript_id_version <- folding_windows$ensembl_transcript_id_version

    for (i in seq_len(nrow(folding_windows))) {
        # rnaplfold's output format is awkward, it stores the files as sequence_0001 to ~sequence_6000 so we have to pad 0s to the filename to keep their format
        padded_zeros <- "000"
        if (i >= 10 && i < 100) {
            padded_zeros <- "00"
        } else if (i >= 100 && i < 1000) {
            padded_zeros <- "0"
        } else if (i >= 1000) {
            padded_zeros <- ""
        }

        rnaplfolds_raw <- read.table(file.path(directories$folds_rnaplfold, experiment_name, paste0("sequence_", padded_zeros, i, "_lunp")), sep = "\t", skip = 2)[-1]

        start_6mer <- folding_windows$rnaplfold_6mer_pos[i]
        start_seed <- start_6mer
        base_9 <- start_6mer + 7
        base_20 <- base_9 + 11

        if (seed_features$seed_binding_type[i] == "7mer-a1" || seed_features$seed_binding_type[i] == "8mer") {
            start_seed <- start_seed - 1 #7mer-a1 / 8mer starts at 1st base
        }


        seed_binding_count <- strtoi(seed_features$seed_binding_count[i])

        rnaplfolds$rnaplfold_seed[i] <- rnaplfolds_raw[start_seed + seed_binding_count, seed_binding_count]
        rnaplfolds$rnaplfold_sup[i] <- rnaplfolds_raw[base_20, 12] # bases 09-20

    return(rnaplfolds)
}

flip_base_pos <- function(pos, len) {
    return(len - pos + 1)
}

mirna_at <- function(pos, seq) {
    return(seq[pos])
}

mrna_at <- function(pos, seq) {
    return(seq[flip_base_pos(pos, length(seq))])
}

is_wobble_paired_by_pos <- function(pos, seq) {
    base_a <- mirna_at(pos, seq)
    base_b <- mrna_at(pos, seq)
    return(
        ((base_a == "U" && base_b == "G") || (base_a == "G" && base_b == "U")) ||
        ((base_a == "A" && base_b == "G") || (base_a == "G" && base_b == "A")) ||
        ((base_a == "U" && base_b == "C") || (base_a == "C" && base_b == "U")) ||
        ((base_a == "A" && base_b == "C") || (base_a == "C" && base_b == "A")))
}

is_paired_by_pos <- function(pos, struct) {
    return(struct[flip_base_pos(pos, length(struct))] == ")")
}

is_perfect_paired_by_struct <- function(pos, struct, seq) {
    return(is_paired_by_pos(pos, struct) && !is_wobble_paired_by_pos(pos, seq))
}

is_wobble_paired_by_base <- function(base_a, base_b) {
    return((base_a == "U" && base_b == "G") || (base_a == "G" && base_b == "U") ||
           (base_a == "A" && base_b == "G") || (base_a == "G" && base_b == "A") ||
           (base_a == "U" && base_b == "C") || (base_a == "C" && base_b == "U") ||
           (base_a == "A" && base_b == "C") || (base_a == "C" && base_b == "A"))
}

is_perfect_paired_by_base <- function(base_a, base_b) {
    return((base_a == "G" && base_b == "C") || (base_a == "C" && base_b == "G") ||
           (base_a == "A" && base_b == "U") || (base_a == "U" && base_b == "A"))
}

extract_seed_features <- function(current_transcript_id, seq, struct) {
    base_count <- length(seq)

    # work out what seed site type we are dealing with
    # aside- 6mer vs 7mer vs 8mer
    # miRNA is 5' -> 3', mRNA is also 5' -> 3'
    # 6mer (mRNA perspective) = ..............||||||.                                 i.e. binding at sites -1 to -6
    # 7mer (mRNA perspective) = ..............||||||A OR .............|||||||         i.e. 6mer plus EITHER an A at 0 OR an extra bind at -7
    # 8mer (mRNA perspective) = .............|||||||A                                 i.e. 6mer plus BOTH an A at 0 AND an extra bind at -7 (both 7mer conditions)
    # we always have at least 6 seed bindings because what we're dealing with wouldn't be considered a target otherwise
    condition_a1 <- seq[base_count] == "A"
    condition_m8 <- if (length(struct) < 8) 0 else as.integer(is_perfect_paired_by_struct(8, struct, seq))

    seed_binding_count <- 6 + condition_a1 + condition_m8
    seed_binding_type <- paste0(seed_binding_count, "mer")
    if (seed_binding_count == 7) {
        if (condition_a1) {
            seed_binding_type <- paste0(seed_binding_type, "-a1")
        } else {
            seed_binding_type <- paste0(seed_binding_type, "-m8")
        }
    }

    perfect_pair_1 <- is_perfect_paired_by_struct(1, struct, seq)

    return(c(current_transcript_id, seed_binding_count, seed_binding_type, perfect_pair_1))
}

extract_single_base_features <- function(current_transcript_id, seq) {

    mirna_1 <- mirna_at(1, seq)
    mrna_1 <- mrna_at(1, seq)
    gu_1 <- (mirna_1 == "G" && mrna_1 == "U") || (mirna_1 == "U" && mrna_1 == "G")


    mirna_8 <- mirna_at(8, seq)
    mrna_8 <- mrna_at(8, seq)
    gu_8 <- (mirna_8 == "G" && mrna_8 == "U") || (mirna_8 == "U" && mrna_8 == "G")

    mirna_9 <- mirna_at(9, seq)
    mrna_9 <- mrna_at(9, seq)
    perfect_pair_9 <- is_perfect_paired_by_base(mirna_9, mrna_9)
    gu_9 <- (mirna_9 == "G" && mrna_9 == "U") || (mirna_9 == "U" && mrna_9 == "G")

    mirna_10 <- mirna_at(10, seq)
    mrna_10 <- mrna_at(10, seq)
    perfect_pair_10 <- is_perfect_paired_by_base(mirna_10, mrna_10)
    gu_10 <- (mirna_10 == "G" && mrna_10 == "U") || (mirna_10 == "U" && mrna_10 == "G")

    return(c(current_transcript_id, gu_1, gu_8, perfect_pair_9, gu_9, perfect_pair_10, gu_10, mirna_1, mirna_8, mrna_8))
}

determine_pairs <- function(struct) {
    # process the mirna in a form such as .((((....(...
    # process the mrna in a form such as ..........)).)))))).

    halfway <- which(struct == "&")

    # first work through the bindings on the mirna side of the fold (the part to the left of the & sign)
    mirna_bindings <- rep(0, halfway - 1)
    base_idx <- 1
    while (base_idx < halfway) {
        if (struct[base_idx] == "(") { # if we have a binding, count it
            mirna_bindings[base_idx] = base_idx
        } else if (as.logical(settings$ignore_second_struct_bind) && struct[base_idx] == ")") { # if we're discounting self binds and we find one, discount the last binding
            mirna_bindings[max(which(mirna_bindings >= 1))] = 0
        }

        base_idx <- base_idx + 1
    }

    # next work through the bindings on the mrna side of the fold (the part to the right of the & sign)
    mrna_bindings <- rep(0, length(struct) - halfway)
    base_idx <- length(struct)
    while (base_idx > halfway) {
        if (struct[base_idx] == ")") { # if we have a binding, count it
            mrna_bindings[base_idx - halfway] <- base_idx
        } else if (as.logical(settings$ignore_second_struct_bind) && struct[base_idx] == "(") { # if we're discounting self binds and we find one, discount the last binding
            mrna_bindings[min(which(mrna_bindings >= 1))] <- 0
        }

        base_idx <- base_idx - 1
    }

    mirna_binds <- which(mirna_bindings >= 1)
    mrna_binds <- which(mrna_bindings >= 1)

    if (length(mrna_binds) == 0 || length(mirna_binds) == 0) {
        return(list(halfway, 0, 0))
    }

    return(list(halfway, mirna_binds, mrna_binds))
}

extract_generic_window_features <- function(current_transcript_id, supplementary_pairs, seq, start, end) {
    halfway <- supplementary_pairs[[1]]
    mirna_binds <- supplementary_pairs[[2]]
    mrna_binds <- supplementary_pairs[[3]]

    if (mirna_binds == 0 && mrna_binds == 0) {
        return(c(current_transcript_id, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
    }

    perfect_pair_count <- 0
    wobble_pair_count <- 0
    any_pair_count <-0

    longest_perfect_sequence <- 0
    longest_any_pair_sequence <- 0
    curr_perfect_sequence <- 0
    curr_any_pair_sequence <- 0
    prev_mirna_index <- 0
    curr_any_sequence_start <- 0
    longest_any_sequence_start <- 0 

    perfect_pair_dist <- 0
    wobble_pair_dist <- 0
    any_pair_dist <- 0
    any_pair_avg_dist <- 0

    gc_count <- 0
    au_count <- 0
    gu_count <- 0
    gc_dist <- 0
    au_dist <- 0
    gu_dist <- 0

    first_base_pos <- 0
    last_base_pos <- 0

    seed_pair_count <- 0
    sup_pair_count <- 0

    for (i in seq_len(length(mirna_binds))) {
        mirna_base_index <- mirna_binds[i]

        if (mirna_base_index < start) {
            next
        } else if (mirna_base_index > end) {
            break
        }

        mrna_base_index <- mrna_binds[length(mrna_binds) - i + 1] # +1 because R counts from 1 rather than 0

        if (abs(mirna_base_index - prev_mirna_index) > 1) {
            curr_perfect_sequence <- 0
            curr_any_pair_sequence <- 0
            curr_any_sequence_start <- 0
        }

        mirna_base <- seq[mirna_base_index]
        mrna_base <- seq[halfway + mrna_base_index]


        if (is_perfect_paired_by_base(mirna_base, mrna_base)) {
            perfect_pair_count <- perfect_pair_count + 1
            any_pair_count <- any_pair_count + 1

            if (mirna_base_index < 9) {
                seed_pair_count <- seed_pair_count + 1
            } else {
                sup_pair_count <- sup_pair_count + 1
            }
            
            if (curr_any_sequence_start == 0) {
                curr_any_sequence_start <- mirna_base_index
            } 

            if (perfect_pair_count > 1) {
                perfect_pair_dist <- perfect_pair_dist + mirna_base_index - prev_mirna_index
            }
            if (any_pair_count > 1) {
                any_pair_dist <- any_pair_dist + mirna_base_index - prev_mirna_index
            }

            curr_perfect_sequence <- curr_perfect_sequence + 1
            if (curr_perfect_sequence > longest_perfect_sequence) {
                longest_perfect_sequence <- curr_perfect_sequence
                longest_any_sequence_start <- curr_any_sequence_start
            }

            curr_any_pair_sequence <- curr_any_pair_sequence + 1
            if (curr_any_pair_sequence > longest_any_pair_sequence) {
                longest_any_pair_sequence <- curr_any_pair_sequence
                longest_any_sequence_start <- curr_any_sequence_start
            }
            
            if (first_base_pos == 0) {
                first_base_pos <- mrna_base_index
            }

            last_base_pos <- mrna_base_index
        } else if (is_wobble_paired_by_base(mirna_base, mrna_base)) {
            wobble_pair_count <- wobble_pair_count + 1
            any_pair_count <- any_pair_count + 1
            
            if (mirna_base_index < 9) {
                seed_pair_count <- seed_pair_count + 1
            } else {
                sup_pair_count <- sup_pair_count + 1
            }

            if (wobble_pair_count > 1) {
                wobble_pair_dist <- wobble_pair_dist + mirna_base_index - prev_mirna_index
            }
            if (any_pair_count > 1) {
                any_pair_dist <- any_pair_dist + mirna_base_index - prev_mirna_index
            }

            curr_perfect_sequence <- 0
            
            if (curr_any_sequence_start == 0) {
                curr_any_sequence_start <- mirna_base_index
            }

            curr_any_pair_sequence <- curr_any_pair_sequence + 1
            if (curr_any_pair_sequence > longest_any_pair_sequence) {
                longest_any_pair_sequence <- curr_any_pair_sequence
                longest_any_sequence_start <- curr_any_sequence_start
            }
            
            if (first_base_pos == 0) {
                first_base_pos <- mrna_base_index
            }

            last_base_pos <- mrna_base_index
        } else {
            curr_perfect_sequence <- 0
            curr_any_pair_sequence <- 0
            curr_any_sequence_start <- 0
        }

        if ((mirna_base == "G" && mrna_base == "C") || (mirna_base == "C" && mrna_base == "G")) {
            gc_count <- gc_count + 1
            gc_dist <- gc_dist + mirna_base_index
        } else if ((mirna_base == "A" && mrna_base == "U") || (mirna_base == "U" && mrna_base == "A")) {
            au_count <- au_count + 1
            au_dist <- au_dist + mirna_base_index
        } else if ((mirna_base == "G" && mrna_base == "U") || (mirna_base == "U" && mrna_base == "G")) {
            gu_count <- gu_count + 1
            gu_dist <- gu_dist + mirna_base_index
        }

        prev_mirna_index <- mirna_base_index
    }

    if (any_pair_dist != 0) {
        any_pair_avg_dist <- any_pair_dist / (any_pair_count - 1)
    }

    mrna_binding_spread <- first_base_pos - last_base_pos
    return(c(current_transcript_id, perfect_pair_count, longest_any_pair_sequence, longest_any_sequence_start, any_pair_avg_dist, gu_count, mrna_binding_spread))
}

extract_au_content_features <- function(current_transcript_id, window_lr, window_rl) {
    au_window_rl <- substr(window_rl, 9, nchar(window_rl)) # counts from 1, so its a 9 rather than 8
    au_window_lr <- substr(window_lr, 1, nchar(window_lr) - 8)
    au_window_sup <- substr(window_rl, 11, 18)

    au_content_3 <- lengths(regmatches(au_window_lr, gregexpr("(A|U)", au_window_lr)))
    au_content_sup <- lengths(regmatches(au_window_sup, gregexpr("(A|U)", au_window_sup)))

    au_content_3 <- au_content_3 / nchar(au_window_lr)
    au_content_sup <- au_content_sup / nchar(au_window_sup)



    weights <- c(1/3, 1/4, 1/5, 1/6, 1/7, 1/8, 1/9, 1/10, 1/11, 1/12, 1/13, 1/14, 1/15, 1/16, 1/17, 1/18, 1/19, 1/20, 1/21, 1/22, 1/23, 1/24, 1/25, 1/26, 1/27, 1/28, 1/29, 1/30, 1/31, 1/32)
    dist_scalars <- weights / sum(weights)

    au_content_5_weighted <- 0
    au_5_match_idx <- gregexpr("(A|U)", au_window_rl)[[1]]
    for (idx in au_5_match_idx) {
        au_content_5_weighted <- au_content_5_weighted + dist_scalars[idx]
    }
    
    
    return(c(current_transcript_id, au_content_3, au_content_sup, au_content_5_weighted))
}



# --- ENTRY POINT ---

annotations <- read.table(file.path(directories$annotations, "annotations.tsv"), sep = "\t", header = TRUE)

experiment_files <- dir(directories$windows, pattern = ".tsv")
experiment_names <- gsub("*.tsv", "", experiment_files)

for (i in seq_len(length(experiment_files))) {
    cat(paste0("Feature extraction ", i, "/", length(experiment_files), "... "))

    experiment_name <- experiment_names[i]
    experiment_filename <- experiment_files[i]
    feature_path <- file.path(directories$features, experiment_filename)

    if (as.logical(settings$use_caching) && file.exists(feature_path)) {
        cat("Loaded from cache.\n")
    } else {
        expanded_binding_sites <- read.table(file.path(directories$bindings, experiment_filename), sep = "\t", header = TRUE)

        
        rnafolds_lr_raw <- read.table(file.path(directories$folds_rnafold_lr, paste0(experiment_name, ".csv")), sep = ",", header = FALSE)
        rnafolds_lr <- do.call("cbind", split(rnafolds_lr_raw, rep(c(1, 2), length.out = nrow(rnafolds_lr_raw))))
        colnames(rnafolds_lr) <- c("mirna_sequence", "rnafold_struct")

        rnafolds_rl_raw <- read.table(file.path(directories$folds_rnafold_rl, paste0(experiment_name, ".csv")), sep = ",", header = FALSE)
        rnafolds_rl <- do.call("cbind", split(rnafolds_rl_raw, rep(c(1, 2), length.out = nrow(rnafolds_rl_raw))))
        colnames(rnafolds_rl) <- c("mirna_sequence", "rnafold_struct")


        rnafolds <- load_rnafolds(experiment_name, expanded_binding_sites, rnafolds_lr, rnafolds_rl)
        rnacofolds <- load_rnacofolds(experiment_name, expanded_binding_sites)

        seed_features <- create_new_frame(
            c("ensembl_transcript_id_version", "seed_binding_count", "seed_binding_type", "perfect_pair_1"), 
            NULL, nrow(expanded_binding_sites))

        single_base_features <- create_new_frame(
            c("ensembl_transcript_id_version", "gu_1", "gu_8", "perfect_pair_9", "gu_9", "perfect_pair_10", "gu_10", "mirna_1", "mirna_8", "mrna_8"),
            NULL, nrow(expanded_binding_sites))

        window_frame_columns <- c("ensembl_transcript_id_version", "perfect_pair_count", "longest_any_sequence", "longest_any_sequence_start", "any_pair_avg_dist", "gu_count", "mrna_binding_spread")
        window_features_12_17 <- create_new_frame(window_frame_columns, NULL, nrow(expanded_binding_sites))
        window_features_09_20 <- create_new_frame(window_frame_columns, NULL, nrow(expanded_binding_sites))
        window_features_full <- create_new_frame(window_frame_columns, NULL, nrow(expanded_binding_sites))
        colnames(window_features_12_17) <- paste(colnames(window_features_12_17), "12_17", sep = "_")
        colnames(window_features_09_20) <- paste(colnames(window_features_09_20), "09_20", sep = "_")
        colnames(window_features_full) <- paste(colnames(window_features_full), "full", sep = "_")
        colnames(window_features_12_17)[1] <- "ensembl_transcript_id_version"
        colnames(window_features_09_20)[1] <- "ensembl_transcript_id_version"
        colnames(window_features_full)[1] <- "ensembl_transcript_id_version"

        au_content_features <- create_new_frame(
            c("ensembl_transcript_id_version", "au_content_3", "au_content_sup", "au_content_5_weighted"),
            NULL, nrow(expanded_binding_sites))

        for (j in seq_len(nrow(expanded_binding_sites))) {
            current_transcript_id <- expanded_binding_sites$ensembl_transcript_id_version[j]

            seq <- strsplit(as.character(rnacofolds$rnacofold_full_sequence[j]), "")[[1]]
            struct <- strsplit(as.character(rnacofolds$rnacofold_full_struct[j]), "")[[1]]

            seed_features[j, ] <- extract_seed_features(current_transcript_id, seq, struct)
            single_base_features[j, ] <- extract_single_base_features(current_transcript_id, seq)

            sequence_pairs <- determine_pairs(struct)
            
            au_content_features[j, ] <- extract_au_content_features(current_transcript_id, rnafolds_lr[j, 1], rnafolds_rl[j, 1])
            window_features_12_17[j, ] <- extract_generic_window_features(current_transcript_id, sequence_pairs, seq, 12, 17)
            window_features_09_20[j, ] <- extract_generic_window_features(current_transcript_id, sequence_pairs, seq, 09, 20)
            window_features_full[j, ] <- extract_generic_window_features(current_transcript_id, sequence_pairs, seq, 01, 99)
        }

        folding_windows <- read.table(file.path(directories$windows, experiment_filename), sep = "\t", header = TRUE)
        rnaplfolds <- load_rnaplfolds(experiment_name, folding_windows, seed_features)

        combined_features <- cbind(
            single_base_features,
            seed_features[, -which(names(seed_features) == "ensembl_transcript_id_version")],
            window_features_12_17[, -which(names(window_features_12_17) == "ensembl_transcript_id_version")],
            window_features_09_20[, -which(names(window_features_09_20) == "ensembl_transcript_id_version")],
            window_features_full[, -which(names(window_features_09_20) == "ensembl_transcript_id_version")],
            expanded_binding_sites[, -which(names(expanded_binding_sites) == "ensembl_transcript_id_version")],
            au_content_features[, -which(names(au_content_features) == "ensembl_transcript_id_version")],
            rnafolds[, c("rnafold_mfe", "rnafold_direction")],
            rnacofolds[, c("rnacofold_full_mfe", "rnacofold_seed_mfe")],
            rnaplfolds[, -which(names(rnaplfolds) == "ensembl_transcript_id_version")]
        )

        # combine utr and cds length data
        combined_features <- merge(combined_features, annotations[, c("ensembl_transcript_id_version", "X3_utr_length", "cds_length")], by = "ensembl_transcript_id_version")

        # determine the closest utr end distance
        combined_features$dist_closest_utr_end <- -1
        for (j in seq_len(nrow(combined_features))) {
            end <- combined_features$X3_utr_length[j]
            binding_pos <- combined_features$binding_site_pos[j]

            dist_to_end <- end - binding_pos

            if (dist_to_end < binding_pos) {
                combined_features$dist_closest_utr_end[j] <- dist_to_end
            } else {
                combined_features$dist_closest_utr_end[j] <- binding_pos
            }
        }

        # for cases of abundance, track which is likely the strongest candidate (most seed bases paired > strongest mfe)
        combined_features$best_abundance <- FALSE
        best_abundance_idx <- -1
        best_abundance_mfe <- 0
        best_abundance_seed_bp <- 6
        prev_transcript_id <- combined_features$ensembl_transcript_id_version[1]
        for (j in seq_len(nrow(combined_features))) {
            abundance <- combined_features$site_abundance_6mer[j]
            if (abundance == 1) {
                combined_features$best_abundance[j] <- TRUE

                if (best_abundance_idx > -1) {
                    combined_features$best_abundance[best_abundance_idx] <- TRUE
                    best_abundance_idx <- -1
                    best_abundance_mfe <- 0
                    best_abundance_seed_bp <- 6
                }
            } else {
                cur_transcript_id <- combined_features$ensembl_transcript_id_version[j]

                if (prev_transcript_id != cur_transcript_id) {
                    if (best_abundance_idx > -1) {
                        combined_features$best_abundance[best_abundance_idx] <- TRUE
                        best_abundance_idx <- -1
                        best_abundance_mfe <- 0
                        best_abundance_seed_bp <- 6
                    }
                }

                seed_bp <- combined_features$seed_binding_count[j]
                mfe <- combined_features$rnacofold_full_mfe[j]

                if (as.numeric(seed_bp) > as.numeric(best_abundance_seed_bp)) {
                    best_abundance_idx <- j
                    best_abundance_mfe <- mfe
                    best_abundance_seed_bp <- seed_bp
                } else if (as.numeric(seed_bp) == as.numeric(best_abundance_seed_bp)) {
                    if (as.numeric(mfe) < as.numeric(best_abundance_mfe)) {
                        best_abundance_idx <- j
                        best_abundance_mfe <- mfe
                    }
                }

                prev_transcript_id <- cur_transcript_id
            }
        }

        if (best_abundance_idx > -1) {
            combined_features$best_abundance[best_abundance_idx] <- TRUE
        }

        # discount the binding of the transcript itself from abundance
        for (j in seq_len(nrow(combined_features))) {
            if (seed_features$seed_binding_type[j] == "6mer") {
                seed_features$site_abundance_6mer[j] <- seed_features$site_abundance_6mer[j] - 1
            } else if (seed_features$seed_binding_type[j] == "7mer-a1") {
                seed_features$site_abundance_7a1[j] <- seed_features$site_abundance_7a1[j] - 1
                seed_features$site_abundance_7mer[j] <- seed_features$site_abundance_7mer[j] - 1
            } else if (seed_features$seed_binding_type[j] == "7mer-m8") {
                seed_features$site_abundance_7m8[j] <- seed_features$site_abundance_7m8[j] - 1
                seed_features$site_abundance_7mer[j] <- seed_features$site_abundance_7mer[j] - 1
            } else if (seed_features$seed_binding_type[j] == "8mer") {
                seed_features$site_abundance_8mer[j] <- seed_features$site_abundance_8mer[j] - 1
            }
        }

        # Output features
        write.table(combined_features, feature_path, quote = FALSE, sep = "\t", row.names = FALSE)

        cat("Done.\n")
    }
}