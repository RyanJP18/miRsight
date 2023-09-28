args <- commandArgs(trailingOnly = TRUE)
suppressMessages(library(jsonlite))
config <- jsonlite::fromJSON(args[1])
settings <- config$settings
directories <- config$directories

suppressMessages(library(dplyr))



# --- FUNCTIONS ---

create_new_frame <- dget("src/functions/create_new_frame.r")
reverse_complement <- dget("src/functions/reverse_complement.r")

# find any 6mer or 6mer offset binding sites and count them up to determine their abundance values, also check the cds
locate_binding_sites <- function(utrs, cds, target_6mer, target_6off, target_7mer_a1, target_7mer_m8, target_8mer) {
    binding_sites <- create_new_frame(
        c("ensembl_transcript_id_version", 
        "site_abundance_6mer", "site_abundance_6off", "site_abundance_7a1", "site_abundance_7m8", "site_abundance_7mer", "site_abundance_8mer",
        "site_abundance_6cds", "site_abundance_7a1cds", "site_abundance_7m8cds", "site_abundance_7cds", "site_abundance_8cds"), NULL, nrow(utrs))
    binding_sites$ensembl_transcript_id_version <- utrs$ensembl_transcript_id_version

    for (i in seq_len(nrow(binding_sites))) {
        # Find any sequences matches within the current utr which match our seed target sequence, multiple is possible
        matches_6mer <- unlist(gregexpr(paste0("(?=", target_6mer, ")"), utrs$X3utr[i], perl = TRUE))
        site_abundance_6mer <- length(matches_6mer)
        if (site_abundance_6mer == 1 && matches_6mer == -1) {
            site_abundance_6mer <- 0
        }
        binding_sites$site_abundance_6mer[i] <- site_abundance_6mer

        # Find any 6mer offset potential binding sites
        matches_6off <- unlist(gregexpr(paste0("(?=", target_6off, ")"), utrs$X3utr[i], perl = TRUE))
        site_abundance_6off <- length(matches_6off)
        if (site_abundance_6off == 1 && matches_6off == -1) {
            site_abundance_6off <- 0
        }
        binding_sites$site_abundance_6off[i] <- site_abundance_6off

        matches_7a1 <- unlist(gregexpr(paste0("(?=", target_7mer_a1, ")"), utrs$X3utr[i], perl = TRUE))
        site_abundance_7a1 <- length(matches_7a1)
        if (site_abundance_7a1 == 1 && matches_7a1 == -1) {
            site_abundance_7a1 <- 0
        }
        binding_sites$site_abundance_7a1[i] <- site_abundance_7a1

        matches_7m8 <- unlist(gregexpr(paste0("(?=", target_7mer_m8, ")"), utrs$X3utr[i], perl = TRUE))
        site_abundance_7m8 <- length(matches_7m8)
        if (site_abundance_7m8 == 1 && matches_7m8 == -1) {
            site_abundance_7m8 <- 0
        }
        binding_sites$site_abundance_7m8[i] <- site_abundance_7m8
        binding_sites$site_abundance_7mer[i] <- site_abundance_7a1 + site_abundance_7m8

        matches_8mer <- unlist(gregexpr(paste0("(?=", target_8mer, ")"), utrs$X3utr[i], perl = TRUE))
        site_abundance_8mer <- length(matches_8mer)
        if (site_abundance_8mer == 1 && matches_8mer == -1) {
            site_abundance_8mer <- 0
        }
        binding_sites$site_abundance_8mer[i] <- site_abundance_8mer

        # repeat for cds/orf
        matches_6cds <- unlist(gregexpr(paste0("(?=", target_6mer, ")"), cds$cds[i], perl = TRUE))
        site_abundance_6cds <- length(matches_6cds)
        if (site_abundance_6cds == 1 && matches_6cds == -1) {
            site_abundance_6cds <- 0
        }
        binding_sites$site_abundance_6cds[i] <- site_abundance_6cds

        matches_7a1cds <- unlist(gregexpr(paste0("(?=", target_7mer_a1, ")"), cds$cds[i], perl = TRUE))
        site_abundance_7a1cds <- length(matches_7a1cds)
        if (site_abundance_7a1cds == 1 && matches_7a1cds == -1) {
            site_abundance_7a1cds <- 0
        }
        binding_sites$site_abundance_7a1cds[i] <- site_abundance_7a1cds

        matches_7m8cds <- unlist(gregexpr(paste0("(?=", target_7mer_m8, ")"), cds$cds[i], perl = TRUE))
        site_abundance_7m8cds <- length(matches_7m8cds)
        if (site_abundance_7m8cds == 1 && matches_7m8cds == -1) {
            site_abundance_7m8cds <- 0
        }
        binding_sites$site_abundance_7m8cds[i] <- site_abundance_7m8cds
        binding_sites$site_abundance_7cds[i] <- site_abundance_7a1cds + site_abundance_7m8cds

        matches_8cds <- unlist(gregexpr(paste0("(?=", target_8mer, ")"), cds$cds[i], perl = TRUE))
        site_abundance_8cds <- length(matches_8cds)
        if (site_abundance_8cds == 1 && matches_8cds == -1) {
            site_abundance_8cds <- 0
        }
        binding_sites$site_abundance_8cds[i] <- site_abundance_8cds
    }

    return(binding_sites)
}

# duplicate any transcripts that were found to have a binding site by the number of times a match was found for it, e.g. 3 sites = 3 rows, so we can process each separately
expand_binding_sites <- function(binding_sites, utrs, target_6mer) {
    # expand out the binding sites and apply the same transformation to the utrs for easy lookup
    # note: we're essentially duplicating any transcript that have more than 1 binding site so we can process them as separate entities
    #--- what we start with
    # transcript1 total_matches=3
    # transcript2 total_matches=1
    #--- what we actually want
    # transcript1 match1 total_matches=3
    # transcript1 match2 total_matches=3
    # transcript1 match3 total_matches=3
    # transcript2 match1 total_matches=1
    expanded_binding_sites <- as.data.frame(lapply(binding_sites, rep, binding_sites$site_abundance_6mer))
    filtered_utrs <- utrs[utrs$ensembl_transcript_id_version %in% binding_sites$ensembl_transcript_id_version, ]
    filtered_utrs <- lapply(filtered_utrs, rep, binding_sites$site_abundance_6mer)

    i <- 1
    while (i <= nrow(expanded_binding_sites)) {
        matches_6mer <- unlist(gregexpr(paste0("(?=", target_6mer, ")"), filtered_utrs$X3utr[i], perl = TRUE))
        for (j in seq_len(length(matches_6mer))) {
            expanded_binding_sites$binding_site_pos[i + j - 1] <- matches_6mer[j]
        }
        i <- i + j
    }

    return(expanded_binding_sites)
}


# --- ENTRY POINT ---

# get utr and cds sequences for all data with a l2fc value
annotations <- read.table(file.path(directories$annotations, "annotations.tsv"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
mane <- read.table(file.path(directories$annotations, "mane.tsv"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
utrs <- read.table(file.path(directories$annotations, "utr_sequences.tsv"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
cds <- read.table(file.path(directories$annotations, "cds_sequences.tsv"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
# apply chromosome filter and also ensure we have annotations for the transcripts
annotations <- annotations[annotations$chromosome_name %in% strsplit(settings$chromosome_filter, ",")[[1]], ]
annotations <- annotations[annotations$ensembl_transcript_id_version %in% mane$ensembl_transcript_id_version, ]
utrs <- utrs[utrs$ensembl_transcript_id_version %in% annotations$ensembl_transcript_id_version, ]
cds <- cds[cds$ensembl_transcript_id_version %in% annotations$ensembl_transcript_id_version, ]
# filter to only transcripts where we have utr sequence annotations
utrs <- utrs[utrs$X3utr != "Sequence unavailable", ]
cds <- cds[cds$ensembl_transcript_id_version %in% utrs$ensembl_transcript_id_version, ]



mirna_sequences <- read.table(file.path(directories$annotations, "mirna_sequences.tsv"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)

if (file.exists(file.path(directories$bindings, "target-sites.tsv"))) {
    target_sites <- read.table(file.path(directories$bindings, "target-sites.tsv"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
} else {
    target_sites <- create_new_frame(c("mirna_id", "site_6mer", "target_6mer", "site_6off", "target_6off", "site_7mer_a1", "target_7mer_a1", "site_7mer_m8",
        "target_7mer_m8", "site_8mer", "target_8mer"), NULL, nrow(mirna_sequences))
}

for (i in seq_len(nrow(mirna_sequences))) {
    mirna <- mirna_sequences[i, ]
    mirna_id <- mirna$mirna_id
    mirna_sequence <- mirna$mirna_sequence

    cat(paste0("Computing target sites ", i, "/", nrow(mirna_sequences), "..."))

    if (!(as.logical(settings$use_caching) && file.exists(file.path(directories$bindings, paste0(mirna_id, '.tsv'))))) {
        # find the current miRNA sequence, search for its 6mer site and RC it to get the 6mer target site
        site_6mer <- substr(mirna_sequence, start = 2, stop = 7) # miRNA is fetched from 5` -> 3`
        target_6mer <- reverse_complement(site_6mer)

        # find the current miRNA sequence, search for its 6mer offset site and RC it to get the 6mer target site
        site_6off <- substr(mirna_sequence, start = 3, stop = 8) # miRNA is fetched from 5` -> 3`
        target_6off <- reverse_complement(site_6off)

        # expand out for 7mers and 8mers
        site_7mer_a1 <- paste0("T", substr(mirna_sequence, start = 2, stop = 7))
        target_7mer_a1 <- reverse_complement(site_7mer_a1)
        site_7mer_m8 <- substr(mirna_sequence, start = 2, stop = 8)
        target_7mer_m8 <- reverse_complement(site_7mer_m8)
        site_8mer <- paste0("T", substr(mirna_sequence, start = 2, stop = 8)) 
        target_8mer <- reverse_complement(site_8mer)

        # locate the binding sites and store them
        binding_sites <- locate_binding_sites(utrs, cds, target_6mer, target_6off, target_7mer_a1, target_7mer_m8, target_8mer)
        write.table(
            binding_sites[order(binding_sites$ensembl_transcript_id_version), ],
            file = file.path(directories$bindings_raw, paste0(mirna_id, ".tsv")), quote = FALSE, sep = "\t",
            row.names = FALSE)

        expanded_binding_sites <- expand_binding_sites(binding_sites, utrs, target_6mer) # handle abundance as separate sites
        write.table(
            expanded_binding_sites[order(expanded_binding_sites$ensembl_transcript_id_version), ],
            file = file.path(directories$bindings, paste0(mirna_id, ".tsv")), quote = FALSE, sep = "\t",
            row.names = FALSE)

        # log binding site and target information for this transfection
        target_sites[i, ] <- c(gsub("*.tsv", "", mirna_id), site_6mer, target_6mer, site_6off, target_6off, site_7mer_a1, target_7mer_a1, 
            site_7mer_m8, target_7mer_m8, site_8mer, target_8mer)

        cat(" Done.\n")
    } else {
        cat(" Loaded from cache.\n")
    }

    write.table(target_sites, file.path(directories$bindings, "target-sites.tsv"), quote = FALSE, sep = "\t", row.names = FALSE)
}