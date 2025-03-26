args <- commandArgs(trailingOnly = TRUE)
suppressMessages(library(jsonlite))
config <- jsonlite::fromJSON(args[1])
settings <- config$settings
directories <- config$directories

suppressMessages(library(dplyr))
suppressMessages(library(parallel))
suppressMessages(library(stringi))


# --- FUNCTIONS ---

create_new_frame <- dget("src/functions/create_new_frame.r")
reverse_complement <- dget("src/functions/reverse_complement.r")

# find any 6mer or 6mer offset binding sites and count them up to determine their abundance values, also check the cds
locate_binding_sites <- function(utrs, cds, target_6mer, target_6off, target_7mer_a1, target_7mer_m8, target_8mer) {
    binding_sites <- create_new_frame(
        c("ensembl_transcript_id_version", 
        "site_abundance_6mer", "site_abundance_6off", "site_abundance_7a1", "site_abundance_7m8", "site_abundance_7mer", "site_abundance_8mer",
        "site_abundance_6cds", "site_abundance_7a1cds", "site_abundance_7m8cds", "site_abundance_7cds", "site_abundance_8cds"), NULL, nrow(utrs))
    
    binding_sites$site_abundance_6mer <- stri_count_regex(utrs$X3utr, paste0("(?=", target_6mer, ")"))
    binding_sites$site_abundance_6off <- stri_count_regex(utrs$X3utr, paste0("(?=", target_6off, ")"))
    binding_sites$site_abundance_7a1 <- stri_count_regex(utrs$X3utr, paste0("(?=", target_7mer_a1, ")"))
    binding_sites$site_abundance_7m8 <- stri_count_regex(utrs$X3utr, paste0("(?=", target_7mer_m8, ")"))
    binding_sites$site_abundance_8mer <- stri_count_regex(utrs$X3utr, paste0("(?=", target_8mer, ")"))
    binding_sites$site_abundance_6cds <- stri_count_regex(cds$cds, paste0("(?=", target_6mer, ")"))
    binding_sites$site_abundance_7a1cds <- stri_count_regex(cds$cds, paste0("(?=", target_7mer_a1, ")"))
    binding_sites$site_abundance_7m8cds <- stri_count_regex(cds$cds, paste0("(?=", target_7mer_m8, ")"))
    binding_sites$site_abundance_8cds <- stri_count_regex(cds$cds, paste0("(?=", target_8mer, ")"))

    binding_sites$site_abundance_7mer <- binding_sites$site_abundance_7a1 + binding_sites$site_abundance_7m8
    binding_sites$site_abundance_7cds <- binding_sites$site_abundance_7a1cds + binding_sites$site_abundance_7m8cds
    binding_sites$ensembl_transcript_id_version <- utrs$ensembl_transcript_id_version

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
    expanded_binding_sites <- binding_sites[rep(1:nrow(binding_sites), binding_sites$site_abundance_6mer), ]
    filtered_utrs <- utrs[utrs$ensembl_transcript_id_version %in% binding_sites$ensembl_transcript_id_version, ]
    filtered_utrs <- filtered_utrs[rep(1:nrow(binding_sites), binding_sites$site_abundance_6mer), ]

    i <- 1
    while (i <= nrow(expanded_binding_sites)) {
        # note: using lookahead assertion (?=pattern) regex to find overlapping matches, e.g. GCAGCAGCA should return a match for GCAGCA at both pos 1 and 4, not just 1 
        matches_6mer <- unlist(gregexpr(paste0("(?=", target_6mer, ")"), filtered_utrs$X3utr[i], perl = TRUE))
        for (j in seq_len(length(matches_6mer))) {
            expanded_binding_sites$binding_site_pos[i + j - 1] <- matches_6mer[j]
        }
        i <- i + j
    }

    return(expanded_binding_sites)
}

process_mirna <- function(i) {

    if (!(as.logical(settings$use_caching) && file.exists(file.path(directories$bindings, paste0(mirna_sequences[i, ]$mirna_id, '.tsv'))))) {
        
        start_time <- Sys.time()
        
        mirna <- mirna_sequences[i, ]
        mirna_id <- mirna$mirna_id
        mirna_sequence <- mirna$mirna_sequence

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

        # handle abundance as separate sites
        expanded_binding_sites <- expand_binding_sites(binding_sites, utrs, target_6mer) 

        # only log binding sites for this transfection if there was at least one 6mer or better
        if (nrow(expanded_binding_sites) > 0) {
            write.table(
                binding_sites[order(binding_sites$ensembl_transcript_id_version), ],
                file = file.path(directories$bindings_raw, paste0(mirna_id, ".tsv")), quote = FALSE, sep = "\t",
                row.names = FALSE)

            write.table(
                expanded_binding_sites[order(expanded_binding_sites$ensembl_transcript_id_version), ],
                file = file.path(directories$bindings, paste0(mirna_id, ".tsv")), quote = FALSE, sep = "\t",
                row.names = FALSE)

            # log binding site and target information for this transfection
            target_sites[i, ] <- c(gsub("*.tsv", "", mirna_id), site_6mer, target_6mer, site_6off, target_6off, site_7mer_a1, target_7mer_a1, 
                site_7mer_m8, target_7mer_m8, site_8mer, target_8mer)
        }

        message(paste0("Locating binding sites ", i, "/", n, " - done.", "-- Elapsed-- ", Sys.time() - start_time))
    } else {
        message(paste0("Locating binding sites ", i, "/", n, " - loaded from cache"))
    }

    write.table(target_sites, file.path(directories$bindings, "target-sites.tsv"), quote = FALSE, sep = "\t", row.names = FALSE)
}


# --- ENTRY POINT ---

# get utr and cds sequences for all data with a l2fc value
annotations <- read.table(file.path(directories$annotations, "annotations.tsv"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
mane <- read.table(file.path(directories$annotations, "mane.tsv"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
utrs <- read.table(file.path(directories$annotations, "utr_sequences.tsv"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
cds <- read.table(file.path(directories$annotations, "cds_sequences.tsv"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)

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
cds <- cds[cds$ensembl_transcript_id_version %in% annotations$ensembl_transcript_id_version, ]
utrs <- utrs[utrs$X3utr != "Sequence unavailable", ]
cds <- cds[cds$ensembl_transcript_id_version %in% utrs$ensembl_transcript_id_version, ]

mirna_sequences <- read.table(file.path(directories$annotations, "mirna_sequences.tsv"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
if (settings$mirna_id_filter != "") {
    mirna_sequences <- mirna_sequences[mirna_sequences$mirna_id %in% strsplit(settings$mirna_id_filter, ",")[[1]], ]
}

target_sites <- create_new_frame(c("mirna_id", "site_6mer", "target_6mer", "site_6off", "target_6off", "site_7mer_a1", "target_7mer_a1", "site_7mer_m8",
    "target_7mer_m8", "site_8mer", "target_8mer"), NULL, nrow(mirna_sequences))


n <- nrow(mirna_sequences)
max_cores <- as.numeric(settings$max_cores)
cores <- ifelse(max_cores == -1, detectCores() - 1, cores <- max_cores)

result <- mclapply(1:n, process_mirna, mc.cores = cores)
