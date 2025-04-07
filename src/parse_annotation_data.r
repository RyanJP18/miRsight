args <- commandArgs(trailingOnly = TRUE)
suppressMessages(library(jsonlite))
config <- jsonlite::fromJSON(args[1])
settings <- config$settings
directories <- config$directories

suppressMessages(library(stringr))
suppressMessages(library(dplyr))



# --- FUNCTIONS ---

create_new_frame <- dget("src/functions/create_new_frame.r")

parse_mirnas <- function() {
    if (as.logical(settings$use_caching) && file.exists(file.path(directories$annotations, "mirna_sequences.tsv"))) {
        cat("Loaded miRNA sequence lookup from cache.\n")
    } else {
        cat("Generating miRNA sequence lookup from miRBase data...")

        raw_mirnas <- read.table(file.path(directories$preload_data, "mirnas.fa"), sep = "\t", header = FALSE, stringsAsFactors = FALSE)

        # extract the sequences
        sequences <- raw_mirnas[seq(2, nrow(raw_mirnas), 2), ]

        # extract the ids from a string of irrelevant information
        raw_ids <- raw_mirnas[seq(1, nrow(raw_mirnas), 2), ]
        ids <- stringr::str_split_fixed(raw_ids, " ", 3)[, 1] # strip everything except the mirna id
        ids <- stringr::str_split_fixed(ids, ">", 2)[, 2] # strip the > prefix

        # compile the ids and sequences into a two-column datatable and log to file
        mirna_sequences <- create_new_frame(c("mirna_id", "mirna_sequence"), NULL, length(ids))
        mirna_sequences$mirna_id <- ids
        mirna_sequences$mirna_sequence <- sequences
        mirna_sequences <- mirna_sequences[grepl("hsa-", mirna_sequences$mirna_id), ]

        # store lookup table
        write.table(mirna_sequences, file.path(directories$annotations, "mirna_sequences.tsv"), quote = FALSE, sep = "\t", row.names = FALSE)

        cat(" Done.\n")
    }
}

parse_annotations <- function() {
    if (as.logical(settings$use_caching) && file.exists(file.path(directories$annotations, "annotations.tsv"))) {
        cat("Loaded annotations from cache.\n")
    } else {
        cat("Parsing annotations from gtf file...")

        # load annotations
        annotations <- read.table(file.path(directories$preload_data, "annotation.gtf"), sep = "\t", header = FALSE, stringsAsFactors = FALSE)

        # get just transcript annotations
        colnames(annotations) <- c("chromosome_name", "source", "feature", "transcript_start", "transcript_end", "score", "strand", "frame", "attribute")
        annotations <- annotations[annotations$feature == "transcript", ]

        # bind attribute property to its own table as it's poorly formatted
        attribute <- as.data.frame(stringr::str_split_fixed(annotations$attribute, ";", 16)[, 1:10])

        # simplify main table to just useful information
        annotations <- annotations[, c("chromosome_name", "transcript_start", "transcript_end", "strand")]

        # parse the poorly formatted attribute data into columns and simplify
        colnames(attribute) <-
            c(
                "ensembl_gene_id", "ensembl_gene_version", "ensembl_transcript_id", "ensembl_transcript_version",
                "external_gene_id", "gene_source", "gene_biotype", "transcript_name", "transcript_source", "transcript_biotype"
            )
        attribute <- attribute[, c("ensembl_gene_id", "ensembl_transcript_id", "ensembl_transcript_version", "external_gene_id", "transcript_biotype")]
        attribute[] <- lapply(attribute, function(x) stringr::str_extract(x, "\\w+$"))

        # merge transcript id columns
        attribute$ensembl_transcript_id_version <- paste(attribute$ensembl_transcript_id, attribute$ensembl_transcript_version, sep = ".")
        attribute <- attribute[, c("ensembl_gene_id", "ensembl_transcript_id_version", "external_gene_id", "transcript_biotype")]

        # bind the two tables together and remove non-protein-coding transcripts
        annotations <- cbind(annotations, attribute)
        annotations <- annotations[annotations$transcript_biotype == "protein_coding", ]
        annotations$transcript_biotype <- NULL

        # order by transcript id and store
        annotations <- annotations[order(annotations$ensembl_transcript_id_version), ]
        annotations <- annotations[, c("ensembl_transcript_id_version", "ensembl_gene_id", "external_gene_id", "chromosome_name", "strand", "transcript_start", "transcript_end")]
        write.table(annotations, file.path(directories$annotations, "annotations.tsv"), quote = FALSE, sep = "\t", row.names = FALSE)

        cat(" Done.\n")
    }
}

parse_mane <- function() {
    if (as.logical(settings$use_caching) && file.exists(file.path(directories$annotations, "mane.tsv"))) {
        cat("Loaded mane from cache.\n")
    } else {
        cat("Parsing mane file...")

        # load mane data
        mane <- read.table(file.path(directories$preload_data, "mane.gtf"), sep = "\t", header = FALSE, stringsAsFactors = FALSE)

        mane <- mane[mane[, 3] == "transcript", ]
        mane <- as.data.frame(stringr::str_split_fixed(mane[, 9], ";", 4)[, 1:2])
        colnames(mane)[1:2] <- c("ensembl_gene_id", "ensembl_transcript_id_version")

        # strip unnnecessary strings, spaces and version from gene ids
        mane$ensembl_transcript_id_version <- gsub("^.{0,1}", "", mane$ensembl_transcript_id_version)
        mane[] <- lapply(mane, function(x) stringr::word(x, 2))
        mane$ensembl_gene_id <- stringr::str_extract(mane$ensembl_gene_id, "[^.]+")

        # order by transcript id and store
        mane <- mane[order(mane$ensembl_transcript_id_version), ]
        mane <- mane[, c("ensembl_transcript_id_version", "ensembl_gene_id")]
        write.table(mane, file.path(directories$annotations, "mane.tsv"), quote = FALSE, sep = "\t", row.names = FALSE)

        cat(" Done.\n")
    }
}



# --- ENTRY POINT ---

parse_mirnas()
parse_annotations()
parse_mane()
