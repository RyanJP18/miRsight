args <- commandArgs(trailingOnly = TRUE)
suppressMessages(library(jsonlite))
config <- jsonlite::fromJSON(args[1])
settings <- config$settings
directories <- config$directories



# --- FUNCTIONS ---

extract_utr_sequences <- function(annotations, ensdb, genome, chromosome_filter) {
    if (as.logical(settings$use_caching) && file.exists(file.path(directories$annotations, "utr_sequences.tsv"))) {
        cat(paste("Loaded 3' utr sequences from cache.\n"))
    } else {
        cat(paste("Extracting 3' utr sequences with EnsemblDB..."))

        utrs <- ensembldb::threeUTRsByTranscript(ensdb, filter = ~ tx_biotype == "protein_coding" & seq_name == chromosome_filter)
        utr_seqs <- GenomicFeatures::extractTranscriptSeqs(genome, utrs)

        utr_table <- data.frame("X3utr" = as.character(utr_seqs), "ensembl_transcript_id" = names(utr_seqs), row.names = NULL)
        utr_table <- merge(utr_table, annotations[, c("ensembl_transcript_id_version", "ensembl_transcript_id")], by = "ensembl_transcript_id", all.y = TRUE)
        utr_table[is.na(utr_table$X3utr), ]$X3utr <- "Sequence unavailable"

        utr_table <- utr_table[, c("ensembl_transcript_id_version", "X3utr")]
        utr_table <- utr_table[order(utr_table$ensembl_transcript_id_version), ]
        write.table(utr_table, file.path(directories$annotations, "utr_sequences.tsv"), quote = FALSE, sep = "\t", row.names = FALSE)

        cat(" Done.\n")
    }
}

extract_sequence_coordinates <- function(annotations, ensdb, genome, chromosome_filter) {
    if (as.logical(settings$use_caching) && file.exists(file.path(directories$annotations, "utr_coords.rds"))) {
        cat(paste("Loaded utr coordinates from cache.\n"))
    } else {
        cat(paste("Extracting utr coordinates with EnsemblDB..."))
    
        filt <- ~ tx_biotype == "protein_coding" & seq_name == chromosome_filter

        utr_infos <- ensembldb::transcripts(ensdb, filter = filt)
        utr_infos_table <- data.frame(utr_infos)
        utr_infos_table <- utr_infos_table[order(utr_infos_table$tx_name), ]

        # outer genomic bounds of the 3'UTR
        plus <- utr_infos_table$strand == "+"
        utr_infos_table$X3_utr_start <- ifelse(plus, utr_infos_table$tx_cds_seq_end + 1, utr_infos_table$start)
        utr_infos_table$X3_utr_end <- ifelse(plus, utr_infos_table$end, utr_infos_table$tx_cds_seq_start - 1)

        # spliced lengths - exons only
        utrs <- ensembldb::threeUTRsByTranscript(ensdb, filter = filt)
        saveRDS(utrs, file.path(directories$annotations, paste0("utr_coords.rds")))
        cds  <- ensembldb::cdsBy(ensdb, by = "tx", filter = filt)

        tx <- utr_infos_table$tx_id
        utr_infos_table$X3_utr_length <- unname(sum(width(utrs))[tx])
        utr_infos_table$cds_length <- unname(sum(width(cds))[tx])

        utr_infos_table$X3_utr_length[is.na(utr_infos_table$X3_utr_length)] <- 0

        annotations <- merge(annotations, utr_infos_table[, c("tx_id", "X3_utr_start", "X3_utr_end", "X3_utr_length", "cds_length")],
            by.x = "ensembl_transcript_id", by.y = "tx_id", all.x = TRUE)
        annotations$ensembl_transcript_id <- NULL
        write.table(annotations, file.path(directories$annotations, "annotations.tsv"), quote = FALSE, sep = "\t", row.names = FALSE)

        cat(" Done.\n")
    }
}

extract_cds_sequences <- function(annotations, ensdb, genome, chromosome_filter) {
    if (as.logical(settings$use_caching) && file.exists(file.path(directories$annotations, "cds_sequences.tsv"))) {
        cat(paste("Loaded cds sequences from cache.\n"))
    } else {
        cat(paste("Extracting cds sequences with EnsemblDB..."))

        cds <- ensembldb::cdsBy(ensdb, by = "tx", filter = ~ tx_biotype == "protein_coding" & seq_name == chromosome_filter)
        cds_seqs <- GenomicFeatures::extractTranscriptSeqs(genome, cds)

        cds_table <- data.frame("cds" = as.character(cds_seqs), "ensembl_transcript_id" = names(cds_seqs), row.names = NULL)
        cds_table <- cds_table[cds_table$ensembl_transcript_id %in% annotations$ensembl_transcript_id, ]

        cds_table <- merge(cds_table, annotations[, c("ensembl_transcript_id_version", "ensembl_transcript_id")], by = "ensembl_transcript_id", all.y = TRUE)
        cds_table[is.na(cds_table$cds), ]$cds <- "Sequence unavailable"

        cds_table <- cds_table[, c("ensembl_transcript_id_version", "cds")]
        cds_table <- cds_table[order(cds_table$ensembl_transcript_id_version), ]
        write.table(cds_table, file.path(directories$annotations, "cds_sequences.tsv"), quote = FALSE, sep = "\t", row.names = FALSE)

        cat(" Done.\n")
    }
}



# --- ENTRY POINT ---

annotations <- read.table(file.path(directories$annotations, "annotations.tsv"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)

if (as.logical(settings$use_caching) && file.exists(file.path(directories$annotations, "utr_coords.rds")) &&
        file.exists(file.path(directories$annotations, "utr_sequences.tsv")) && file.exists(file.path(directories$annotations, "cds_sequences.tsv"))) {
    cat(paste("Loaded sequences from cache.\n"))
} else {
    suppressMessages(library(ensembldb))
    suppressMessages(library(AnnotationHub))
    suppressMessages(library(BSgenome.Hsapiens.NCBI.GRCh38))

    annotations$ensembl_transcript_id <- stringr::str_extract(annotations$ensembl_transcript_id_version, "[^.]+")

    gtfdb <- ensDbFromGtf(file.path(directories$preload_data, "annotation.gtf"), organism = "Homo_sapiens", version = settings$ensembl_release, genomeVersion = "GRCh38")
    ensdb <- EnsDb(gtfdb)

    genome <- BSgenome.Hsapiens.NCBI.GRCh38::BSgenome.Hsapiens.NCBI.GRCh38

    chromosome_filter <- strsplit(settings$chromosome_filter, ",")[[1]]

    extract_utr_sequences(annotations, ensdb, genome, chromosome_filter)
    extract_sequence_coordinates(annotations, ensdb, genome, chromosome_filter)
    extract_cds_sequences(annotations, ensdb, genome, chromosome_filter)
}
