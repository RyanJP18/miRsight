args <- commandArgs(trailingOnly = TRUE)
suppressMessages(library(jsonlite))
config <- jsonlite::fromJSON(args[1])
settings <- config$settings
directories <- config$directories



# --- FUNCTIONS ---

generate_conservation_scores <- function(track, filename, transcripts) {
    conservation_name <- gsub("*.tsv", "", filename)
    if (as.logical(settings$use_caching) && file.exists(file.path(directories$conservation, filename))) {
        cat(paste("Loaded", conservation_name, "from cache.\n"))
    } else {
        cat(paste("Generating", conservation_name, "scores..."))

        suppressMessages(library(GenomicScores))
        suppressMessages(library(phyloP100way.UCSC.hg38))

        sink(file.path(directories$conservation, filename))
        for (i in seq_len(nrow(transcripts))) {
            if (is.na(transcripts$X3_utr_start[i]) || is.na(transcripts$X3_utr_end[i])) {
                cat(transcripts$ensembl_transcript_id[i])
                cat(" ")
                cat(NA)
                cat("\n")
                next
            }

            cat(transcripts$ensembl_transcript_id[i])
            cat(" ")
            cat(GenomicScores::gscores(track, GRanges(seqnames = transcripts$chrom[i], IRanges(start = transcripts$X3_utr_start[i]:transcripts$X3_utr_end[i], width = 1)))$default)
            cat("\n")
        }
        sink()
        cat(" Done.\n")
    }
}



# --- ENTRY POINT ---

annotations <- read.table(file.path(directories$annotations, "annotations.tsv"), sep = "\t", header = TRUE)
generate_conservation_scores(phyloP100way.UCSC.hg38, "phylo100.tsv", annotations)