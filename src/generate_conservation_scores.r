#
# Use a GenomicScores conservation track to produce a genome-wide set of base scores that are easy to extract and store.
# In cases where a track library does not exist, AnnotationHub is used to download necessary data and build a track package.
#

args <- commandArgs(trailingOnly = TRUE)
suppressMessages(library(jsonlite))
config <- jsonlite::fromJSON(args[1])
settings <- config$settings
directories <- config$directories



# --- FUNCTIONS ---

ah_download <- function(ah, ids, max_retries = 5) {
    attempt <- 1
    success <- rep(FALSE, length(ids))
    downloaded_files <- character(length(ids))

    while (attempt <= max_retries && any(!success)) {
        for (i in which(!success)) {
            tryCatch(
                {
                    downloaded_files[i] <- cache(ah[ids[i]])
                    success[i] <- TRUE
                },
                error = function(e) { }
            )
        }

        if (any(!success)) {
            attempt <- attempt + 1
            message(sprintf("Retrying annotationHub download attempt %d for %s.", attempt, ids[i]))
            Sys.sleep(5) # delay retrying a few seconds
        }
    }

    return(downloaded_files[success])
}

build_gscores_pkg <- function(resource) {
    suppressMessages(library(GenomicScores))
    suppressMessages(library(AnnotationHub))

    message("Downloading ", resource, " data from AnnotationHub...")

    # This is a reimplementation of the GenomicScores getGScores function https://rdrr.io/bioc/GenomicScores/src/R/getGScores.R
    # which adds a workaround for an issue where there is a baked in concurrent download limit of 10 in Docker/Python/R
    # --- start of getGScores function
    ah <- AnnotationHub()
    ah <- query(ah, resource)

    ## use the AnnotationHub metadata to figure out the correspondence
    ## between downloaded filenames and object names
    mdah <- mcols(ah)
    objnames <- mdah$title
    ahids <- rownames(mdah)

    # start of custom miRsight code: fix download limit by fetching files in smaller batches
    batch_size <- 10
    fnames <- character(0)

    for (i in seq(1, length(ahids), by = batch_size)) {
        subset_ids <- ahids[i:min(i + batch_size - 1, length(ahids))]
        files <- ah_download(ah, subset_ids)
        fnames <- c(fnames, files)
    }
    # end of custom miRsight code:

    ## in BioC 3.8 names from 'cache()' became 'AHXXXX : YYYY'
    ## so we need to override those names.
    names(fnames) <- ahids
    serializedobjs <- basename(fnames[ahids])
    names(serializedobjs) <- objnames

    ## load the first object to get the metadata
    obj <- ah[[1]]
    mdobj <- metadata(obj)
    gsco <- GScores(
        provider = mdobj$provider,
        provider_version = mdobj$provider_version,
        download_url = mdobj$download_url,
        download_date = mdobj$download_date,
        reference_genome = mdobj$reference_genome,
        data_pkgname = mdobj$data_pkgname,
        data_dirpath = getAnnotationHubOption("CACHE"),
        data_serialized_objnames = serializedobjs
    )
    scorlelist <- get(mdobj$data_pkgname, envir = gsco@.data_cache)
    scorlelist[[mdobj$seqname]] <- obj
    assign(mdobj$data_pkgname, scorlelist, envir = gsco@.data_cache)
    # --- end of getGScores function

    message("Creating ", resource, " package...")
    makeGScoresPackage(gsco, maintainer = "miRsight <mirsight@mirsight.info>", author = "miRsight", version = "1.0.0")
    suppressMessages(install.packages(file.path(".", resource), repos = NULL, type = "source", lib = "/rlib"))
    message("Installed ", resource, " package.")
}

generate_conservation_scores <- function(track, filename, transcripts) {
    conservation_name <- gsub("*.tsv", "", filename)
    if (as.logical(settings$use_caching) && file.exists(file.path(directories$conservation, filename))) {
        message(paste("Conservation track", conservation_name, "- loaded from cache."))
    } else {  
        message(paste("Generating", conservation_name, "scores... This will take a while..."))

        if (!as.logical(settings$use_precompiled_conservation) || !file.exists(file.path(directories$conservation, tsv_name))) {
            build_gscores_pkg(resource)
        }
        suppressMessages(library(resource, character.only = TRUE))

        # per-base genomic positions in transcript (5'->3') order, exons only
        u <- unlist(utrs)
        minus <- as.character(strand(u)) == "-"
        bases <- Map(function(s, e, m) if (m) e:s else s:e, start(u), end(u), minus)
        tx_pos <- split(unlist(bases, use.names = FALSE), rep(names(u), lengths(bases)))

        con <- file(file.path(directories$conservation, filename), "w")
        on.exit(close(con))

        chunk_size <- 2000L
        ids <- transcripts$ensembl_transcript_id          
        key <- sub("\\..*$", "", ids)                     
        ok  <- key %in% names(tx_pos)
        ok[ok] <- lengths(tx_pos[key[ok]]) > 0

        for (block in split(seq_len(nrow(transcripts)), ceiling(seq_len(nrow(transcripts)) / chunk_size))) {
            valid <- block[ok[block]]
            if (length(valid)) {
                pos <- tx_pos[key[valid]]
                n <- lengths(pos)
                res <- GenomicScores::gscores(track, GRanges(
                    seqnames = rep(transcripts$chrom[valid], n),
                    ranges   = IRanges(unlist(pos), width = 1L)
                ))
                stopifnot(all(start(res) == unlist(pos)))   # scores returned in input order
                sc <- split(res$default, factor(rep(seq_along(valid), n), levels = seq_along(valid)))
            }

            j <- 0L
            for (i in block) {
                cat(ids[i], file = con)
                cat(" ", file = con)
                if (ok[i]) { j <- j + 1L; cat(sc[[j]], file = con) } else cat(NA, file = con)
                cat("\n", file = con)
            }
            flush(con)
        }

        message(paste("Generating", conservation_name, "scores - done."))
    }
}



# --- ENTRY POINT ---

annotations <- read.table(file.path(directories$annotations, "annotations.tsv"), sep = "\t", header = TRUE)
utrs <- readRDS(file.path(directories$annotations, "utr-coords.rds"))
generate_conservation_scores(phastCons30way.UCSC.hg38, "phast30.tsv", annotations)
generate_conservation_scores(phastCons100way.UCSC.hg38, "phast100.tsv", annotations)
generate_conservation_scores(phyloP100way.UCSC.hg38, "phylo100.tsv", annotations)
