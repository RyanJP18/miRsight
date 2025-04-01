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
            tryCatch({
                downloaded_files[i] <- cache(ah[ids[i]])
                success[i] <- TRUE
            }, error = function(e) { })
        }
        
        if (any(!success)) {
            attempt <- attempt + 1
            message(sprintf("Retrying annotationHub download attempt %d for %s.", attempt, ids[i]))
            Sys.sleep(5) # delay retrying a few seconds
        }
    }
    
    return(downloaded_files[success])
}

load_phylo_pkg <- function() {

    if (!as.logical(settings$use_precompiled_conservation) || !file.exists(file.path(directories$conservation, "phylo100.tsv"))) {
        suppressMessages(library(GenomicScores))
        suppressMessages(library(AnnotationHub))

        message("Downloading phyloP100way.UCSC.hg38 data from AnnotationHub...")
        # This is a reimplementation of the GenomicScores getGScores function https://rdrr.io/bioc/GenomicScores/src/R/getGScores.R
        # which adds a workaround for an issue where there is a baked in concurrent download limit of 10 in Docker/Python/R 
        # ---
        ah <- AnnotationHub()
        ah <- query(ah, 'phyloP100way.UCSC.hg38')

        ## use the AnnotationHub metadata to figure out the correspondence
        ## between downloaded filenames and object names
        mdah <- mcols(ah)
        objnames <- mdah$title
        ahids <- rownames(mdah)

        # Custom miRsight code: fix download limit by getting files in smaller batches
        batch_size <- 10 
        fnames <- character(0) 

        # for (i in seq(1, length(ahids), by = batch_size)) {
        #     subset_ids <- ahids[i:min(i + batch_size - 1, length(ahids))]
        #     fnames <- c(fnames, cache(ah[subset_ids]))
        # }
        for (i in seq(1, length(ahids), by = batch_size)) {
            subset_ids <- ahids[i:min(i + batch_size - 1, length(ahids))]
            files <- ah_download(ah, subset_ids)
            fnames <- c(fnames, files)
        }
        # End of custom miRsight code:

        ## in BioC 3.8 names from 'cache()' became 'AHXXXX : YYYY'
        ## so we need to override those names.
        names(fnames) <- ahids
        serializedobjs <- basename(fnames[ahids])
        names(serializedobjs) <- objnames

        ## load the first object to get the metadata
        obj <- ah[[1]]
        mdobj <- metadata(obj)
        gsco <- GScores(provider=mdobj$provider,
                        provider_version=mdobj$provider_version,
                        download_url=mdobj$download_url,
                        download_date=mdobj$download_date,
                        reference_genome=mdobj$reference_genome,
                        data_pkgname=mdobj$data_pkgname,
                        data_dirpath=getAnnotationHubOption("CACHE"),
                        data_serialized_objnames=serializedobjs)
        scorlelist <- get(mdobj$data_pkgname, envir=gsco@.data_cache)
        scorlelist[[mdobj$seqname]] <- obj
        assign(mdobj$data_pkgname, scorlelist, envir=gsco@.data_cache)
        # ---
        
        message("Creating phyloP100way.UCSC.hg38 package...")
        makeGScoresPackage(gsco, maintainer="miRsight <mirsight@mirsight.info>", author="miRsight", version="1.0.0")
        suppressMessages(install.packages("./phyloP100way.UCSC.hg38", repos=NULL, type="source", lib="/rlib"))
        message("Installed phyloP100way.UCSC.hg38 package.")
    }

    suppressMessages(library(phyloP100way.UCSC.hg38))
}

generate_conservation_scores <- function(track, filename, transcripts) {
    conservation_name <- gsub("*.tsv", "", filename)
    if (as.logical(settings$use_caching) && file.exists(file.path(directories$conservation, filename))) {
        message(paste("Conservation track", conservation_name, "- loaded from cache.\n"))
    } else {
        load_phylo_pkg()

        message(paste("Generating", conservation_name, "scores... This will take a while..."))

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
        
        message(paste("Generating", conservation_name, "scores - done."))
    }
}



# --- ENTRY POINT ---

annotations <- read.table(file.path(directories$annotations, "annotations.tsv"), sep = "\t", header = TRUE)
generate_conservation_scores(phyloP100way.UCSC.hg38, "phylo100.tsv", annotations)