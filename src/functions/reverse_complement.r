# get the reverse complement of the given sequence
reverse_complement <- function(sequence) {
    sequence <- toupper(sequence)
    sequence <- gsub("T", "a", sequence)
    sequence <- gsub("U", "a", sequence)
    sequence <- gsub("A", "t", sequence)
    sequence <- gsub("G", "c", sequence)
    sequence <- gsub("C", "g", sequence)
    sequence <- toupper(sequence)
    return(intToUtf8(rev(utf8ToInt(sequence))))
}
