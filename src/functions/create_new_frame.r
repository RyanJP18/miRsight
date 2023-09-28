# generic method for creating a new fixed size data frame with named columes in one line
create_new_frame <- function(col_names, row_names = NULL, row_count = 0) {
    if (!is.null(row_names)) {
        row_count <- length(row_names)
    }

    mat <- matrix(nrow = row_count, ncol = length(col_names))
    colnames(mat) <- col_names

    if (is.null(row_names)) {
        rownames(mat) <- row_names
    }

    return(data.frame(mat))
}