
setOldClass("xts")


#' @export
MatrixXTS <- setClass("MatrixXTS",
    slots = list(
        data = "xts",
        dims = "numeric",
        dimnames = "list" 
    )
)

#' @export
create_matrix_xts <- function(covMatList) {
    dims <- dim(covMatList[[1]])
    dates <- as.Date(names(covMatList))
    flat_data <- do.call(rbind, lapply(covMatList, c))
    xts_data <- xts(flat_data, order.by = dates)
    
    matrix_dimnames <- dimnames(covMatList[[1]])
    
    new("MatrixXTS", 
        data = xts_data,
        dims = dims,
        dimnames = matrix_dimnames)
}

#' @export
setMethod("[", "MatrixXTS", function(x, i, j, ..., drop = TRUE) {
    # Get raw data using xts indexing
    raw_data <- x@data[i]
    
    # Always return a matrix with proper dimensions and names
    result <- matrix(as.numeric(raw_data), 
                    nrow = x@dims[1], 
                    ncol = x@dims[2])
    dimnames(result) <- x@dimnames
    result
})

#' @export
setMethod("show", "MatrixXTS", function(object) {
    cat("MatrixXTS object with", NROW(object@data), 
        "matrices", object@dims[1], "x", object@dims[2], "\n",
        "Symbols:", paste(object@dimnames[[1]], collapse=", "), "\n",
        "Date range:", format(range(index(object@data))), "\n")
})

#' @export
setMethod("to.period", "MatrixXTS", function(x, period = "months", ...) {
    new_data <- to.period(x@data, period, ...)
    new("MatrixXTS", data = new_data, dims = x@dims, dimnames = x@dimnames)
})

#' @export
setMethod("lag", "MatrixXTS", function(x, k = 1, ...) {
    new_data <- lag(x@data, k, ...)
    new("MatrixXTS", data = new_data, dims = x@dims, dimnames = x@dimnames)
})

#' @export
setMethod("index", "MatrixXTS", function(x) index(x@data))

#' @export
setMethod("start", "MatrixXTS", function(x) start(x@data))

#' @export
setMethod("end", "MatrixXTS", function(x) end(x@data))

#' @export
`%*%.matrix_xts` <- function(x, y) {
    dims <- attr(x, "matrix_dims")
    result <- NextMethod("%*%")
    if (inherits(result, "xts")) {
        attr(result, "matrix_dims") <- dims
        class(result) <- c("matrix_xts", class(result)[-1])
    }
    result
}

#' @export
setMethod("coredata", "MatrixXTS", function(x, ...) {
    raw_data <- coredata(x@data)
    if (NROW(raw_data) == 1) {
        matrix(raw_data, nrow=x@dims[1], ncol=x@dims[2])
    } else {
        raw_data
    }
})

#' @export
setMethod("Ops", "MatrixXTS", function(e1, e2) {
    if (inherits(e2, "MatrixXTS")) {
        e2 <- coredata(e2)
    }
    new_data <- callGeneric(coredata(e1), e2)
    new("MatrixXTS", data = new_data, dims = e1@dims, dimnames = e1@dimnames)
})