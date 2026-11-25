#' Read an external vector library
#'
#' Read selected rows from a little-endian, row-major float32 vector library.
#' The library is stored outside the package. Row numbers are one-based and
#' output rows retain the requested order, including duplicates.
#'
#' @param path Directory containing `meta.json`, `ids.tsv`, and `embeddings.f32`.
#' @param rows Integer row indices; `NULL` reads every row.
#' @return A list with numeric matrix `vectors`, corresponding `ids`, and
#'   the metadata. Float32 values are represented exactly as R doubles.
#' @examples
#' \dontrun{lib <- read_vector_library("vector-library", 1:3)}
#' @seealso [readBin()], [jsonlite::fromJSON()]
#' @export
read_vector_library <- function(path, rows = NULL) {
  meta <- jsonlite::fromJSON(file.path(path, "meta.json"), simplifyVector = TRUE)
  shape <- as.double(meta$shape)
  if (length(shape) != 2L || any(!is.finite(shape)) || any(shape < 1) ||
      any(shape != floor(shape)) || !identical(meta$dtype, "float32") ||
      !identical(meta$byte_order, "little") || !identical(meta$order, "C"))
    stop("Invalid or unsupported vector metadata", call. = FALSE)
  bin <- file.path(path, "embeddings.f32")
  if (is.na(file.info(bin)$size) || file.info(bin)$size != prod(shape) * 4)
    stop("Vector file size does not match metadata", call. = FALSE)
  ids <- utils::read.delim(file.path(path, "ids.tsv"), check.names = FALSE,
                           stringsAsFactors = FALSE, quote = "", comment.char = "")
  if (nrow(ids) != shape[1]) stop("ID count does not match vector rows", call. = FALSE)
  if (is.null(rows)) rows <- seq_len(shape[1])
  if (!is.numeric(rows) || anyNA(rows) || any(!is.finite(rows)) ||
      any(rows != floor(rows)) || any(rows < 1 | rows > shape[1]))
    stop("rows must be valid one-based indices", call. = FALSE)
  con <- file(bin, "rb"); on.exit(close(con))
  result <- matrix(numeric(length(rows) * shape[2]), nrow = length(rows), ncol = shape[2])
  for (i in seq_along(rows)) {
    seek(con, where = (rows[i] - 1) * shape[2] * 4, origin = "start")
    v <- readBin(con, what = "numeric", n = shape[2], size = 4, endian = "little")
    if (length(v) != shape[2]) stop("Truncated vector file", call. = FALSE)
    result[i, ] <- v
  }
  list(vectors = result, ids = tibble::as_tibble(ids[rows, , drop = FALSE]), meta = meta)
}
