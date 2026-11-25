# Reproduce the compact JSON convention used for recorded model input. R's JSON
# encoder rounds doubles and drops the .0 that Python emits for float values.
.msi_python_float <- function(x) {
  if (!is.finite(x)) stop("Nonfinite JSON number")
  if (x == 0) return(if (is.infinite(1 / x) && 1 / x < 0) "-0.0" else "0.0")
  short <- NULL
  for (digits in seq_len(17L)) {
    candidate <- sprintf(paste0("%.", digits, "g"), x)
    if (identical(as.numeric(candidate), x)) { short <- candidate; break }
  }
  if (is.null(short)) stop("Cannot represent JSON number")
  negative <- startsWith(short, "-")
  if (negative) short <- substring(short, 2L)
  parts <- strsplit(tolower(short), "e", fixed = TRUE)[[1L]]
  exp10 <- if (length(parts) == 2L) as.integer(parts[2L]) else 0L
  mantissa <- parts[1L]
  dot <- regexpr(".", mantissa, fixed = TRUE)[1L]
  before <- if (dot < 0L) nchar(mantissa) else dot - 1L
  digits <- gsub(".", "", mantissa, fixed = TRUE)
  leading <- regexpr("[^0]", digits)[1L] - 1L
  digits <- substring(digits, leading + 1L)
  exponent <- before - leading - 1L + exp10
  digits <- sub("0+$", "", digits)
  if (!nzchar(digits)) digits <- "0"
  if (exponent < -4L || exponent >= 16L) {
    value <- paste0(substr(digits, 1L, 1L), if (nchar(digits) > 1L) paste0(".", substring(digits, 2L)) else "",
                    "e", if (exponent < 0L) "-" else "+", sprintf("%02d", abs(exponent)))
  } else {
    point <- exponent + 1L
    value <- if (point <= 0L) paste0("0.", strrep("0", -point), digits)
      else if (point >= nchar(digits)) paste0(digits, strrep("0", point - nchar(digits)), ".0")
      else paste0(substr(digits, 1L, point), ".", substring(digits, point + 1L))
  }
  paste0(if (negative) "-" else "", value)
}
.msi_python_json <- function(x) {
  if (is.null(x)) return("null")
  if (is.list(x)) {
    values <- vapply(x, .msi_python_json, character(1))
    if (is.null(names(x))) return(paste0("[", paste(values, collapse = ","), "]"))
    if (!length(x)) return("{}")
    keys <- vapply(names(x), function(y) as.character(jsonlite::toJSON(y, auto_unbox = TRUE)), character(1))
    return(paste0("{", paste(paste0(keys, ":", values), collapse = ","), "}"))
  }
  if (length(x) != 1L || is.na(x)) stop("Expected scalar JSON value")
  if (is.character(x)) return(as.character(jsonlite::toJSON(x, auto_unbox = TRUE)))
  if (is.logical(x)) return(if (x) "true" else "false")
  if (is.integer(x)) return(as.character(x))
  if (is.double(x)) return(.msi_python_float(x))
  stop("Unsupported JSON type")
}
