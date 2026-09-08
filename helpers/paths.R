# helpers/paths.R
#
# Location-independent paths. Source this at the top of a script instead of
# hard-coding setwd() to a machine-specific directory.
#
#   source(file.path(dirname(sys.frame(1)$ofile), "..", "..", "helpers", "paths.R"))
#
# or simply, from the repository root:
#
#   source("helpers/paths.R")
#
# Provides:
#   script_dir()  the directory of the script being run
#   out_dir()     where figures and tables should be written.
#                 Defaults to the script's own directory; override with FIG_OUTDIR.
#   data_dir()    where large inputs live. Set FIG_DATADIR, since these are not
#                 in the repository. Errors with a clear message if unset and needed.

script_dir <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) return(normalizePath(dirname(f)))
  for (i in seq_len(sys.nframe())) {
    ofile <- sys.frame(i)$ofile
    if (!is.null(ofile)) return(normalizePath(dirname(ofile)))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable())
    return(normalizePath(dirname(rstudioapi::getActiveDocumentContext()$path)))
  getwd()
}

out_dir <- function() {
  d <- Sys.getenv("FIG_OUTDIR", unset = "")
  if (nzchar(d)) { dir.create(d, showWarnings = FALSE, recursive = TRUE); return(d) }
  script_dir()
}

data_dir <- function(what = "the input data") {
  d <- Sys.getenv("FIG_DATADIR", unset = "")
  if (!nzchar(d))
    stop("This script needs ", what, ", which is not included in this repository ",
         "because of its size. Download it (see README, 'Getting the data') and set ",
         "the FIG_DATADIR environment variable to the folder containing it.\n",
         "  export FIG_DATADIR=/path/to/data", call. = FALSE)
  if (!dir.exists(d)) stop("FIG_DATADIR does not exist: ", d, call. = FALSE)
  d
}

# Cross-platform font handling. windowsFonts() exists only on Windows.
FONT_FAMILY <- local({
  if (.Platform$OS.type == "windows") {
    try(grDevices::windowsFonts(Arial = grDevices::windowsFont("Arial")), silent = TRUE)
    "Arial"
  } else if (requireNamespace("systemfonts", quietly = TRUE) &&
             "Arial" %in% systemfonts::system_fonts()$family) {
    "Arial"
  } else ""
})
