# helpers/paths.R
#
# Location-independent paths. Source this near the top of a script instead of
# hard-coding setwd() to a machine-specific directory.
#
# Paths are resolved ONCE, at source time, and stored as absolute paths. A later
# setwd() therefore cannot break them, which matters because most scripts call
# setwd(out_dir()) and then read their inputs.
#
# Provides:
#   script_dir()  directory of the script being run
#   repo_root()   top of this repository
#   out_dir()     where to write. Defaults to the script's directory; FIG_OUTDIR overrides.
#   data_dir(f)   where to read a large input from. Looks in the bundled data/ folder
#                 first, then FIG_DATADIR. Errors with an actionable message if absent.
#   FONT_FAMILY   "Arial" where available, "" otherwise. windowsFonts() is Windows-only.

.PATHS <- local({
  sd <- {
    a <- commandArgs(trailingOnly = FALSE)
    f <- sub("^--file=", "", a[grep("^--file=", a)])
    d <- NULL
    if (length(f) && nzchar(f[1])) d <- dirname(f[1])
    if (is.null(d)) for (i in seq_len(sys.nframe())) {
      of <- sys.frame(i)$ofile
      if (!is.null(of)) { d <- dirname(of); break }
    }
    if (is.null(d) && requireNamespace("rstudioapi", quietly = TRUE) &&
        rstudioapi::isAvailable())
      d <- dirname(rstudioapi::getActiveDocumentContext()$path)
    if (is.null(d) || !nzchar(d)) d <- getwd()
    normalizePath(d, mustWork = FALSE)      # absolute, before any setwd()
  }
  rr <- sd
  while (!file.exists(file.path(rr, "helpers", "paths.R")) && dirname(rr) != rr)
    rr <- dirname(rr)
  list(script = sd, root = rr)
})

script_dir <- function() .PATHS$script
repo_root  <- function() .PATHS$root

out_dir <- function() {
  d <- Sys.getenv("FIG_OUTDIR", unset = "")
  if (nzchar(d)) {
    dir.create(d, showWarnings = FALSE, recursive = TRUE)
    return(normalizePath(d, mustWork = FALSE))
  }
  script_dir()
}

data_dir <- function(what = "the input data") {
  bundled <- file.path(repo_root(), "data")
  if (dir.exists(bundled) && file.exists(file.path(bundled, what))) return(bundled)
  d <- Sys.getenv("FIG_DATADIR", unset = "")
  if (nzchar(d) && file.exists(file.path(d, what))) return(normalizePath(d, mustWork = FALSE))
  if (nzchar(d) && dir.exists(d)) return(normalizePath(d, mustWork = FALSE))
  stop("This script needs ", what, ", which is not bundled in this repository ",
       "because of its size. See README, 'Getting the data', then point FIG_DATADIR ",
       "at the folder containing it:\n  export FIG_DATADIR=/path/to/data", call. = FALSE)
}

FONT_FAMILY <- local({
  if (.Platform$OS.type == "windows") {
    try(grDevices::windowsFonts(Arial = grDevices::windowsFont("Arial")), silent = TRUE)
    "Arial"
  } else if (requireNamespace("systemfonts", quietly = TRUE) &&
             "Arial" %in% systemfonts::system_fonts()$family) {
    "Arial"
  } else ""
})

# Report missing packages clearly rather than installing into someone else's library.
require_pkgs <- function(...) {
  pkgs <- unlist(list(...))
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing))
    stop("Missing R packages: ", paste(missing, collapse = ", "),
         "\nInstall them, then re-run. Versions used for the paper are in ",
         "helpers/package_versions.csv.", call. = FALSE)
  invisible(lapply(pkgs, library, character.only = TRUE))
}

# ---- source data export (Science "data S1") ----------------------------
# Writes the data behind a saved figure as CSV, next to the figure. Science
# requires the tabulated values underlying every panel, including P values.
#
#   p <- ggplot(...) ; ggsave("Fig1B.pdf", p) ; save_source_data(p, "Fig1B.pdf")
#
# Set FIG_NO_SOURCE_DATA=1 to skip.
save_source_data <- function(plot_obj, pdf_path) {
  if (nzchar(Sys.getenv("FIG_NO_SOURCE_DATA", unset = ""))) return(invisible(NULL))

  clean <- function(d) {
    if (!is.data.frame(d) || !nrow(d)) return(NULL)
    d[, !vapply(d, is.list, logical(1)), drop = FALSE]
  }

  # 1. the plot-level data, where there is any
  d <- clean(tryCatch(ggplot2::ggplot_build(plot_obj)$plot$data, error = function(e) NULL))
  if (is.null(d)) d <- clean(tryCatch(plot_obj$data, error = function(e) NULL))

  # 2. otherwise the per-layer data. ggplot() is often called empty, with each
  #    geom supplying its own data, in which case the plot-level slot is a waiver.
  if (is.null(d)) {
    layers <- tryCatch(plot_obj$layers, error = function(e) NULL)
    parts <- list()
    for (i in seq_along(layers)) {
      ld <- clean(tryCatch(layers[[i]]$data, error = function(e) NULL))
      if (!is.null(ld)) {
        ld$.layer <- i
        ld$.geom  <- tryCatch(class(layers[[i]]$geom)[1], error = function(e) NA_character_)
        parts[[length(parts) + 1]] <- ld
      }
    }
    if (length(parts)) {
      cols <- Reduce(intersect, lapply(parts, names))
      d <- if (length(cols) > 2) do.call(rbind, lapply(parts, function(x) x[, cols, drop = FALSE]))
           else parts[[1]]
    }
  }

  if (is.null(d)) {
    message("source data: nothing extractable for ", basename(pdf_path))
    return(invisible(NULL))
  }
  out <- sub("\\.(pdf|png|svg|tiff)$", "_source_data.csv", pdf_path)
  if (identical(out, pdf_path)) out <- paste0(pdf_path, "_source_data.csv")
  utils::write.csv(d, out, row.names = FALSE)
  message("source data: ", basename(out), "  (", nrow(d), " rows x ", ncol(d), " cols)")
  invisible(out)
}
