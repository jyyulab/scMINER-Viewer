#!/usr/bin/env Rscript
# Safely (re)install scminerViewer without tripping the
# "lazy-load database is corrupt" / "interrupted promise evaluation"
# pitfall.
#
# Why this script exists
# ----------------------
# R memory-maps a package's R/<pkg>.rdb file the first time `library()`
# resolves a promise from it. R does NOT release that file handle on
# reinstall — so if any R process (including the current one) has
# `library(scminerViewer)` active when you run `R CMD INSTALL`, the new
# .rdb gets written over the old one while the OS still serves the old
# bytes via the existing mmap. The next promise evaluation reads
# garbage and you get
#
#   Error: lazy-load database '.../scminerViewer.rdb' is corrupt
#   Warning: restarting interrupted promise evaluation
#   Error in R_decompress1
#
# Strategy here
# -------------
# 1. Detach + unload the namespace in the current R session.
# 2. gc() to release mmap pages held by this process.
# 3. Wipe the install dir (any partial old write).
# 4. Run `R CMD INSTALL` in a fresh callr child process — the child
#    has NEVER loaded scminerViewer, so it cannot hold a stale mmap.
#
# What this script CANNOT do
# --------------------------
# Unload the package from OTHER R sessions. If you have a second
# RStudio open with `library(scminerViewer)` already done, you must
# quit (not just restart) that R process first.

safe_install <- function(pkg = "scminerViewer",
                         lib = NULL,
                         install_deps = FALSE,
                         verbose = TRUE) {
  pkg_name <- basename(pkg)
  if (is.null(lib)) lib <- .libPaths()[1]
  lib <- normalizePath(lib, mustWork = FALSE)
  dir.create(lib, showWarnings = FALSE, recursive = TRUE)

  msg <- function(...) if (isTRUE(verbose)) message(...)

  # 1. Detach + unload in THIS session, if loaded.
  ns_search <- paste0("package:", pkg_name)
  if (ns_search %in% search()) {
    msg("Detaching package:", pkg_name)
    try(detach(ns_search, unload = TRUE, character.only = TRUE),
        silent = TRUE)
  }
  if (pkg_name %in% loadedNamespaces()) {
    msg("Unloading namespace:", pkg_name)
    try(unloadNamespace(pkg_name), silent = TRUE)
  }

  # 2. Force GC so any mmap'd .rdb pages held by THIS process are
  #    released before we wipe the install dir.
  gc(verbose = FALSE, full = TRUE)

  # 3. Wipe stale install (any half-written .rdb / Meta state).
  installed_at <- file.path(lib, pkg_name)
  if (dir.exists(installed_at)) {
    msg("Removing existing install at ", installed_at)
    unlink(installed_at, recursive = TRUE, force = TRUE)
  }

  # 4. Run R CMD INSTALL in a fresh subprocess. callr handles this
  #    cleanly; if callr isn't installed, fall back to system().
  if (!requireNamespace("callr", quietly = TRUE)) {
    msg("callr not installed; falling back to system R CMD INSTALL")
    cmd <- sprintf('"%s" CMD INSTALL --library="%s" %s',
                   file.path(R.home("bin"), "R"), lib, shQuote(pkg))
    status <- system(cmd)
    if (status != 0L) {
      stop("R CMD INSTALL failed (status ", status, ")")
    }
  } else {
    msg("Installing in callr subprocess...")
    callr::rcmd_safe(
      cmd  = "INSTALL",
      cmdargs = c(paste0("--library=", lib), pkg),
      show = isTRUE(verbose),
      echo = isTRUE(verbose)
    )
  }

  msg("Done. Reload with: library(", pkg_name, ", lib.loc = '", lib, "')")
  invisible(file.path(lib, pkg_name))
}

# CLI entry point: `Rscript install.R [pkg_path] [lib_path]`
if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  pkg  <- if (length(args) >= 1) args[1] else "scminerViewer"
  lib  <- if (length(args) >= 2) args[2] else NULL
  safe_install(pkg = pkg, lib = lib)
}
