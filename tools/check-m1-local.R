# Local macOS arm64 check script for BayesRTMB.
#
# Run from the package root directory:
#   source("tools/check-m1-local.R")
#
# Or from a terminal:
#   Rscript tools/check-m1-local.R
#
# The script builds a source tarball, installs it into a temporary local
# library, runs the rtmb_fa regression that failed on CRAN M1mac, runs the
# installed help example, and then runs R CMD check on the built tarball.

options(warn = 1)

pkg_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
desc_path <- file.path(pkg_root, "DESCRIPTION")
if (!file.exists(desc_path)) {
  stop("Run this script from the BayesRTMB package root directory.", call. = FALSE)
}

desc <- read.dcf(desc_path)
pkg_name <- unname(desc[, "Package"])
pkg_version <- unname(desc[, "Version"])
if (!identical(pkg_name, "BayesRTMB")) {
  stop("This does not look like the BayesRTMB package root.", call. = FALSE)
}

timestamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
out_dir <- file.path(dirname(pkg_root), paste0("BayesRTMB-m1-local-check-", timestamp))
lib_dir <- file.path(out_dir, "Rlib")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(lib_dir, recursive = TRUE, showWarnings = FALSE)

message("BayesRTMB local M1 check")
message("Package root: ", pkg_root)
message("Version: ", pkg_version)
message("Output dir: ", out_dir)

writeLines(c(
  "=== Sys.info() ===",
  capture.output(print(Sys.info())),
  "",
  "=== R.version ===",
  capture.output(print(R.version.string)),
  capture.output(print(R.version$platform)),
  capture.output(print(R.version$arch)),
  "",
  "=== macOS ===",
  capture.output(system2("uname", "-a", stdout = TRUE, stderr = TRUE)),
  capture.output(system2("sw_vers", stdout = TRUE, stderr = TRUE)),
  capture.output(system2("sysctl", c("-n", "machdep.cpu.brand_string"), stdout = TRUE, stderr = TRUE))
), file.path(out_dir, "system-info.txt"))

if (!identical(Sys.info()[["sysname"]], "Darwin")) {
  warning("This is not macOS. The script will run, but it is not an M1mac check.")
}
if (!grepl("arm64|aarch64", paste(R.version$arch, R.version$platform), ignore.case = TRUE)) {
  warning("This R session does not appear to be arm64/aarch64. It may be running under Rosetta.")
}

old_libs <- .libPaths()
.libPaths(c(lib_dir, old_libs))
on.exit(.libPaths(old_libs), add = TRUE)

repos <- getOption("repos")
if (is.null(repos) || identical(unname(repos["CRAN"]), "@CRAN@")) {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}

install_missing_deps <- TRUE

parse_dep_names <- function(x) {
  if (length(x) == 0 || all(is.na(x))) return(character())
  x <- paste(x[!is.na(x)], collapse = ",")
  x <- gsub("\\s*\\([^)]*\\)", "", x)
  x <- trimws(unlist(strsplit(x, ",")))
  x[nzchar(x) & x != "R"]
}

if (isTRUE(install_missing_deps)) {
  dep_fields <- intersect(c("Depends", "Imports", "LinkingTo", "Suggests"), colnames(desc))
  direct_deps <- parse_dep_names(desc[, dep_fields, drop = TRUE])
  direct_deps <- setdiff(direct_deps, c("BayesRTMB", "R"))

  if (length(direct_deps) > 0) {
    db <- available.packages()
    recursive <- tools::package_dependencies(
      direct_deps,
      db = db,
      which = c("Depends", "Imports", "LinkingTo"),
      recursive = TRUE
    )
    deps <- unique(c(direct_deps, unlist(recursive, use.names = FALSE)))
    deps <- setdiff(deps, c("R", rownames(installed.packages())))

    if (length(deps) > 0) {
      message("Installing missing dependencies into: ", lib_dir)
      message(paste(deps, collapse = ", "))
      install.packages(deps, lib = lib_dir, dependencies = FALSE)
    } else {
      message("No missing dependencies detected.")
    }
  }
}

r_bin <- file.path(R.home("bin"), "R")
rscript_bin <- file.path(R.home("bin"), "Rscript")
r_libs_user <- paste(c(lib_dir, Sys.getenv("R_LIBS_USER")), collapse = .Platform$path.sep)
cmd_env <- c(
  paste0("R_LIBS_USER=", r_libs_user),
  "R_KEEP_PKG_SOURCE=yes",
  "_R_CHECK_FORCE_SUGGESTS_=false",
  "_R_CHECK_INSTALL_DEPENDS_=true",
  "_R_CHECK_DEPRECATED_DEFUNCT_=true",
  "_R_CHECK_CODE_ASSIGN_TO_GLOBALENV_=true",
  "_R_CHECK_CODE_DATA_INTO_GLOBALENV_=true",
  "_R_CHECK_SCREEN_DEVICE_=warn",
  "_R_CHECK_S3_METHODS_NOT_REGISTERED_=true",
  "_R_CHECK_OVERWRITE_REGISTERED_S3_METHODS_=true",
  "_R_CHECK_NATIVE_ROUTINE_REGISTRATION_=true",
  "_R_CHECK_FF_CALLS_=registration",
  "_R_CHECK_COMPILATION_FLAGS_=true",
  "_R_CHECK_THINGS_IN_TEMP_DIR_=true",
  "_R_CHECK_MATRIX_DATA_=TRUE",
  "_R_CHECK_ORPHANED_=true",
  "_R_CHECK_BROWSER_NONINTERACTIVE_=true",
  "_R_CHECK_MBCS_CONVERSION_FAILURE_=true",
  "_R_CHECK_LIMIT_CORES_=true",
  "_R_CHECK_TIMINGS_=10",
  "_R_CHECK_TESTS_NLINES_=0",
  "_R_CHECK_VIGNETTES_NLINES_=0",
  "R_BROWSER=false",
  "R_PDFVIEWER=false",
  "RGL_USE_NULL=true",
  "NOAWT=1"
)

run_cmd <- function(label, command, args, wd = out_dir, env = cmd_env) {
  log_file <- file.path(out_dir, paste0(label, ".log"))
  message("Running ", label, " ...")
  old_wd <- setwd(wd)
  on.exit(setwd(old_wd), add = TRUE)

  status <- system2(
    command,
    args = args,
    stdout = log_file,
    stderr = log_file,
    env = env
  )

  if (!identical(status, 0L)) {
    message("FAILED: ", label)
    message("Log file: ", log_file)
    message(paste(tail(readLines(log_file, warn = FALSE), 80), collapse = "\n"))
    stop(label, " failed.", call. = FALSE)
  }

  message("OK: ", label)
  invisible(log_file)
}

run_cmd(
  "01-R-CMD-build",
  r_bin,
  c("CMD", "build", shQuote(pkg_root)),
  wd = out_dir
)

tarballs <- list.files(
  out_dir,
  pattern = paste0("^", pkg_name, "_.*[.]tar[.]gz$"),
  full.names = TRUE
)
if (length(tarballs) != 1L) {
  stop("Expected exactly one source tarball in ", out_dir, call. = FALSE)
}
tarball <- normalizePath(tarballs, winslash = "/", mustWork = TRUE)
message("Built tarball: ", tarball)

tar_files <- utils::untar(tarball, list = TRUE)
writeLines(tar_files, file.path(out_dir, "tarball-files.txt"))
inst_doc <- grep(paste0("^", pkg_name, "/inst/doc/"), tar_files, value = TRUE)
writeLines(inst_doc, file.path(out_dir, "inst-doc-files.txt"))
if (length(inst_doc) == 0L) {
  stop("The built tarball does not contain inst/doc files.", call. = FALSE)
}

run_cmd(
  "02-R-CMD-INSTALL",
  r_bin,
  c("CMD", "INSTALL", paste0("--library=", shQuote(lib_dir)), shQuote(tarball)),
  wd = out_dir
)

regression_script <- file.path(out_dir, "03-installed-regression.R")
writeLines(c(
  sprintf(".libPaths(c(%s, .libPaths()))", deparse(lib_dir)),
  "library(BayesRTMB)",
  "cat('Installed BayesRTMB version:', as.character(packageVersion('BayesRTMB')), '\\n')",
  sprintf("stopifnot(packageVersion('BayesRTMB') >= %s)", deparse(pkg_version)),
  "cat('=== Session info ===\\n')",
  "print(sessionInfo())",
  "cat('=== CRAN M1mac rtmb_fa regression ===\\n')",
  "data(BigFive, package = 'BayesRTMB')",
  "fa_data <- BigFive[, 1:10]",
  "fit_fa2 <- rtmb_fa(data = fa_data, nfactors = 2, score = TRUE)",
  "stopifnot(inherits(fit_fa2, 'RTMB_Model'))",
  "code <- capture.output(fit_fa2$print_code())",
  "stopifnot(!any(grepl('rowSums\\\\(L_raw', code)))",
  "stopifnot(!any(grepl('L_raw\\\\^2', code)))",
  "stopifnot(any(grepl('for \\\\(k in 1:min\\\\(j, K\\\\)\\\\)', code)))",
  "map_fa2 <- fit_fa2$optimize(num_estimate = 1, se_method = 'none')",
  "print(map_fa2)",
  "print(map_fa2$summary(max_rows = 10))",
  "cat('=== Post-hoc Promax rotation ===\\n')",
  "map_fa2$fa_rotate(target = 'L', scores = 'score', rotate = 'promax')",
  "print(map_fa2$summary('L_promax', max_rows = 10))",
  "cat('=== Full installed help example ===\\n')",
  "example(rtmb_fa, package = 'BayesRTMB', ask = FALSE)",
  "cat('installed regression OK\\n')"
), regression_script)

run_cmd(
  "03-installed-regression",
  rscript_bin,
  shQuote(regression_script),
  wd = out_dir
)

run_cmd(
  "04-R-CMD-check",
  r_bin,
  c("CMD", "check", "--as-cran", "--no-manual", shQuote(tarball)),
  wd = out_dir
)

summary_file <- file.path(out_dir, "SUMMARY.txt")
writeLines(c(
  "BayesRTMB local M1 check: OK",
  paste("Package:", pkg_name),
  paste("Version:", pkg_version),
  paste("Tarball:", tarball),
  paste("Output dir:", out_dir),
  "",
  "Completed checks:",
  "* R CMD build",
  "* tarball contains inst/doc",
  "* R CMD INSTALL into temporary local library",
  "* installed rtmb_fa nfactors = 2, score = TRUE regression",
  "* post-hoc Promax rotation",
  "* installed example(rtmb_fa)",
  "* R CMD check --as-cran --no-manual"
), summary_file)

message("")
message("All checks completed successfully.")
message("Summary: ", summary_file)
message("Please send the entire output directory back if possible:")
message(out_dir)
