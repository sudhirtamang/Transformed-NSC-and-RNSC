# ============================================================
# Run the Table 1 simulation at several contamination levels
#
# This script sources functions.R; keep both files in
# the same directory.
#
# Run from a terminal:
#   Rscript simulation.R
#
# Or open this file in RStudio and click "Source".
#
# For every value in CONTAMINATION_LEVELS the full simulation is run
# (same models, error types, methods and replications) and saved in its
# own folder. A combined summary over all levels is written at the end.
# A level whose results already exist is skipped, so an interrupted run
# can be restarted without redoing finished levels.
# ============================================================

# ---- User settings ----

FUNCTIONS_FILE <- "functions.R"

# Percentage (0-100) of training entries replaced by t(3) draws.
# 0 = no contamination (the original simulation).
CONTAMINATION_LEVELS <- c(0, 5, 10, 15)

N_REPLICATIONS <- 200
MODELS <- c("A", "B")
ERROR_TYPES <- 1:2

PREDICTORS <- 2000
TRAIN_PER_CLASS <- 100

# The paper's Section 4 does not specify the test-set size.
# This is an implementation choice; change it if you have
# additional information from the authors' code or supplement.
TEST_PER_CLASS <- 1000

RANDOM_SEED <- 123
N_FOLDS <- 10

# Number of parallel workers (1 = no parallelism). Leave one core free.
N_CORES <- max(1L, parallel::detectCores(logical = FALSE) - 1L, na.rm = TRUE)

# Illustrative tuning grids: the paper does not report exact grids.
RNSC_H_GRID <- c(1, 2, 3, 5, 10, Inf)
RNSC_HUBER_GRID <- c(1, 2, 3, 5)
RNSC_LAMBDA_GRID <- c(0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 6, 8)

OUTPUT_ROOT <- "contamination_simulation_output"


# ---- Locate and load the functions file ----

# If running via Rscript, first try the directory containing this script.
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)

if (length(file_arg) > 0) {
  script_path <- normalizePath(
    sub("^--file=", "", file_arg[1]),
    winslash = "/",
    mustWork = FALSE
  )
  script_dir <- dirname(script_path)

  candidate <- file.path(script_dir, FUNCTIONS_FILE)

  if (file.exists(candidate)) {
    FUNCTIONS_FILE <- candidate
  }
}

if (!file.exists(FUNCTIONS_FILE)) {
  stop(
    "Cannot find ", FUNCTIONS_FILE, ". Put it in the same directory ",
    "as this script, or update FUNCTIONS_FILE at the top."
  )
}

source(FUNCTIONS_FILE)

if (!exists("run_table1_simulation") ||
    !exists("summarize_table1") ||
    !exists("run_one_replication_contaminated")) {
  stop("The functions file did not define the expected simulation functions.")
}


# ---- Check required packages ----

required_packages <- c("pamr", "glmnet")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace,
          quietly = TRUE, FUN.VALUE = logical(1))
]

if (length(missing_packages) > 0) {
  stop(
    "Install the following packages before running: ",
    paste(missing_packages, collapse = ", "),
    "\nUse: install.packages(c(",
    paste(sprintf('"%s"', missing_packages), collapse = ", "),
    "))"
  )
}


# ---- Check contamination levels and create the output root ----

if (any(!is.finite(CONTAMINATION_LEVELS)) ||
    any(CONTAMINATION_LEVELS < 0) || any(CONTAMINATION_LEVELS > 100)) {
  stop("CONTAMINATION_LEVELS must be percentages between 0 and 100.")
}

if (!dir.exists(OUTPUT_ROOT)) {
  dir.create(OUTPUT_ROOT, recursive = TRUE)
}


# ---- Run the simulation at every contamination level ----

cat("\nStarting simulation\n")
cat("Contamination levels (%):", paste(CONTAMINATION_LEVELS, collapse = ", "), "\n")
cat("Replications per setting:", N_REPLICATIONS, "\n")
cat("Settings:", paste(
  as.vector(outer(MODELS, ERROR_TYPES, paste, sep = ".")),
  collapse = ", "
), "\n")

all_summaries <- list()

for (CONTAM in CONTAMINATION_LEVELS) {

  level_dir <- file.path(OUTPUT_ROOT, paste0("contam_", CONTAM, "pct"))
  if (!dir.exists(level_dir)) dir.create(level_dir, recursive = TRUE)

  raw_file <- file.path(level_dir, "table1_raw_results.rds")

  cat("\n==================================================\n")
  cat("Contamination level:", CONTAM, "%\n")
  cat("==================================================\n")

  if (file.exists(raw_file)) {

    cat("Results already exist; skipping the run and reusing them.\n")
    simulation_results <- readRDS(raw_file)

  } else {

    settings <- list(
      contamination = CONTAM,
      replications = N_REPLICATIONS,
      models = MODELS,
      error_types = ERROR_TYPES,
      predictors = PREDICTORS,
      training_observations_per_class = TRAIN_PER_CLASS,
      test_observations_per_class = TEST_PER_CLASS,
      seed = RANDOM_SEED,
      folds = N_FOLDS,
      cores = N_CORES,
      rnsc_H_grid = RNSC_H_GRID,
      rnsc_h_grid = RNSC_HUBER_GRID,
      rnsc_lambda_grid = RNSC_LAMBDA_GRID,
      functions_file = FUNCTIONS_FILE,
      started = Sys.time()
    )

    saveRDS(settings, file.path(level_dir, "simulation_settings.rds"))

    simulation_results <- run_table1_simulation(
      B = N_REPLICATIONS,
      models = MODELS,
      error_types = ERROR_TYPES,
      p = PREDICTORS,
      n_per_class = TRAIN_PER_CLASS,
      test_n_per_class = TEST_PER_CLASS,
      seed = RANDOM_SEED,
      rnsc_H_grid = RNSC_H_GRID,
      rnsc_h_grid = RNSC_HUBER_GRID,
      rnsc_lambda_grid = RNSC_LAMBDA_GRID,
      nfold = N_FOLDS,
      progress = TRUE,
      n_cores = N_CORES,
      checkpoint_file = file.path(level_dir, "table1_partial_results.rds"),
      contamination = CONTAM
    )

    saveRDS(simulation_results, file = raw_file)

    settings$finished <- Sys.time()
    saveRDS(settings, file.path(level_dir, "simulation_settings.rds"))

    writeLines(
      capture.output(sessionInfo()),
      con = file.path(level_dir, "sessionInfo.txt")
    )
  }

  # Summarize and export this level
  table1_summary <- summarize_table1(simulation_results)

  write.csv(
    table1_summary,
    file = file.path(level_dir, "table1_summary.csv"),
    row.names = FALSE
  )

  print(table1_summary)

  table1_summary$Contamination <- CONTAM
  all_summaries[[as.character(CONTAM)]] <- table1_summary
}


# ---- Combined summary over all contamination levels ----

combined <- do.call(rbind, all_summaries)
combined <- combined[, c("Contamination",
                         setdiff(names(combined), "Contamination"))]
rownames(combined) <- NULL

write.csv(
  combined,
  file = file.path(OUTPUT_ROOT, "contamination_summary_all.csv"),
  row.names = FALSE
)

cat("\nSimulation completed.\n")
cat("Output directory:", normalizePath(OUTPUT_ROOT), "\n")
cat("One folder per level, e.g.:",
    file.path(OUTPUT_ROOT, paste0("contam_", CONTAMINATION_LEVELS[1], "pct")), "\n")
cat("Combined summary:",
    file.path(OUTPUT_ROOT, "contamination_summary_all.csv"), "\n")
