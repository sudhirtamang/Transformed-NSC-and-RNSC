# ---- User settings ----

FUNCTIONS_FILE <- "functions.R"

N_REPLICATIONS <- 200 
MODELS <- c("A", "B")
ERROR_TYPES <- 1:4

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

OUTPUT_DIR <- "table1_simulation_output"


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
    !exists("summarize_table1")) {
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


# ---- Create output directory ----

if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}


# ---- Record settings ----

settings <- list(
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

saveRDS(settings, file.path(OUTPUT_DIR, "simulation_settings.rds"))


# ---- Run simulation ----

set.seed(RANDOM_SEED)

cat("\nStarting simulation\n")
cat("Replications per setting:", N_REPLICATIONS, "\n")
cat("Settings:", paste(
  as.vector(outer(MODELS, ERROR_TYPES, paste, sep = ".")),
  collapse = ", "
), "\n\n")

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
  checkpoint_file = file.path(OUTPUT_DIR, "table1_partial_results.rds")
)


# ---- Save raw replication-level results ----

saveRDS(
  simulation_results,
  file = file.path(OUTPUT_DIR, "table1_raw_results.rds")
)


# ---- Summarize and export results ----

table1_summary <- summarize_table1(simulation_results)

write.csv(
  table1_summary,
  file = file.path(OUTPUT_DIR, "table1_summary.csv"),
  row.names = FALSE
)

print(table1_summary)


# ---- Save session information and finish time ----

settings$finished <- Sys.time()
saveRDS(settings, file.path(OUTPUT_DIR, "simulation_settings.rds"))

writeLines(
  capture.output(sessionInfo()),
  con = file.path(OUTPUT_DIR, "sessionInfo.txt")
)

cat("\nSimulation completed.\n")
cat("Output directory:", normalizePath(OUTPUT_DIR), "\n")
cat("Summary:", file.path(OUTPUT_DIR, "table1_summary.csv"), "\n")
cat("Raw results:", file.path(OUTPUT_DIR, "table1_raw_results.rds"), "\n")
