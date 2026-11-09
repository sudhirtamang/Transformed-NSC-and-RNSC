# ============================================================
# Functions for reproducing selected methods from Table 1
# Ren & Mai (2022): TF-RNSC, RNSC, NSC, Med, and L1-penalized Logistic
#
# IMPORTANT:
# The paper specifies p = 2000, n = 100K, 200 replications,
# the class means, and four error distributions. It does not
# fully specify the test-set size, random seeds, or tuning grids.
# Defaults below for these unspecified choices are illustrative.
#
# Required packages:
# install.packages(c("pamr", "glmnet"))
# ============================================================


# -------------------------
# 1. Data-generation function
# -------------------------

generate_data <- function(model = c("A", "B"),
                          error_type = 1,
                          p = 2000,
                          n_per_class = 100) {

  model <- match.arg(model)

  if (length(error_type) != 1L ||
      !error_type %in% 1:4) {
    stop("error_type must be 1, 2, 3, or 4.")
  }

  if (p < 15) stop("p must be at least 15.")

  if (model == "A") {
    K <- 2L
    mu <- matrix(0, nrow = p, ncol = K)
    mu[1:10, 2] <- 0.5 * c(
      -1.9, -1.8, -1.7, -1.6, -1.5,
       2.5,  2.6,  2.7,  2.8,  2.9
    )
  } else {
    K <- 4L
    mu <- matrix(0, nrow = p, ncol = K)
    mu[1:5, 2] <- c(1.9, 1.3, -1.6, 1.4, -0.7)
    mu[6:10, 3] <- c(-2, 1.4, -1.9, -1.5, 0.8)
    mu[11:15, 4] <- c(-1.8, -2, -1.7, -2, -0.7)
  }

  y <- factor(rep(seq_len(K), each = n_per_class))
  X <- matrix(0, nrow = length(y), ncol = p)

  for (k in seq_len(K)) {
    idx <- which(y == levels(y)[k])
    m <- length(idx)

    if (error_type == 1) {
      # Model X.1: standard normal errors
      E <- matrix(rnorm(m * p), nrow = m, ncol = p)

    } else if (error_type == 2) {
      # Model X.2: t distribution with 3 degrees of freedom
      E <- matrix(rt(m * p, df = 3), nrow = m, ncol = p)

    } else if (error_type == 3) {
      # Model X.3: 0.99 N(0,1) + 0.01 N(10,1)
      E <- matrix(rnorm(m * p), nrow = m, ncol = p)
      outlier <- matrix(rbinom(m * p, 1, 0.01),
                        nrow = m, ncol = p)
      E <- E + 10 * outlier

    } else {
      # Model X.4: t(4) - 9 * Binomial(0.01)
      E <- matrix(rt(m * p, df = 4), nrow = m, ncol = p)
      Z <- matrix(rbinom(m * p, 1, 0.01),
                  nrow = m, ncol = p)
      E <- E - 9 * Z
    }

    X[idx, ] <- sweep(E, 2, mu[, k], FUN = "+")
  }

  colnames(X) <- paste0("X", seq_len(p))

  true_variables <- which(
    apply(mu, 1, function(z) length(unique(z)) > 1)
  )

  list(
    X = X,
    y = y,
    mu = mu,
    true_variables = true_variables
  )
}


# -------------------------
# 2. Huber loss and location
# -------------------------

huber_loss <- function(x, r) {
  if (is.infinite(r)) return(x^2)
  if (length(r) != 1L || is.na(r) || r <= 0) {
    stop("r must be a positive number or Inf.")
  }

  ifelse(abs(x) <= r, x^2, 2 * r * abs(x) - r^2)
}


huber_location <- function(x, H) {
  x <- x[is.finite(x)]

  if (length(x) == 0L) return(NA_real_)
  if (is.infinite(H)) return(mean(x))
  if (length(H) != 1L || is.na(H) || H <= 0) {
    stop("H must be a positive number or Inf.")
  }

  sx <- sd(x)
  if (!is.finite(sx) || sx == 0) return(stats::median(x))

  objective <- function(a) sum(huber_loss(x - a, H))

  lower <- min(x) - 5 * sx
  upper <- max(x) + 5 * sx

  optimize(objective, interval = c(lower, upper))$minimum
}


# Vectorized Huber location estimate for every column of X (IRLS).
# H may be one number (same for all columns) or one value per column.
# Uses the paper's convention: psi is clipped at H (half the derivative of f_H).
huber_location_mat <- function(X, H, tol = 1e-8, maxit = 200, start = NULL) {
  X <- as.matrix(X)
  nr <- nrow(X)

  if (length(H) == 1L && is.infinite(H)) return(colMeans(X))

  Hmat <- if (length(H) == 1L) H else rep(H, each = nr)

  mu <- if (is.null(start)) apply(X, 2, median) else start

  for (it in seq_len(maxit)) {
    R <- sweep(X, 2, mu)
    W <- pmin(Hmat / pmax(abs(R), 1e-12), 1)
    mu_new <- colSums(W * X) / colSums(W)

    if (max(abs(mu_new - mu)) < tol) {
      mu <- mu_new
      break
    }
    mu <- mu_new
  }
  mu
}


# Tuning-free Huber estimation (Wang, Zheng, Zhou and Zhou, 2021).
# For each column, solve jointly
#   sum_i psi_r(Z_i - a) = 0   and   (1/m) sum_i min{(Z_i - a)^2 / r^2, 1} = z/m,
# with z = log(m). The second equation is solved for r by bisection on
# log(r^2); the first is the Huber location estimate for that r.
# (The text of Ren & Mai's eq. (14) lost the r^2 denominator in extraction;
#  this scale-equivariant form is the one used by Wang et al.)
tf_scale <- function(E2, z, iters = 40) {
  m <- nrow(E2)
  hi2 <- pmax(colSums(E2) / z, apply(E2, 2, max), 1e-12)
  lo2 <- hi2 * 1e-10

  for (i in seq_len(iters)) {
    mid2 <- sqrt(lo2 * hi2)
    g <- colSums(pmin(E2 / rep(mid2, each = m), 1))
    too_small <- g > z
    lo2 <- ifelse(too_small, mid2, lo2)
    hi2 <- ifelse(too_small, hi2, mid2)
  }

  (lo2 * hi2)^(1 / 4)
}


huber_tf_mat <- function(X, z = log(nrow(X)), maxit = 30, tol = 1e-6) {
  X <- as.matrix(X)
  alpha <- apply(X, 2, median)

  for (it in seq_len(maxit)) {
    r <- tf_scale(sweep(X, 2, alpha)^2, z)
    alpha_new <- huber_location_mat(X, r, start = alpha)
    change <- max(abs(alpha_new - alpha) / r)
    alpha <- alpha_new
    if (change < tol) break
  }

  alpha
}


# -------------------------
# 3. RNSC fitting and prediction
# -------------------------
# The expensive part (Huber centroids and second moments) depends only on
# the data and H, not on h or lambda. It is computed once by rnsc_stats()
# and reused by rnsc_fit_from_stats() for every (h, lambda).

rnsc_stats <- function(X, y, H, type = c("huber", "tf")) {

  type <- match.arg(type)
  X <- as.matrix(X)
  y <- factor(y)

  classes <- levels(y)
  K <- length(classes)
  n <- nrow(X)
  p <- ncol(X)

  nk <- as.numeric(table(y))
  prior <- nk / n

  centroids <- matrix(0, nrow = p, ncol = K)
  second_moments <- matrix(0, nrow = p, ncol = K)

  for (k in seq_len(K)) {
    Xk <- X[y == classes[k], , drop = FALSE]

    if (type == "huber") {
      centroids[, k] <- huber_location_mat(Xk, H)
      second_moments[, k] <- huber_location_mat(Xk^2, H)
    } else {
      centroids[, k] <- huber_tf_mat(Xk)
      second_moments[, k] <- huber_tf_mat(Xk^2)
    }
  }

  list(
    centroids = centroids,
    second_moments = second_moments,
    nk = nk,
    prior = prior,
    classes = classes,
    n = n,
    p = p
  )
}


rnsc_fit_from_stats <- function(st, h, lambda, H = NA_real_) {

  if (lambda < 0) stop("lambda must be non-negative.")

  K <- length(st$classes)
  p <- st$p
  n <- st$n

  grand_mean <- as.vector(st$centroids %*% st$prior)

  pooled_variance <- rowSums(
    sweep(st$second_moments, 2, st$prior, FUN = "*") -
      sweep(st$centroids^2, 2, st$prior, FUN = "*")
  )

  # Guard against small negative values from numerical error.
  pooled_variance[!is.finite(pooled_variance) |
                    pooled_variance <= 0] <- NA_real_
  scale <- sqrt(pooled_variance)

  shrunken <- matrix(0, nrow = p, ncol = K)

  for (k in seq_len(K)) {
    mk <- sqrt(1 / st$nk[k] - 1 / n)

    d <- (st$centroids[, k] - grand_mean) / (mk * scale)
    d[!is.finite(d)] <- 0

    d_shrunk <- sign(d) * pmax(abs(d) - lambda, 0)
    shrunken[, k] <- grand_mean + mk * scale * d_shrunk
  }

  # A feature is selected if at least one shrunken class centroid
  # differs from the overall centroid.
  deviation <- sweep(shrunken, 1, grand_mean, FUN = "-")
  selected <- which(apply(abs(deviation) > 1e-8, 1, any))

  list(
    centroids = shrunken,
    grand_mean = grand_mean,
    scale = scale,
    prior = st$prior,
    classes = st$classes,
    selected = selected,
    H = H,
    h = h,
    lambda = lambda
  )
}


rnsc_fit <- function(X, y, H, h, lambda) {

  X <- as.matrix(X)
  y <- factor(y)

  if (nrow(X) != length(y)) stop("X and y have incompatible dimensions.")
  if (any(!is.finite(X))) stop("X contains missing or non-finite values.")

  st <- rnsc_stats(X, y, H, type = "huber")
  rnsc_fit_from_stats(st, h = h, lambda = lambda, H = H)
}


rnsc_predict <- function(fit, newX) {

  newX <- as.matrix(newX)

  if (ncol(newX) != nrow(fit$centroids)) {
    stop("The number of columns in newX does not match the fitted model.")
  }

  nnew <- nrow(newX)
  K <- length(fit$classes)
  scores <- matrix(0, nrow = nnew, ncol = K)

  sel <- fit$selected

  if (length(sel) > 0L) {
    Xs <- newX[, sel, drop = FALSE]
    sc <- fit$scale[sel]

    for (k in seq_len(K)) {
      Z <- sweep(Xs, 2, fit$centroids[sel, k], FUN = "-")
      Z <- sweep(Z, 2, sc, FUN = "/")
      scores[, k] <- rowSums(huber_loss(Z, fit$h))
    }
  }

  for (k in seq_len(K)) {
    scores[, k] <- scores[, k] - 2 * log(fit$prior[k])
  }

  factor(fit$classes[max.col(-scores)], levels = fit$classes)
}


# -------------------------
# 4. Stratified fold assignment
# -------------------------

make_stratified_folds <- function(y, nfold = 10) {
  y <- factor(y)

  if (nfold < 2 || nfold > length(y)) {
    stop("nfold must be between 2 and the number of observations.")
  }

  fold_id <- integer(length(y))

  for (cl in levels(y)) {
    ids <- which(y == cl)
    fold_id[ids] <- sample(rep(seq_len(nfold), length.out = length(ids)))
  }

  fold_id
}


# -------------------------
# 5. RNSC / TF-RNSC tuning by cross-validation
# -------------------------
# type = "huber": two-step CV as in the paper (H first, then h and lambda).
# type = "tf":    tuning-free Huber estimates; only (h, lambda) are tuned.
# Huber estimates are cached per (training fold, H), so they are computed
# once instead of once per (h, lambda).

tune_rnsc <- function(X, y,
                      H_grid = c(1, 2, 3, 5, 10, Inf),
                      h_grid = c(1, 2, 3, 5),
                      lambda_grid = c(0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 6, 8),
                      nfold = 10,
                      fold_id = NULL,
                      type = c("huber", "tf")) {

  type <- match.arg(type)
  X <- as.matrix(X)
  y <- factor(y)

  if (is.null(fold_id)) {
    fold_id <- make_stratified_folds(y, nfold)
  }

  if (length(fold_id) != length(y)) {
    stop("fold_id must have one entry per observation.")
  }

  fold_levels <- sort(unique(fold_id))
  stats_cache <- list()

  get_fold_stats <- function(H) {
    key <- if (type == "tf") "tf" else as.character(H)

    if (is.null(stats_cache[[key]])) {
      stats_cache[[key]] <<- lapply(fold_levels, function(f) {
        tr <- which(fold_id != f)
        rnsc_stats(X[tr, , drop = FALSE], y[tr], H, type = type)
      })
    }

    stats_cache[[key]]
  }

  cv_error <- function(H, h, lambda) {
    fs <- get_fold_stats(H)
    errors <- numeric(length(fold_levels))

    for (i in seq_along(fold_levels)) {
      valid_idx <- which(fold_id == fold_levels[i])

      fit <- rnsc_fit_from_stats(fs[[i]], h = h, lambda = lambda, H = H)
      pred <- rnsc_predict(fit, X[valid_idx, , drop = FALSE])

      errors[i] <- mean(pred != y[valid_idx])
    }

    mean(errors)
  }

  grid1 <- NULL

  if (type == "huber") {
    # Step 1: select H and lambda using squared loss (h = Inf).
    grid1 <- expand.grid(H = H_grid, lambda = lambda_grid)
    grid1$error <- mapply(
      function(H, lambda) cv_error(H, Inf, lambda),
      grid1$H, grid1$lambda
    )

    best1 <- grid1[grid1$error <= min(grid1$error) + 1e-12, , drop = FALSE]
    best1 <- best1[which.max(best1$lambda), , drop = FALSE]
    best_H <- best1$H
  } else {
    best_H <- NA_real_
  }

  # Step 2: fix H and select h and lambda using Huber loss.
  grid2 <- expand.grid(h = h_grid, lambda = lambda_grid)
  grid2$error <- mapply(
    function(h, lambda) cv_error(best_H, h, lambda),
    grid2$h, grid2$lambda
  )

  best2 <- grid2[grid2$error <= min(grid2$error) + 1e-12, , drop = FALSE]
  best2 <- best2[which.max(best2$lambda), , drop = FALSE]

  st_full <- rnsc_stats(X, y, best_H, type = type)
  final_fit <- rnsc_fit_from_stats(
    st_full, h = best2$h, lambda = best2$lambda, H = best_H
  )

  list(
    fit = final_fit,
    H = best_H,
    h = best2$h,
    lambda = best2$lambda,
    cv_step1 = grid1,
    cv_step2 = grid2
  )
}


# -------------------------
# 5b. Componentwise median classifier (Med)
# -------------------------
# Hall, Titterington and Xue (2009): assign a new observation to the class
# with the smallest L1 distance to the class componentwise medians.
#
# ASSUMPTION: Hall et al.'s classifier has no variable selection, but
# Table 1 of Ren & Mai reports C and IC for Med. Here selection is done in
# the NSC style: class medians are soft-thresholded toward the overall
# median, using a MAD-based standard error, with lambda chosen by CV.
# This is a reconstruction, not the authors' code, so Med results may not
# match the paper's Med column.

med_stats <- function(X, y) {

  X <- as.matrix(X)
  y <- factor(y)

  classes <- levels(y)
  K <- length(classes)
  n <- nrow(X)
  p <- ncol(X)

  nk <- as.numeric(table(y))
  prior <- nk / n

  med <- matrix(0, nrow = p, ncol = K)
  mad_k <- matrix(0, nrow = p, ncol = K)

  for (k in seq_len(K)) {
    Xk <- X[y == classes[k], , drop = FALSE]
    med[, k] <- apply(Xk, 2, median)
    mad_k[, k] <- 1.4826 * apply(abs(sweep(Xk, 2, med[, k])), 2, median)
  }

  overall <- apply(X, 2, median)

  scale <- sqrt(rowSums(sweep(mad_k^2, 2, prior, FUN = "*")))
  scale[!is.finite(scale) | scale <= 0] <- NA_real_

  list(
    med = med,
    overall = overall,
    scale = scale,
    nk = nk,
    prior = prior,
    classes = classes,
    n = n,
    p = p
  )
}


med_fit_from_stats <- function(st, lambda) {

  K <- length(st$classes)
  shrunken <- matrix(0, nrow = st$p, ncol = K)

  for (k in seq_len(K)) {
    mk <- sqrt(1 / st$nk[k] - 1 / st$n)

    # 1.2533 = sqrt(pi / 2): asymptotic SE inflation of the sample median.
    step <- 1.2533 * mk * st$scale
    step[!is.finite(step)] <- 0

    d <- (st$med[, k] - st$overall) / step
    d[!is.finite(d)] <- 0

    d_shrunk <- sign(d) * pmax(abs(d) - lambda, 0)
    shrunken[, k] <- st$overall + step * d_shrunk
  }

  selected <- which(apply(abs(shrunken - st$overall) > 1e-8, 1, any))

  list(
    centroids = shrunken,
    selected = selected,
    classes = st$classes,
    lambda = lambda
  )
}


med_predict <- function(fit, newX) {

  newX <- as.matrix(newX)
  K <- length(fit$classes)
  scores <- matrix(0, nrow = nrow(newX), ncol = K)

  sel <- fit$selected

  if (length(sel) > 0L) {
    Xs <- newX[, sel, drop = FALSE]

    for (k in seq_len(K)) {
      scores[, k] <- rowSums(abs(sweep(Xs, 2, fit$centroids[sel, k])))
    }
  }

  # Deterministic tie-breaking (does not touch the random number stream).
  factor(
    fit$classes[max.col(-scores, ties.method = "first")],
    levels = fit$classes
  )
}


tune_med <- function(X, y,
                     lambda_grid = c(0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 6, 8),
                     nfold = 10,
                     fold_id = NULL) {

  X <- as.matrix(X)
  y <- factor(y)

  if (is.null(fold_id)) {
    fold_id <- make_stratified_folds(y, nfold)
  }

  fold_levels <- sort(unique(fold_id))

  fold_stats <- lapply(fold_levels, function(f) {
    tr <- which(fold_id != f)
    med_stats(X[tr, , drop = FALSE], y[tr])
  })

  cv_error <- function(lambda) {
    errors <- numeric(length(fold_levels))

    for (i in seq_along(fold_levels)) {
      valid_idx <- which(fold_id == fold_levels[i])
      fit <- med_fit_from_stats(fold_stats[[i]], lambda)
      pred <- med_predict(fit, X[valid_idx, , drop = FALSE])
      errors[i] <- mean(pred != y[valid_idx])
    }

    mean(errors)
  }

  cv <- vapply(lambda_grid, cv_error, numeric(1))

  best_lambda <- max(lambda_grid[cv <= min(cv) + 1e-12])

  final_fit <- med_fit_from_stats(med_stats(X, y), best_lambda)

  list(
    fit = final_fit,
    lambda = best_lambda,
    cv_error = data.frame(lambda = lambda_grid, error = cv)
  )
}


# -------------------------
# 6. NSC using pamr
# -------------------------

fit_nsc <- function(X_train, y_train, X_test,
                    nfold = 10) {

  if (!requireNamespace("pamr", quietly = TRUE)) {
    stop("Install package 'pamr' first: install.packages('pamr')")
  }

  y_train <- factor(y_train)

  train_data <- list(
    x = t(as.matrix(X_train)),
    y = as.character(y_train)
  )

  fit <- pamr::pamr.train(train_data)
  cv <- pamr::pamr.cv(fit, train_data, nfold = nfold)

  # Match the CV error to the threshold values returned by pamr.
  thresholds <- fit$threshold
  cv_error <- cv$error

  if (is.null(thresholds) || length(thresholds) != length(cv_error)) {
    stop(
      "Could not match pamr CV errors to thresholds. ",
      "Inspect str(fit) and str(cv) for your installed pamr version."
    )
  }

  best_threshold <- max(thresholds[cv_error <= min(cv_error) + 1e-12])

  prediction <- pamr::pamr.predict(
    fit,
    t(as.matrix(X_test)),
    threshold = best_threshold,
    type = "class"
  )

  cen <- pamr::pamr.predict(
    fit,
    t(as.matrix(X_train)),
    threshold = best_threshold,
    type = "centroid"
  )
  cen <- as.matrix(cen)   # expected p x K
  stopifnot(nrow(cen) == ncol(X_train))
  selected <- which(apply(cen, 1, function(r) diff(range(r))) > 1e-8)

  list(
    prediction = factor(as.vector(prediction), levels = levels(y_train)),
    selected = selected,
    fit = fit,
    cv = cv,
    threshold = best_threshold
  )
}


# -------------------------
# 7. L1-penalized logistic regression using glmnet
# -------------------------

fit_logistic <- function(X_train, y_train, X_test,
                         nfold = 10,
                         lambda_choice = c("lambda.min", "lambda.1se")) {

  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("Install package 'glmnet' first: install.packages('glmnet')")
  }

  lambda_choice <- match.arg(lambda_choice)
  X_train <- as.matrix(X_train)
  X_test <- as.matrix(X_test)
  y_train <- factor(y_train)

  if (nlevels(y_train) < 2) stop("At least two response classes are required.")

  if (nlevels(y_train) == 2) {
    fit <- glmnet::cv.glmnet(
      x = X_train,
      y = y_train,
      family = "binomial",
      alpha = 1,
      nfolds = nfold,
      type.measure = "class"
    )
  } else {
    fit <- glmnet::cv.glmnet(
      x = X_train,
      y = y_train,
      family = "multinomial",
      alpha = 1,
      nfolds = nfold,
      type.measure = "class"
    )
  }

  prediction <- predict(
    fit, newx = X_test, s = lambda_choice, type = "class"
  )

  coefficients <- coef(fit, s = lambda_choice)

  if (is.list(coefficients)) {
    selected <- sort(unique(unlist(lapply(coefficients, function(b) {
      which(as.vector(b)[-1] != 0)
    }))))
  } else {
    selected <- which(as.vector(coefficients)[-1] != 0)
  }

  list(
    prediction = factor(as.vector(prediction), levels = levels(y_train)),
    selected = selected,
    fit = fit
  )
}


# -------------------------
# 8. Prediction and selection metrics
# -------------------------

calculate_metrics <- function(prediction, y_test,
                              selected, true_variables) {

  prediction <- as.character(prediction)
  y_test <- as.character(y_test)

  PE <- mean(prediction != y_test) * 100
  C <- sum(selected %in% true_variables)
  IC <- sum(!(selected %in% true_variables))

  c(PE = PE, C = C, IC = IC)
}


# -------------------------
# 9. Run one replication
# -------------------------

run_one_replication <- function(model = c("A", "B"),
                                error_type = 1,
                                p = 2000,
                                n_per_class = 100,
                                test_n_per_class = 1000,
                                rnsc_H_grid = c(1, 2, 3, 5, 10, Inf),
                                rnsc_h_grid = c(1, 2, 3, 5),
                                rnsc_lambda_grid = c(0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 6, 8),
                                nfold = 10) {

  model <- match.arg(model)

  train <- generate_data(
    model, error_type, p, n_per_class
  )

  test <- generate_data(
    model, error_type, p, test_n_per_class
  )

  # One set of CV folds shared by RNSC, TF-RNSC and Med.
  fold_id <- make_stratified_folds(train$y, nfold)

  # RNSC
  rnsc <- tune_rnsc(
    train$X, train$y,
    H_grid = rnsc_H_grid,
    h_grid = rnsc_h_grid,
    lambda_grid = rnsc_lambda_grid,
    nfold = nfold,
    fold_id = fold_id,
    type = "huber"
  )

  met_rnsc <- calculate_metrics(
    rnsc_predict(rnsc$fit, test$X), test$y,
    rnsc$fit$selected, train$true_variables
  )

  # TF-RNSC (tuning-free Huber estimates; h and lambda by CV)
  tf <- tune_rnsc(
    train$X, train$y,
    h_grid = rnsc_h_grid,
    lambda_grid = rnsc_lambda_grid,
    nfold = nfold,
    fold_id = fold_id,
    type = "tf"
  )

  met_tf <- calculate_metrics(
    rnsc_predict(tf$fit, test$X), test$y,
    tf$fit$selected, train$true_variables
  )

  # Med
  med <- tune_med(
    train$X, train$y,
    lambda_grid = rnsc_lambda_grid,
    nfold = nfold,
    fold_id = fold_id
  )

  met_med <- calculate_metrics(
    med_predict(med$fit, test$X), test$y,
    med$fit$selected, train$true_variables
  )

  # NSC
  nsc <- fit_nsc(
    train$X, train$y, test$X, nfold = nfold
  )

  met_nsc <- calculate_metrics(
    nsc$prediction, test$y,
    nsc$selected, train$true_variables
  )

  # Logistic
  logistic <- fit_logistic(
    train$X, train$y, test$X, nfold = nfold
  )

  met_logistic <- calculate_metrics(
    logistic$prediction, test$y,
    logistic$selected, train$true_variables
  )

  out <- c(met_tf, met_rnsc, met_nsc, met_med, met_logistic)
  names(out) <- paste0(
    rep(c("TF_RNSC", "RNSC", "NSC", "Med", "Logistic"), each = 3),
    "_", c("PE", "C", "IC")
  )
  out
}


# -------------------------
# 10. Reproducible per-replication seeds and seeded wrapper
# -------------------------
# Each replication gets its own L'Ecuyer-CMRG stream, so results are
# identical whether you use 1 core or many, and in any scheduling order.

make_rep_seeds <- function(n_seeds, seed) {
  old_kind <- RNGkind()
  on.exit(RNGkind(old_kind[1], old_kind[2], old_kind[3]), add = TRUE)

  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  s <- get(".Random.seed", envir = globalenv())

  seeds <- vector("list", n_seeds)
  for (i in seq_len(n_seeds)) {
    seeds[[i]] <- s
    s <- parallel::nextRNGStream(s)
  }
  seeds
}

# Runs one replication with a given RNG stream. Errors are caught and
# returned as an NA row (with the message stored in attr "error"), so one
# failed replication does not kill a long run.
run_one_replication_seeded <- function(seed_vec, ...) {
  RNGkind("L'Ecuyer-CMRG")
  assign(".Random.seed", seed_vec, envir = globalenv())

  tryCatch(
    as.numeric(run_one_replication(...)),
    error = function(e) {
      structure(rep(NA_real_, 15), error = conditionMessage(e))
    }
  )
}


# -------------------------
# 11. Run all simulation settings (optionally in parallel)
# -------------------------

run_table1_simulation <- function(B = 200,
                                  models = c("A", "B"),
                                  error_types = 1:4,
                                  p = 2000,
                                  n_per_class = 100,
                                  test_n_per_class = 1000,
                                  seed = 123,
                                  rnsc_H_grid = c(1, 2, 3, 5, 10, Inf),
                                  rnsc_h_grid = c(1, 2, 3, 5),
                                  rnsc_lambda_grid = c(0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 6, 8),
                                  nfold = 10,
                                  progress = TRUE,
                                  n_cores = 1L,
                                  checkpoint_file = NULL) {

  col_names <- paste0(
    rep(c("TF_RNSC", "RNSC", "NSC", "Med", "Logistic"), each = 3),
    "_", c("PE", "C", "IC")
  )

  n_settings <- length(models) * length(error_types)
  all_seeds <- make_rep_seeds(B * n_settings, seed)

  # Start a cluster once and reuse it for all settings.
  cl <- NULL
  if (n_cores > 1L) {
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    # Send every function defined in the global environment to the workers
    # (packages are called with pamr:: / glmnet:: so no library() is needed).
    fn_names <- Filter(
      function(nm) is.function(get(nm, envir = globalenv())),
      ls(globalenv())
    )
    parallel::clusterExport(cl, varlist = fn_names, envir = globalenv())
  }

  all_results <- list()
  counter <- 1L

  for (model in models) {
    for (err in error_types) {

      setting <- paste0(model, ".", err)
      idx <- ((counter - 1L) * B + 1L):(counter * B)
      rep_seeds <- all_seeds[idx]

      if (progress) {
        message("Running setting ", setting, " (", B,
                " replications, ", n_cores, " core(s))")
      }
      t0 <- Sys.time()

      res_list <- if (is.null(cl)) {
        lapply(
          rep_seeds, run_one_replication_seeded,
          model = model,
          error_type = err,
          p = p,
          n_per_class = n_per_class,
          test_n_per_class = test_n_per_class,
          rnsc_H_grid = rnsc_H_grid,
          rnsc_h_grid = rnsc_h_grid,
          rnsc_lambda_grid = rnsc_lambda_grid,
          nfold = nfold
        )
      } else {
        parallel::parLapplyLB(
          cl, rep_seeds, run_one_replication_seeded,
          model = model,
          error_type = err,
          p = p,
          n_per_class = n_per_class,
          test_n_per_class = test_n_per_class,
          rnsc_H_grid = rnsc_H_grid,
          rnsc_h_grid = rnsc_h_grid,
          rnsc_lambda_grid = rnsc_lambda_grid,
          nfold = nfold
        )
      }

      errs <- unlist(lapply(res_list, function(r) attr(r, "error")))
      if (length(errs) > 0) {
        warning(
          "Setting ", setting, ": ", length(errs), " of ", B,
          " replications failed. First error: ", errs[1],
          call. = FALSE
        )
      }

      reps <- matrix(
        unlist(res_list), nrow = B, byrow = TRUE,
        dimnames = list(NULL, col_names)
      )

      all_results[[counter]] <- list(
        model = setting,
        results = reps
      )

      if (!is.null(checkpoint_file)) {
        saveRDS(all_results, checkpoint_file)
      }

      if (progress) {
        message("  Finished ", setting, " in ",
                round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1),
                " min")
      }

      counter <- counter + 1L
    }
  }

  all_results
}


# -------------------------
# 12. Summarize means and standard errors
# -------------------------

summarize_table1 <- function(results) {

  rows <- list()
  counter <- 1L

  for (item in results) {
    x <- item$results
    setting <- item$model

    for (method in c("TF_RNSC", "RNSC", "NSC", "Med", "Logistic")) {
      for (metric in c("PE", "C", "IC")) {

        col <- paste0(method, "_", metric)
        if (!col %in% colnames(x)) next

        values <- x[, col]
        values <- values[is.finite(values)]

        rows[[counter]] <- data.frame(
          Model = setting,
          Method = method,
          Metric = metric,
          Mean = mean(values),
          SE = if (length(values) > 1) {
            stats::sd(values) / sqrt(length(values))
          } else {
            NA_real_
          }
        )

        counter <- counter + 1L
      }
    }
  }

  do.call(rbind, rows)
}


# -------------------------
# 13. Optional helper: save summary as CSV
# -------------------------

save_table1_summary <- function(results,
                                file = "Table1_simulation_results.csv") {
  summary_table <- summarize_table1(results)
  utils::write.csv(summary_table, file, row.names = FALSE)
  invisible(summary_table)
}

