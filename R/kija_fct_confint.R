#' Confidence set for functions of model parameters
#'
#' Computes confidence sets of functions of model parameters by computing a
#' confidence set of the model parameters and returning the codomain of the
#' provided function given the confidence set of model parameters as domain.
#'
#' @param object A fitted model object.
#' @param f A function taking the parameter vector as its single argument, and
#'   returning a numeric vector.
#' @param which_parm Either a logical vector the same length as the coefficient
#'   vector, with `TRUE` indicating a coefficient is used by `f`, or an integer
#'   vector with the indices of the coefficients used by `f`.
#' @param level The confidence level required.
#' @param ... Additional argument(s) passed to methods.
#'
#' @returns A tibble with columns estimate, conf.low, and conf.high or if
#'   return_beta is `TRUE`, a list with the tibble and the beta values on the
#'   boundary used to calculate the confidence limits.
#'
#' @details Assume the response Y and predictors X are given by a generalized
#'   linear model, that is, they fulfill the assumptions
#'   \deqn{E(Y|X)=\mu(X^T\beta)}{E(Y|X)=\mu(X^T\beta)}
#'   \deqn{V(Y|X)=\psi \nu(\mu(X^T\beta))}{V(Y|X)=\psi \nu(\mu(X^T\beta))}
#'   \deqn{Y|X\sim\varepsilon(\theta,\nu_{\psi}).}{%
#'   Y|X ~ \varepsilon(\theta,\nu_(\psi)).}
#'   Here \eqn{\mu} is the mean value function, \eqn{\nu} is the variance
#'   function, and \eqn{\psi} is the dispersion parameter in the exponential
#'   dispersion model
#'   \eqn{\varepsilon(\theta,\nu_{\psi})}{\varepsilon(\theta,\nu_(\psi))}, where
#'   \eqn{\theta} is the canonical parameter and \eqn{\nu_{\psi}}{\nu_(\psi)} is
#'   the structure measure. Then it follows from the central limit theorem that
#'   \deqn{hat(\beta)\sim N(\beta, (X^TWX)^{-1})}{%
#'   hat(\beta)~N(\beta, (X^TWX)^(-1))}
#'   will be a good approximation in large samples, where \eqn{X^TWX} is the
#'   Fisher information of the exponential dispersion model.
#'
#'   From this, the combinant
#'   \deqn{(hat(\beta)-\beta)^TX^TWX(hat(\beta)-\beta)}{%
#'   (hat(\beta)-\beta)^TX^TWX(hat(\beta)-\beta)}
#'   is an approximate pivot, with a \eqn{\chi_p^2} distribution. Then
#'   \deqn{C_{\beta}=\{\beta|(hat(\beta)-\beta)^TX^TWX(hat(\beta)-\beta)<\chi_p^2(1-\alpha)\}}{%
#'   C_(\beta)=\{\beta|(hat(\beta)-\beta)^TX^TWX(hat(\beta)-\beta)<\chi_p^2(1-\alpha)\}}
#'   is an approximate \eqn{(1-\alpha)}-confidence set for the parameter vector
#'   \eqn{\beta}. Similarly, confidence sets for sub-vectors of \eqn{\beta} can
#'   be obtained by the fact that marginal distributions of normal distributions
#'   are again normally distributed, where the mean vector and covariance matrix
#'   are appropriate subvectors and submatrices.
#'
#'   Finally, a confidence set for the transformed parameters \eqn{f(\beta)}
#'   is obtained as
#'   \deqn{\{f(\beta)|\beta\in C_{\beta}\}}{\{f(\beta)|\beta\in C_(\beta)\}.}
#'   Note this is a conservative confidence set, since parameters outside the
#'   confidence set of \eqn{\beta} can be mapped to the confidence set of the
#'   transformed parameter.
#'
#'   To determine \eqn{C_{\beta}}, `fct_confint()` finds the boundary by taking
#'   a number of points around \eqn{hat(\beta)} and projecting them onto the
#'   boundary. Therefore, the confidence set of the transformed parameter will
#'   only be valid if the boundary of \eqn{C_{\beta}} is mapped to the boundary
#'   of the confidence set for the transformed parameter.
#'
#'   The points projected to the boundary are either laid out in a grid around
#'   \eqn{hat(\beta)}, with the number of points in each direction determined
#'   by `n_grid`, or uniformly at random on a hypersphere, with the number of
#'   points determined by `k`. The radius of the grid/sphere is determined by
#'   `len`.
#'
#' @author
#' KIJA
#'
#' @examples
#' 1+1
#'
#' @export

fct_confint <- function(
    object,
    f,
    which_parm = rep(TRUE, length(coef(object))),
    level = 0.95,
    ...
) {
  UseMethod("fct_confint")
}


#' @rdname fct_confint
#' @param len numeric, the radius of the sphere or box used to define directions
#'   in which to look for boundary points of the parameter confidence set.
#' @param n_grid Either `NULL` or an integer vector of length 1 or the number of
#'   `TRUE`/indices in which_parm. Specifies the number of grid points in
#'   each dimension of a grid with endpoints defined by len. If `NULL` or `0L`,
#'   will instead sample k points uniformly on a sphere.
#' @param k If n_grid is `NULL` or `0L`, the number of points to sample
#'   uniformly from a sphere.
#' @param parallel Logical, if `TRUE` parallel computing is used when solving
#'   for points on the boundary of the parameter confidence set.
#' @param n_cores An integer specifying the number of threads to use for
#'   parallel computing.
#' @param return_beta Logical, if `TRUE` returns both the confidence limits and
#'   the parameter values used from the boundary of the parameter confidence
#'   set.
#' @param verbose Logical, if `TRUE` prints information about the number of
#'   points on the boundary used to calculate the confidence limits.
#'
#' @export

fct_confint.lm <- function(
    object,
    f,
    which_parm = rep(TRUE, length(coef(object))),
    level = 0.95,
    len = 0.1,
    n_grid = 0L,
    k = 1000L,
    parallel = FALSE,
    n_cores = 10,
    return_beta = FALSE,
    verbose = FALSE,
    ...
) {
  if (parallel) {
    nCores <- min(
      parallel::detectCores(),
      n_cores
    )
    cluster <- parallel::makeCluster(nCores)
    registerDoParallel(cluster)
    on.exit(parallel::stopCluster(cluster))
  }
  # convert which_parm to a logical vector
  if (!rlang::is_logical(which_parm)) {
    if (rlang::is_double(which_parm)) {
      which_parm_ind <- as.integer(round(which_parm))
      rlang::warn(
        "A double was provided and was rounded and converted to an integer.",
        "It is strongly recommended to use an integer vector with indices to",
        "avoid unexpected behavior."
      )
    } else if (
      length(which_parm) <= length(object$coefficients) &&
      rlang::is_integer(which_parm) &&
      all(which_parm >= 1L) &&
      all(which_parm <= length(object$coefficients))
    ) {
      which_parm_ind <- which_parm
    } else {
      rlang::abort(
        paste0(
          "which_parm must either be a logical vector the same length as the",
          "coefficients or an integer vector with indices of the coefficient",
          "vector."
        )
      )
    }
    which_parm <- rep(FALSE, length(object$coefficients))
    which_parm[which_parm_ind] <- TRUE
  }
  # extract MLE parameters and standard deviation
  beta_hat <- coef(object)
  # contruct points around the estimated parameter vector defining directions
  # to look for points on the boundary of the confidence set.
  if (is.null(n_grid) || any(n_grid) == 0L) {
    # create k point uniformly on a sum(which_parm)-dimensional sphere
    delta <- matrix(rnorm(sum(which_parm) * k), ncol = sum(which_parm))
    norm <- apply(delta, 1, function(y) sqrt(sum(y^2)))
    delta <- t(delta / rep(norm / len, ncol(delta)))
    rownames(delta) <- names(beta_hat)[which_parm]
  } else if (!rlang::is_integer(n_grid)) {
    rlang::abort("n_grid must be an integer vector.")
  } else if (!(match(length(n_grid), c(1L, sum(which_parm))))) {
    rlang::abort(glue::glue("n_grid must have length 1 or {sum(which_parm)}."))
  } else {
    # create sum(which_parm)-dimensional grid with n_grid points in each direction
    if (length(n_grid) == 1L) {
      delta <- as.matrix(expand.grid(
        rep(list(seq(-len, len, length.out = n_grid)), sum(which_parm))
      ))
    } else if(length(n_grid) == sum(which_parm)) {
      delta <- as.matrix(expand.grid(
        lapply(n_grid, function(y) seq(-len, len, length.out = y))
      ))
    } else {
      rlang::abort(glue::glue("n_grid must be length 1 or {sum(which_parm)}"))
    }
    delta <- t(delta)
    rownames(delta) <- names(beta_hat)[which_parm]
  }
  # design matrix
  X <- model.matrix(object)[, which_parm, drop = FALSE]
  # weight
  rdf <- object$df.residual
  if (is.null(object$weights)) {
    rss <- sum(object$residuals^2)
  } else {
    rss <- sum(object$weights * object$residuals^2)
  }
  W <- rss / rdf # equivalent to summary(object)$sigma^2
  # solve equation for scaling of points to end up on boundary of confidence set
  xtx_inv <- W * solve(crossprod(X))
  xtx_red <- solve(xtx_inv[which_parm, which_parm])
  if(parallel) {
    fn <- function(y) {
      t(delta[, y, drop = FALSE]) %*% xtx_red %*% delta[, y, drop = FALSE]
    }
    a <- foreach(y = seq_len(ncol(delta))) %dopar% fn(y)
    a <- unlist(a)
  } else {
    a <- vapply(
      seq_len(ncol(delta)),
      function(y) {
        t(delta[, y, drop = FALSE]) %*% xtx_red %*% delta[, y, drop = FALSE]
      },
      double(1L)
    )
  }
  c <- rep(
    qchisq(level, length(beta_hat[which_parm])),
    ncol(delta)
  )
  suppressWarnings(
    alpha1 <- -sqrt(c / a)
  )
  suppressWarnings(
    alpha2 <- sqrt(c / a)
  )
  # remove delta values that lead to complex solutions
  which_alpha <- which(!is.nan(alpha1))
  alpha1 <- alpha1[which_alpha]
  alpha2 <- alpha2[which_alpha]
  delta_red <- delta[, which_alpha, drop = FALSE]
  if (verbose) cat("points on the boundary:", 2 * length(which_alpha), "\n")
  # determine coefficient values on boundary of confidence set
  beta_1 <- matrix(
    rep(beta_hat[which_parm], length(which_alpha)),
    byrow = TRUE,
    ncol = length(beta_hat[which_parm])
  )
  beta_1 <- beta_1 +
    t(
      matrix(
        rep(alpha1, nrow(delta_red)),
        byrow = TRUE,
        nrow = nrow(delta_red)
      ) * delta_red
    )
  colnames(beta_1) <- names(beta_hat)[which_parm]
  beta_2 <- matrix(
    rep(beta_hat[which_parm], length(which_alpha)),
    byrow = TRUE,
    ncol = length(beta_hat[which_parm])
  )
  beta_2 <- beta_2 +
    t(
      matrix(
        rep(alpha2, nrow(delta_red)),
        byrow = TRUE,
        nrow = nrow(delta_red)
      ) * delta_red
    )
  colnames(beta_2) <- names(beta_hat)[which_parm]
  suppressMessages(
    ci_data_1 <- lapply(
      seq_len(nrow(beta_1)),
      function(x) f(beta_1[x, ])
    ) |>
      structure(
        names = paste0("ci_bound_", seq_len(nrow(beta_1)))
      ) |>
      tibble::as_tibble() |>
      dplyr::bind_cols(
        tibble::tibble(
          estimate = f(beta_hat[which_parm])
        )
      )
  )
  suppressMessages(
    ci_data_2 <- lapply(
      seq_len(nrow(beta_2)),
      function(x) f(beta_2[x, ])
    ) |>
      structure(
        names = paste0("ci_bound_", seq_len(nrow(beta_2)) + nrow(beta_1))
      ) |>
      tibble::as_tibble()
  )
  # combine negative and positive solutions
  ci_data <- dplyr::bind_cols(ci_data_1, ci_data_2) |>
    dplyr::select("estimate", tidyselect::everything())

  ci_data$conf.low <- do.call(pmin, ci_data[, -1])
  ci_data$conf.high <- do.call(pmax, ci_data[, -1])

  if(return_beta) {
    ci_data <- ci_data |>
      dplyr::select("estimate", "conf.low", "conf.high", tidyselect::everything())
    return(list(ci_data = ci_data, beta = dplyr::bind_rows(beta_1, beta_2)))
  } else {
    ci_data <- ci_data |>
      dplyr::select("estimate", "conf.low", "conf.high")
    return(ci_data)
  }
}

#' @rdname fct_confint
#' @export

fct_confint.glm <- function(
    object,
    f,
    which_parm = rep(TRUE, length(coef(object))),
    level = 0.95,
    len = 0.1,
    n_grid = 0L,
    k = 1000L,
    parallel = FALSE,
    n_cores = 10,
    return_beta = FALSE,
    verbose = FALSE,
    ...
) {
  if (parallel) {
    nCores <- min(
      parallel::detectCores(),
      n_cores
    )
    cluster <- parallel::makeCluster(nCores)
    registerDoParallel(cluster)
    on.exit(parallel::stopCluster(cluster))
  }
  # convert which_parm to a logical vector
  if (!rlang::is_logical(which_parm)) {
    if (rlang::is_double(which_parm)) {
      which_parm_ind <- as.integer(round(which_parm))
      rlang::warn(
        "A double was provided and was rounded and converted to an integer.",
        "It is strongly recommended to use an integer vector with indices to",
        "avoid unexpected behavior."
      )
    } else if (
      length(which_parm) <= length(object$coefficients) &&
      rlang::is_integer(which_parm) &&
      all(which_parm >= 1L) &&
      all(which_parm <= length(object$coefficients))
    ) {
      which_parm_ind <- which_parm
    } else {
      rlang::abort(
        paste0(
          "which_parm must either be a logical vector the same length as the",
          "coefficients or an integer vector with indices of the coefficient",
          "vector."
        )
      )
    }
    which_parm <- rep(FALSE, length(object$coefficients))
    which_parm[which_parm_ind] <- TRUE
  }
  # extract MLE parameters and standard deviation
  beta_hat <- coef(object)
  # contruct points around the estimated parameter vector defining directions
  # to look for points on the boundary of the confidence set.
  if (is.null(n_grid) || any(n_grid) == 0L) {
    # create k point uniformly on a sum(which_parm)-dimensional sphere
    delta <- matrix(rnorm(sum(which_parm) * k), ncol = sum(which_parm))
    norm <- apply(delta, 1, function(x) sqrt(sum(x^2)))
    delta <- t(delta / rep(norm / len, ncol(delta)))
    rownames(delta) <- names(beta_hat)[which_parm]
  } else if (!rlang::is_integer(n_grid)) {
    rlang::abort("n_grid must be an integer vector.")
  } else if (!(length(n_grid) %in% c(1L, sum(which_parm)))) {
    rlang::abort(glue::glue("n_grid must have length 1 or {sum(which_parm)}."))
  } else {
    # create sum(which_parm)-dimensional grid with n_grid points in each direction
    if (length(n_grid) == 1L) {
      delta <- as.matrix(expand.grid(
        rep(list(seq(-len, len, length.out = n_grid)), sum(which_parm))
      ))
    } else if(length(n_grid) == sum(which_parm)) {
      delta <- as.matrix(expand.grid(
        lapply(n_grid, function(x) seq(-len, len, length.out = x))
      ))
    } else {
      rlang::abort(glue::glue("n_grid must be length 1 or {sum(which_parm)}"))
    }
    delta <- t(delta)
    rownames(delta) <- names(beta_hat)[which_parm]
  }
  # design matrix
  X <- model.matrix(object)[, which_parm, drop = FALSE]
  # weight
  W <- object$weights
  # solve equation for scaling of points to end up on boundary of confidence set
  xtx_inv <- solve(crossprod(X * sqrt(W)))
  xtx_red <- solve(xtx_inv[which_parm, which_parm])
  if(parallel) {
    fn <- function(y) {
      t(delta[, y, drop = FALSE]) %*% xtx_red %*% delta[, y, drop = FALSE]
    }
    a <- a <- foreach(y = seq_len(ncol(delta))) %dopar% fn(y)
    a <- unlist(a)
  } else {
    a <- vapply(
      seq_len(ncol(delta)),
      function(y) {
        t(delta[, y, drop = FALSE]) %*% xtx_red %*% delta[, y, drop = FALSE]
      },
      double(1L)
    )
  }
  c <- rep(
    qchisq(level, length(beta_hat[which_parm])),
    ncol(delta)
  )
  suppressWarnings(
    alpha1 <- -sqrt(c / a)
  )
  suppressWarnings(
    alpha2 <- sqrt(c / a)
  )
  # remove delta values that lead to complex solutions
  which_alpha <- which(!is.nan(alpha1))
  alpha1 <- alpha1[which_alpha]
  alpha2 <- alpha2[which_alpha]
  delta_red <- delta[, which_alpha, drop = FALSE]
  if (verbose) cat("points on the boundary:", 2 * length(which_alpha), "\n")
  # determine coefficient values on boundary of confidence set
  beta_1 <- matrix(
    rep(beta_hat[which_parm], length(which_alpha)),
    byrow = TRUE,
    ncol = length(beta_hat[which_parm])
  )
  beta_1 <- beta_1 +
    t(
      matrix(
        rep(alpha1, nrow(delta_red)),
        byrow = TRUE,
        nrow = nrow(delta_red)
      ) * delta_red
    )
  colnames(beta_1) <- names(beta_hat)[which_parm]
  beta_2 <- matrix(
    rep(beta_hat[which_parm], length(which_alpha)),
    byrow = TRUE,
    ncol = length(beta_hat[which_parm])
  )
  beta_2 <- beta_2 +
    t(
      matrix(
        rep(alpha2, nrow(delta_red)),
        byrow = TRUE,
        nrow = nrow(delta_red)
      ) * delta_red
    )
  colnames(beta_2) <- names(beta_hat)[which_parm]
  suppressMessages(
    ci_data_1 <- lapply(
      seq_len(nrow(beta_1)),
      function(x) f(beta_1[x, ])
    ) |>
      structure(
        names = paste0("ci_bound_", seq_len(nrow(beta_1)))
      ) |>
      tibble::as_tibble() |>
      dplyr::bind_cols(
        tibble::tibble(
          estimate = f(beta_hat[which_parm])
        )
      )
  )
  suppressMessages(
    ci_data_2 <- lapply(
      seq_len(nrow(beta_2)),
      function(x) f(beta_2[x, ])
    ) |>
      structure(
        names = paste0("ci_bound_", seq_len(nrow(beta_2)) + nrow(beta_1))
      ) |>
      tibble::as_tibble()
  )
  # combine negative and positive solutions
  ci_data <- dplyr::bind_cols(ci_data_1, ci_data_2) |>
    dplyr::select("estimate", tidyselect::everything())

  ci_data$conf.low <- do.call(pmin, ci_data[, -1])
  ci_data$conf.high <- do.call(pmax, ci_data[, -1])

  if(return_beta) {
    ci_data <- ci_data |>
      dplyr::select("estimate", "conf.low", "conf.high", tidyselect::everything())
    return(list(ci_data = ci_data, beta = dplyr::bind_rows(beta_1, beta_2)))
  } else {
    ci_data <- ci_data |>
      dplyr::select("estimate", "conf.low", "conf.high")
    return(ci_data)
  }
}

#' @rdname fct_confint
#' @export
# NOTE: the lms method is currently a copy of the lm method.
fct_confint.lms <- function(
    object,
    f,
    which_parm = rep(TRUE, length(coef(object))),
    level = 0.95,
    len = 0.1,
    n_grid = 0L,
    k = 1000L,
    parallel = FALSE,
    n_cores = 10,
    return_beta = FALSE,
    verbose = FALSE,
    ...
) {
  if (parallel) {
    nCores <- min(
      parallel::detectCores(),
      n_cores
    )
    cluster <- parallel::makeCluster(nCores)
    registerDoParallel(cluster)
    on.exit(parallel::stopCluster(cluster))
  }
  # convert which_parm to a logical vector
  if (!rlang::is_logical(which_parm)) {
    if (rlang::is_double(which_parm)) {
      which_parm_ind <- as.integer(round(which_parm))
      rlang::warn(
        "A double was provided and was rounded and converted to an integer.",
        "It is strongly recommended to use an integer vector with indices to",
        "avoid unexpected behavior."
      )
    } else if (
      length(which_parm) <= length(object$coefficients) &&
      rlang::is_integer(which_parm) &&
      all(which_parm >= 1L) &&
      all(which_parm <= length(object$coefficients))
    ) {
      which_parm_ind <- which_parm
    } else {
      rlang::abort(
        paste0(
          "which_parm must either be a logical vector the same length as the",
          "coefficients or an integer vector with indices of the coefficient",
          "vector."
        )
      )
    }
    which_parm <- rep(FALSE, length(object$coefficients))
    which_parm[which_parm_ind] <- TRUE
  }
  # extract MLE parameters and standard deviation
  beta_hat <- coef(object)
  # contruct points around the estimated parameter vector defining directions
  # to look for points on the boundary of the confidence set.
  if (is.null(n_grid) || any(n_grid) == 0L) {
    # create k point uniformly on a sum(which_parm)-dimensional sphere
    delta <- matrix(rnorm(sum(which_parm) * k), ncol = sum(which_parm))
    norm <- apply(delta, 1, function(y) sqrt(sum(y^2)))
    delta <- t(delta / rep(norm / len, ncol(delta)))
    rownames(delta) <- names(beta_hat)[which_parm]
  } else if (!rlang::is_integer(n_grid)) {
    rlang::abort("n_grid must be an integer vector.")
  } else if (!(match(length(n_grid), c(1L, sum(which_parm))))) {
    rlang::abort(glue::glue("n_grid must have length 1 or {sum(which_parm)}."))
  } else {
    # create sum(which_parm)-dimensional grid with n_grid points in each direction
    if (length(n_grid) == 1L) {
      delta <- as.matrix(expand.grid(
        rep(list(seq(-len, len, length.out = n_grid)), sum(which_parm))
      ))
    } else if(length(n_grid) == sum(which_parm)) {
      delta <- as.matrix(expand.grid(
        lapply(n_grid, function(y) seq(-len, len, length.out = y))
      ))
    } else {
      rlang::abort(glue::glue("n_grid must be length 1 or {sum(which_parm)}"))
    }
    delta <- t(delta)
    rownames(delta) <- names(beta_hat)[which_parm]
  }
  # design matrix
  X <- model.matrix(object)[, which_parm, drop = FALSE]
  # weight
  rdf <- object$df.residual
  if (is.null(object$weights)) {
    rss <- sum(object$residuals^2)
  } else {
    rss <- sum(object$weights * object$residuals^2)
  }
  W <- rss / rdf # equivalent to summary(object)$sigma^2
  # solve equation for scaling of points to end up on boundary of confidence set
  xtx_inv <- W * solve(crossprod(X))
  xtx_red <- solve(xtx_inv[which_parm, which_parm])
  if(parallel) {
    fn <- function(y) {
      t(delta[, y, drop = FALSE]) %*% xtx_red %*% delta[, y, drop = FALSE]
    }
    a <- foreach(y = seq_len(ncol(delta))) %dopar% fn(y)
    a <- unlist(a)
  } else {
    a <- vapply(
      seq_len(ncol(delta)),
      function(y) {
        t(delta[, y, drop = FALSE]) %*% xtx_red %*% delta[, y, drop = FALSE]
      },
      double(1L)
    )
  }
  c <- rep(
    qchisq(level, length(beta_hat[which_parm])),
    ncol(delta)
  )
  suppressWarnings(
    alpha1 <- -sqrt(c / a)
  )
  suppressWarnings(
    alpha2 <- sqrt(c / a)
  )
  # remove delta values that lead to complex solutions
  which_alpha <- which(!is.nan(alpha1))
  alpha1 <- alpha1[which_alpha]
  alpha2 <- alpha2[which_alpha]
  delta_red <- delta[, which_alpha, drop = FALSE]
  if (verbose) cat("points on the boundary:", 2 * length(which_alpha), "\n")
  # determine coefficient values on boundary of confidence set
  beta_1 <- matrix(
    rep(beta_hat[which_parm], length(which_alpha)),
    byrow = TRUE,
    ncol = length(beta_hat[which_parm])
  )
  beta_1 <- beta_1 +
    t(
      matrix(
        rep(alpha1, nrow(delta_red)),
        byrow = TRUE,
        nrow = nrow(delta_red)
      ) * delta_red
    )
  colnames(beta_1) <- names(beta_hat)[which_parm]
  beta_2 <- matrix(
    rep(beta_hat[which_parm], length(which_alpha)),
    byrow = TRUE,
    ncol = length(beta_hat[which_parm])
  )
  beta_2 <- beta_2 +
    t(
      matrix(
        rep(alpha2, nrow(delta_red)),
        byrow = TRUE,
        nrow = nrow(delta_red)
      ) * delta_red
    )
  colnames(beta_2) <- names(beta_hat)[which_parm]
  suppressMessages(
    ci_data_1 <- lapply(
      seq_len(nrow(beta_1)),
      function(x) f(beta_1[x, ])
    ) |>
      structure(
        names = paste0("ci_bound_", seq_len(nrow(beta_1)))
      ) |>
      tibble::as_tibble() |>
      dplyr::bind_cols(
        tibble::tibble(
          estimate = f(beta_hat[which_parm])
        )
      )
  )
  suppressMessages(
    ci_data_2 <- lapply(
      seq_len(nrow(beta_2)),
      function(x) f(beta_2[x, ])
    ) |>
      structure(
        names = paste0("ci_bound_", seq_len(nrow(beta_2)) + nrow(beta_1))
      ) |>
      tibble::as_tibble()
  )
  # combine negative and positive solutions
  ci_data <- dplyr::bind_cols(ci_data_1, ci_data_2) |>
    dplyr::select("estimate", tidyselect::everything())

  ci_data$conf.low <- do.call(pmin, ci_data[, -1])
  ci_data$conf.high <- do.call(pmax, ci_data[, -1])

  if(return_beta) {
    ci_data <- ci_data |>
      dplyr::select("estimate", "conf.low", "conf.high", tidyselect::everything())
    return(list(ci_data = ci_data, beta = dplyr::bind_rows(beta_1, beta_2)))
  } else {
    ci_data <- ci_data |>
      dplyr::select("estimate", "conf.low", "conf.high")
    return(ci_data)
  }
}

#' @rdname fct_confint
#' @export
fct_confint.default <- function(
    object,
    f,
    which_parm = rep(TRUE, length(coef(object))),
    level = 0.95,
    ...
) {
  rlang::abort("Class is not supported. Object must have class lm, glm or lms.")
}
