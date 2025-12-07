#' Pool results with Robins-Wang variance
#' 
#' @param object A with_rw object from with_rw()
#' @param ... Additional arguments (currently unused)
#' @return A pool_rw object with pooled estimates and Robins-Wang SEs
#' 
#' @examples
#' \dontrun{
#' library(mice)
#' imp <- mice(nhanes, method = "norm", m = 5, tasks = "train", print = FALSE)
#' fit <- with_rw(imp, lm(bmi ~ age + hyp))
#' pooled <- pool_rw(fit)
#' summary(pooled)
#' }
#' @export
pool_rw <- function(object, ...) {
  if (!inherits(object, "with_rw")) {
    cli::cli_abort(
      "{.arg object} must be a {.cls with_rw} object, not {.cls {class(object)}}."
    )
  }
  
  m <- object$m
  n <- object$n
  results <- object$results
  
  # Extract coefficients
  coefs <- vapply(results, function(r) stats::coef(r$model), 
                  numeric(length(stats::coef(results[[1]]$model))))
  est <- rowMeans(coefs)
  
  # Compute RW variance components
  U_sum <- Reduce(`+`, lapply(results, function(r) r$U))
  tau_sum <- Reduce(`+`, lapply(results, function(r) r$tau))
  
  u_bar <- U_sum / m
  omega <- crossprod(u_bar) / n
  
  kappa_sum <- 0
  alpha_sum <- 0
  d_bar_sum <- 0
  
  for (p in seq_len(m)) {
    U <- results[[p]]$U
    S_mis_imp <- results[[p]]$S_mis_imp
    d <- results[[p]]$d
    
    kappa_sum <- kappa_sum + crossprod(U, S_mis_imp)
    alpha_sum <- alpha_sum + crossprod(d)
    d_bar_sum <- d_bar_sum + d
  }
  
  kappa <- kappa_sum / (n * m)
  alpha <- alpha_sum / (n * m)
  d_bar <- d_bar_sum / m
  
  delta <- omega + kappa %*% alpha %*% t(kappa) +
    (1 / n) * (kappa %*% t(d_bar) %*% u_bar + t(kappa %*% t(d_bar) %*% u_bar))
  
  tau <- tau_sum / (m * n)
  tau_inv <- solve(tau)
  GAMMA <- (1 / n) * tau_inv %*% delta %*% t(tau_inv)
  
  se <- sqrt(diag(GAMMA))
  
  # Create output similar to mice::pool
  out <- data.frame(
    term = names(est),
    estimate = est,
    std.error = se,
    statistic = est / se,
    p.value = 2 * stats::pnorm(-abs(est / se)),
    conf.low = est - 1.96 * se,
    conf.high = est + 1.96 * se
  )
  rownames(out) <- NULL
  
  structure(
    list(
      pooled = out,
      variance = GAMMA,
      m = m,
      n = n,
      call = match.call()
    ),
    class = "pool_rw"
  )
}

# S3 methods
#' @export
print.pool_rw <- function(x, ...) {
  cli::cli_h1("Robins-Wang Pooled Results")
  cli::cli_text("Number of imputations: {x$m}")
  cli::cli_text("Sample size: {x$n}")
  cli::cli_text("")
  print(x$pooled, digits = 4)
  invisible(x)
}

#' @export
summary.pool_rw <- function(object, ...) {
  print(object, ...)
}

#' @export
coef.pool_rw <- function(object, ...) {
  stats::setNames(object$pooled$estimate, object$pooled$term)
}

#' @export
vcov.pool_rw <- function(object, ...) {
  object$variance
}