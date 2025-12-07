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
  
  est <- pool_estimates(results)
  GAMMA <- compute_rw_variance(results, m, n)
  se <- sqrt(diag(GAMMA))
  
  out <- construct_pooled_output(est, se)
  
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