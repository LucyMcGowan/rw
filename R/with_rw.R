#' Fit analysis model with Robins-Wang variance calculation
#' 
#' @param data A mids object from mice with parametric imputation methods
#' @param expr Expression to evaluate on each imputed dataset (typically a model formula)
#' @param ... Additional arguments passed to the analysis function
#' @return A with_rw object containing fitted models and RW components
#' 
#' @details
#' IMPORTANT: This function requires parametric imputation methods that produce
#' score functions and information matrices. Supported methods include:
#' - norm (for continuous variables)
#' - logreg (for binary variables)
#' 
#' 
#' @examples
#' \dontrun{
#' library(mice)
#' imp <- mice(nhanes, method = "norm", m = 5, tasks = "train", print = FALSE)
#' fit <- with_rw(imp, lm(bmi ~ age + hyp))
#' pool_rw(fit)
#' }
#' @export
with_rw <- function(data, expr, ...) {
  call <- match.call()
  
  validate_mids_input(data)
  
  m <- data$m
  n <- nrow(data$data)
  
  imputed_vars <- names(data$models)
  model_list <- extract_imputation_models(data, imputed_vars, m)
  
  results <- process_imputed_datasets(data, substitute(expr), model_list, 
                                      imputed_vars, m, n, parent.frame())  
  structure(
    list(
      call = call,
      results = results,
      m = m,
      n = n,
      expr = substitute(expr),
      mids = data
    ),
    class = "with_rw"
  )
}