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
  
  # Input validation
  if (!inherits(data, "mids")) {
    cli::cli_abort(
      "{.arg data} must be a {.cls mids} object from {.pkg mice}, not {.cls {class(data)}}."
    )
  }
  
  if (is.null(data$models)) {
    cli::cli_abort(c(
      "{.arg data} must contain imputation models.",
      "i" = "Run {.fn mice} with {.code tasks = 'train'} to save models."
    ))
  }
  
  # Check that all imputation methods are supported (parametric)
  supported_methods <- c("norm", "logreg")
  
  imputed_vars <- names(data$models)
  methods_used <- data$method[imputed_vars]
  
  unsupported <- setdiff(methods_used, supported_methods)
  if (length(unsupported) > 0) {
    cli::cli_abort(c(
      "Unsupported imputation method{?s}: {.val {unsupported}}",
      "i" = "Robins-Wang variance currently only supports: {.val {supported_methods}}",
      "x" = "Other methods may not store the required variance information."
    ))
  }
  
  m <- data$m
  n <- nrow(data$data)
  
  # Extract imputation model info for each iteration
  model_list <- vector("list", m)
  for (p in seq_len(m)) {
    model_list[[p]] <- lapply(imputed_vars, function(var) {
      mod <- data$models[[var]][[p]]
      
      # Extract coefficients
      beta <- stats::setNames(as.numeric(mod$beta.dot), mod$xnames)
      
      # Determine family from imputation method
      imp_method <- data$method[var]
      
      if (imp_method == "logreg") {
        list(family = "binomial", coefficients = beta)
      } else if (imp_method == "norm") {
        # Extract sigma2
        sigma2 <- if (!is.null(mod$sigma2)) {
          as.numeric(mod$sigma2)
        } else if (!is.null(mod$sigma.dot)) {
          as.numeric(mod$sigma.dot)^2
        } else if (!is.null(mod$sigma)) {
          as.numeric(mod$sigma)^2
        } else {
          cli::cli_abort(c(
            "Cannot extract {.field sigma2} from imputation model for {.field {var}}.",
            "i" = "Model object is missing {.field sigma2}, {.field sigma.dot}, and {.field sigma}."
          ))
        }
        list(family = "gaussian", coefficients = beta, sigma2 = sigma2)
      } else {
        cli::cli_abort(
          "Unsupported imputation method {.val {imp_method}} for variable {.field {var}}."
        )
      }
    })
    names(model_list[[p]]) <- imputed_vars
  }
  
  # Process each imputed dataset
  results <- vector("list", m)
  
  for (p in seq_len(m)) {
    data.i <- mice::complete(data, action = p)
    
    # Convert factors back to numeric 0/1 for logreg variables
    for (var in imputed_vars) {
      if (is.factor(data.i[[var]])) {
        # logreg returns factors - convert to numeric 0/1
        levels_char <- levels(data.i[[var]])
        if (all(levels_char %in% c("0", "1"))) {
          # Levels are "0" and "1" - convert via character
          data.i[[var]] <- as.numeric(as.character(data.i[[var]]))
        } else {
          # Generic binary: first level = 0, second = 1
          data.i[[var]] <- as.numeric(data.i[[var]]) - 1
        }
      }
    }
    
    # Track which observations were imputed
    imputed_mask <- Reduce(`|`, lapply(imputed_vars, function(var) {
      as.integer(is.na(data$data[[var]]))
    }), init = rep(0, n))
    data.i$.imputed <- imputed_mask
    
    # Fit analysis model
    mod_analysis <- eval(expr = substitute(expr), envir = data.i, enclos = parent.frame())
    if (is.expression(mod_analysis)) {
      mod_analysis <- eval(expr = mod_analysis, envir = data.i, enclos = parent.frame())
    }
    
    if (!inherits(mod_analysis, c("lm", "glm"))) {
      cli::cli_abort(c(
        "Analysis model must be {.cls lm} or {.cls glm}, not {.cls {class(mod_analysis)}}.",
        "i" = "Use {.fn lm} or {.fn glm} for analysis models."
      ))
    }
    
    # Compute imputation components (S_mis_imp and d)
    S_u_list <- lapply(seq_along(imputed_vars), function(i) {
      compute_score(data.i, model_list[[p]][[i]], imputed_vars[i])
    })
    S_u <- do.call(cbind, S_u_list)
    
    ImputedMat <- matrix(data.i$.imputed == 1, nrow(S_u), ncol(S_u), byrow = FALSE)
    S_mis_imp <- S_u * ImputedMat
    S_orig <- S_u * (1 - ImputedMat)
    
    S2_list <- lapply(seq_along(imputed_vars), function(i) {
      compute_information(data.i, model_list[[p]][[i]], imputed_vars[i], data.i$.imputed)
    })
    
    # Build block diagonal S2 matrix
    dims <- vapply(S2_list, ncol, integer(1))
    n_tot <- sum(dims)
    S2 <- matrix(0, n_tot, n_tot)
    idx <- 1
    for (i in seq_along(S2_list)) {
      rng <- idx:(idx + dims[i] - 1)
      S2[rng, rng] <- S2_list[[i]]
      idx <- idx + dims[i]
    }
    
    Dmat <- solve(S2)
    d <- t((-1) * Dmat %*% t(S_orig))
    
    # Compute analysis components (U and tau)
    X_analysis <- stats::model.matrix(mod_analysis)
    
    if (inherits(mod_analysis, "lm") && !inherits(mod_analysis, "glm")) {
      y_analysis <- stats::model.response(stats::model.frame(mod_analysis))
      pred <- stats::fitted(mod_analysis)
      resid <- y_analysis - pred
      U_analysis <- X_analysis * resid
      
      tau <- - crossprod(X_analysis)
    } else {
      # GLM
      y_analysis <- stats::model.response(stats::model.frame(mod_analysis))
      pred <- stats::fitted(mod_analysis)
      U_analysis <- X_analysis * (y_analysis - pred)
      
      mu <- stats::predict(mod_analysis, type = "response")
      w <- stats::family(mod_analysis)$mu.eta(
        stats::predict(mod_analysis, type = "link"))^2 /
        stats::family(mod_analysis)$variance(mu)
      tau <- - crossprod(X_analysis * sqrt(w))
    }
    
    # Extend U back to full dataset
    U_imp <- matrix(0, n, ncol(U_analysis))
    colnames(U_imp) <- names(stats::coef(mod_analysis))
    
    mf <- stats::model.frame(mod_analysis)
    analysis_ids <- as.numeric(rownames(mf))
    if (is.null(analysis_ids) || any(is.na(analysis_ids))) {
      analysis_ids <- match(rownames(mf), rownames(data.i))
    }
    
    U_imp[analysis_ids, ] <- U_analysis

    results[[p]] <- list(
      model = mod_analysis,
      U = U_imp,
      tau = tau,
      S_mis_imp = S_mis_imp,
      d = d,
      n_analysis = nrow(mf)
    )
  }
  
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
