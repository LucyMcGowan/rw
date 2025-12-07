compute_omega <- function(u_bar, n) {
  crossprod(u_bar) / n
}

compute_variance_components <- function(results, m, n) {
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
  
  list(
    kappa = kappa_sum / (n * m),
    alpha = alpha_sum / (n * m),
    d_bar = d_bar_sum / m
  )
}

compute_rw_variance <- function(results, m, n) {
  U_sum <- Reduce(`+`, lapply(results, function(r) r$U))
  tau_sum <- Reduce(`+`, lapply(results, function(r) r$tau))
  
  u_bar <- U_sum / m
  omega <- compute_omega(u_bar, n)
  
  components <- compute_variance_components(results, m, n)
  
  delta <- omega + 
    components$kappa %*% components$alpha %*% t(components$kappa) +
    (1 / n) * (components$kappa %*% t(components$d_bar) %*% u_bar + 
                 t(components$kappa %*% t(components$d_bar) %*% u_bar))
  
  tau <- tau_sum / (m * n)
  tau_inv <- solve(tau)
  (1 / n) * tau_inv %*% delta %*% t(tau_inv)
}

pool_estimates <- function(results) {
  coefs <- vapply(results, function(r) stats::coef(r$model), 
                  numeric(length(stats::coef(results[[1]]$model))))
  rowMeans(coefs)
}

construct_pooled_output <- function(est, se, df, family) {
  test_stat <- est / se
  
  if (family == "gaussian") {
    crit_val <- stats::qt(0.975, df = df)
    p_value <- 2 * stats::pt(-abs(test_stat), df = df)
  } else {
    crit_val <- stats::qnorm(0.975)  
    p_value <- 2 * stats::pnorm(-abs(test_stat))
  }
  out <- data.frame(
    term = names(est),
    estimate = est,
    std.error = se,
    statistic = test_stat,
    p.value = p_value,
    conf.low = est - crit_val * se,
    conf.high = est + crit_val * se
  )
  out
}