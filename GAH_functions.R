library(GIGrvg)
library(crch)
library(MASS)
library(actuar)

g <- function(gamma) {
  2 * exp(gamma^2 / 2 + pnorm(-abs(gamma), log.p = TRUE))
}

gammaINT <- function(tau0) {
  gammaU <- uniroot(function(x) { g(x) - tau0 }, c(0, 300))$root
  gammaL <- uniroot(function(x) { g(x) - 1 + tau0}, c(-300, 0))$root
  return(c(gammaL, gammaU))
}

f_GAH <- function(y, mu, gamma, eta, rho2, tau0, log.f = FALSE) {
  tau <- as.numeric(I(gamma < 0) + (tau0 - I(gamma < 0)) / g(gamma))
  taugamma1 <- as.numeric(tau - I(gamma > 0)) # plus
  taugamma2 <- as.numeric(tau - I(gamma < 0)) # minus
  func <- function(y, eta, rho2) {
    logP <- function(y, t) {
      gen_lg1 <- function(y) {
        temp <- pnorm(-taugamma1 * y / abs(gamma) + taugamma2 / taugamma1 * abs(gamma), log.p = TRUE) - 
          pnorm(taugamma2 / taugamma1 * abs(gamma), log.p = TRUE)
        index1 <- which(temp <= 600)
        index2 <- which(temp > 600)
        lg1res <- numeric(length(y))
        lg1res[index1] <- pnorm(taugamma2 / taugamma1 * abs(gamma), log.p = TRUE) +
          log(exp(pnorm(-taugamma1 * y[index1] / abs(gamma) + taugamma2 / taugamma1 * abs(gamma), log.p = TRUE) - 
                    pnorm(taugamma2 / taugamma1 * abs(gamma), log.p = TRUE)) - 1)
        lg1res[index2] <- pnorm(-taugamma1 * y[index2] / abs(gamma) + taugamma2 / taugamma1 * abs(gamma), log.p = TRUE)
        return(lg1res)
      }
      ystar <- y / t
      lg1 <- gen_lg1(ystar)
      lg2 <- gamma^2 * taugamma2^2 / taugamma1^2 / 2 - taugamma2 * ystar
      lg3 <- pnorm(-abs(gamma) + taugamma1 * ystar / abs(gamma), log.p = T)
      lg4 <- gamma^2 / 2 - taugamma1 * ystar
      return(lg1 + lg2 + log(1 + exp(lg3 + lg4 - lg1 - lg2)))
    }
    unlist(lapply(y, function(x) { 
      integrate(function(t) { 
        exp(logP(x, t) + dgig(t, 0.5, eta * rho2, eta / rho2, log = T)) },
        0, Inf, stop.on.error = F)$value }))
  }
  
  fy <- numeric(length(y))
  index1 <- which((y - mu) / gamma <= 0)
  index2 <- which((y - mu) / gamma > 0)
  if (gamma == 0) {
    fy <- as.numeric(tau * (1 - tau) * eta / rho2 / (eta + 1) *
                       exp(eta - sqrt(eta * (eta + 2 * (y - mu) * (tau - I(y < mu)) / rho2))))
  } else { 
    fy[index1] <- exp(gamma^2 / 2 + pnorm(-abs(gamma), log.p = TRUE) + eta - sqrt(eta * (eta + 2 * (y[index1] - mu) / rho2 * taugamma1)))
    fy[index2] <- func(y[index2] - mu, eta, rho2)
    fy <- fy * 2 * tau * (1 - tau) * eta / rho2 / (eta + 1)
  }
  if (log.f) { return(log(fy)) } else { return(fy) }
}

logGAH <- function(y, X, beta, gamma, eta, rho2, tau0) {
  tau <- as.numeric(I(gamma < 0) + (tau0 - I(gamma < 0)) / g(gamma))
  taugamma1 <- as.numeric(tau - I(gamma > 0)) # plus
  taugamma2 <- as.numeric(tau - I(gamma < 0)) # minus
  func <- function(y, eta, rho2) {
    logP <- function(y, t) {
      gen_lg1 <- function(y) {
        temp <- pnorm(-taugamma1 * y / abs(gamma) + taugamma2 / taugamma1 * abs(gamma), log.p = TRUE) - 
          pnorm(taugamma2 / taugamma1 * abs(gamma), log.p = TRUE)
        index1 <- which(temp <= 600)
        index2 <- which(temp > 600)
        lg1res <- numeric(length(y))
        lg1res[index1] <- pnorm(taugamma2 / taugamma1 * abs(gamma), log.p = TRUE) +
          log(exp(pnorm(-taugamma1 * y[index1] / abs(gamma) + taugamma2 / taugamma1 * abs(gamma), log.p = TRUE) - 
                    pnorm(taugamma2 / taugamma1 * abs(gamma), log.p = TRUE)) - 1)
        lg1res[index2] <- pnorm(-taugamma1 * y[index2] / abs(gamma) + taugamma2 / taugamma1 * abs(gamma), log.p = TRUE)
        return(lg1res)
      }
      ystar <- y / t
      lg1 <- gen_lg1(ystar)
      lg2 <- gamma^2 * taugamma2^2 / taugamma1^2 / 2 - taugamma2 * ystar
      lg3 <- pnorm(-abs(gamma) + taugamma1 * ystar / abs(gamma), log.p = T)
      lg4 <- gamma^2 / 2 - taugamma1 * ystar
      return(lg1 + lg2 + log(1 + exp(lg3 + lg4 - lg1 - lg2)))
    }
    unlist(lapply(y, function(x) { 
      integrate(function(t) { 
        exp(logP(x, t) + dgig(t, 0.5, eta * rho2, eta / rho2, log = T)) },
        0, Inf, stop.on.error = F)$value }))
  }
  
  ynew <- y - X %*% beta
  fynew <- numeric(length(ynew))
  index1 <- which(ynew / gamma <= 0)
  index2 <- which(ynew / gamma > 0)
  if (gamma == 0) {
    fynew <- as.numeric(tau * (1 - tau) * eta / rho2 / (eta + 1) *
                          exp(eta - sqrt(eta * (eta + 2 * ynew * (tau - I(ynew < 0)) / rho2))))
  } else { 
    fynew[index1] <- exp(gamma^2 / 2 + pnorm(-abs(gamma), log.p = TRUE) + eta - sqrt(eta * (eta + 2 * ynew[index1] / rho2 * taugamma1)))
    fynew[index2] <- func(ynew[index2], eta, rho2)
    fynew <- fynew * 2 * tau * (1 - tau) * eta / rho2 / (eta + 1)
  }
  return(log(fynew))
}

BQR_GAH <- function(y, X, tau0, nburn, niter, thin, method = "alasso",
                    beta0, gamma0 = 0, eta0 = 1, rho20 = 1,
                    a = 1, b = 1, r = 1, delta = 1,
                    model.fit = FALSE) {
  # y: response variable
  # X: predictor variable (including an intercept term)
  # tau0: expected quantile level
  # method: "alasso" or "lasso"
  # model.fit: whether return DIC and LPML
  
  logit <- function(x, gamINT) {
    return(log((x - gamINT[1]) / (gamINT[2] - x)))
  }
  inv.logit <- function(x, gamINT) {
    return((gamINT[2] * exp(x) + gamINT[1]) / (exp(x) + 1))
  }
  
  Gibbs_beta <- function(y, X, s, z, sigma, gamma, rho2, v) {
    n <- nrow(X); p <- ncol(X)
    tau <- as.numeric(I(gamma < 0) + (tau0 - I(gamma < 0)) / g(gamma))
    A <- (1 - 2 * tau) / tau / (1 - tau)
    B <- 2 / tau / (1 - tau)
    C <- as.numeric(1 / (I(gamma > 0) - tau))
    V_inv <- diag(1 / B / sigma / z, n, n)
    v_inv <- diag(1 / v, p, p)
    Sigma <- solve(t(X) %*% V_inv %*% X + v_inv / rho2)
    mu <- Sigma %*% t(X) %*% V_inv %*% (y - C * abs(gamma) * sigma * s - A * z)
    beta_new <- mvrnorm(1, mu, Sigma)
    return(beta_new)
  }
  Gibbs_s <- function(y, X, beta, z, sigma, gamma) {
    n <- nrow(X); p <- ncol(X)
    tau <- as.numeric(I(gamma < 0) + (tau0 - I(gamma < 0)) / g(gamma))
    A <- (1 - 2 * tau) / tau / (1 - tau)
    B <- 2 / tau / (1 - tau)
    C <- as.numeric(1 / (I(gamma > 0) - tau))
    mus <- C * abs(gamma) * (y - X %*% beta - A * z) / (B * z + C^2 * gamma^2 * sigma)
    sigmas <- sqrt(B * z / (B * z + C^2 * gamma^2 * sigma))
    s_new <- rtnorm(n, mus, sigmas, 0, Inf)
    s_new[s_new == Inf] <- 0
    return(s_new)
  }
  Gibbs_z <- function(y, X, beta, s, sigma, gamma) {
    n <- nrow(X); p <- ncol(X)
    tau <- as.numeric(I(gamma < 0) + (tau0 - I(gamma < 0)) / g(gamma))
    A <- (1 - 2 * tau) / tau / (1 - tau)
    B <- 2 / tau / (1 - tau)
    C <- as.numeric(1 / (I(gamma > 0) - tau))
    psi <- (A^2 + 2 * B) / B / sigma
    chi <- (y - X %*% beta - C * abs(gamma) * sigma * s)^2 / B / sigma
    if (length(psi) == n) {
      z_new <- unlist(lapply(1:n, function(i) rgig(1, 1/2, chi[i], psi[i])))
    } else {
      z_new <- unlist(lapply(1:n, function(i) rgig(1, 1/2, chi[i], psi)))
    }
    return(z_new)
  }
  Gibbs_sigma <- function(y, X, beta, s, z, gamma, eta, rho2) {
    n <- nrow(X); p <- ncol(X)
    tau <- as.numeric(I(gamma < 0) + (tau0 - I(gamma < 0)) / g(gamma))
    A <- (1 - 2 * tau) / tau / (1 - tau)
    B <- 2 / tau / (1 - tau)
    C <- as.numeric(1 / (I(gamma > 0) - tau))
    psi <- C^2 * gamma^2 * s^2 / B / z + eta / rho2
    chi <- (y - X %*% beta - A * z)^2 / B / z + 2 * z + eta * rho2
    sigma_new <- unlist(lapply(1:n, function(i) rgig(1, 0, chi[i], psi[i])))
    return(sigma_new)
  }
  RWMH_gamma <- function(y, X, beta, s, z, sigma, gamma, suggest) {
    logfgamma <- function(y, X, beta, s, z, sigma, gamma) {
      tau <- as.numeric(I(gamma < 0) + (tau0 - I(gamma < 0)) / g(gamma))
      A <- (1 - 2 * tau) / tau / (1 - tau)
      B <- 2 / tau / (1 - tau)
      C <- as.numeric(1 / (I(gamma > 0) - tau))
      sum(dnorm(y, mean = X %*% beta + C * abs(gamma) * sigma * s + A * z,
                sd = sqrt(B * sigma * z), log = TRUE))
    }
    n <- nrow(X); p <- ncol(X)
    theta_now <- logit(gamma, gamINT)
    theta_star <- rnorm(1, theta_now, suggest)
    gamma_star <- inv.logit(theta_star, gamINT)
    propose <- logfgamma(y, X, beta, s,  z, sigma, gamma_star) + 
      log(gamINT[2] - gamINT[1]) + theta_star - 2 * log(1 + exp(theta_star))
    origin <- logfgamma(y, X, beta, s, z, sigma, gamma)  +
      log(gamINT[2] - gamINT[1]) + theta_now - 2 * log(1 + exp(theta_now))
    log_acc_prob <- min(0, propose - origin)
    accept <- log(runif(1)) < log_acc_prob
    new_gamma <- accept * gamma_star + (1 - accept) * gamma
    return(new_gamma)
  }
  Gibbs_v <- function(beta, rho2, lambda, method = "alasso") {
    p <- length(beta)
    psi <- lambda^2
    chi <- beta^2 / rho2
    if (method == "alasso") {
      v_new <- unlist(lapply(1:p, function(i) rgig(1, 1/2, chi[i], psi[i])))
    } else {
      v_new <- unlist(lapply(1:p, function(i) rgig(1, 1/2, chi[i], psi)))
    }
    return(v_new)
  }
  Gibbs_rho2 <- function(beta, sigma, eta, v) {
    n <- length(sigma); p <- length(beta)
    v_inv <- diag(1 / v, p, p)
    psi <- eta * sum(1 / sigma)
    chi <- eta * sum(sigma) + t(beta) %*% v_inv %*% beta
    rho2_new <- rgig(1, -(3 * n + p) / 2, chi, psi)
    return(rho2_new)
  }
  RWMH_eta <- function(sigma, rho2, a, b, eta, suggest) {
    n <- length(sigma)
    P <- sum(sigma / rho2 + rho2 / sigma) / 2 + b - n
    eta_star <- exp(rnorm(1, log(eta), suggest))
    propose <- (a + 3 * n / 2) * log(eta_star) - n * log(eta_star + 1) - eta_star * P
    origin <- (a + 3 * n / 2) * log(eta) - n * log(eta + 1) - eta * P
    log_acc_prob <- min(0, propose - origin)
    accept <- log(runif(1)) < log_acc_prob
    eta_new <- accept * eta_star + (1 - accept) * eta
    return(eta_new)
  }
  Gibbs_lambda <- function(v, r, delta, method = "alasso") {
    if (method == "alasso") {
      p <- length(v)
      lambda_new <- sqrt(rgamma(p, r + 1, delta + v / 2))
    } else {
      lambda_new <- sqrt(rgamma(1, r + p, delta + sum(v / 2)))
    }
    return(lambda_new)
  }
  
  n <- nrow(X); p <- ncol(X); gamINT <- gammaINT(tau0)
  beta <- matrix(0, niter, p); beta[1,] <- ifelse(missing(beta0), 0, beta0)
  gamma <- rep(0, niter); gamma[1] <- gamma0
  rho2 <- rep(0, niter); rho2[1] <- rho20
  eta <- rep(0, niter); eta[1] <- eta0
  s <- rep(1, n)
  z <- rep(1, n)
  sigma <- rep(1, n)
  v <- rep(1, p)
  if (method == "alasso") { 
    lambda <- rep(1, p)
  } else {
    if (method == "lasso") {
      lambda <- 1
    } else {
      print("Method not supported!")
      return()
    }
  }
  accept_gamma <- accept_eta <- rep(0, nburn)
  scf_gamma <- scf_eta <- rep(2.38, niter); DD <- 0.0238
  
  for (brn in 2:nburn) {
    beta[brn,] <- Gibbs_beta(y, X, s, z, sigma, gamma[brn-1], rho2[brn-1], v)
    s <- Gibbs_s(y, X, beta[brn,], z, sigma, gamma[brn-1])
    z <- Gibbs_z(y, X, beta[brn,], s, sigma, gamma[brn-1])
    sigma <- Gibbs_sigma(y, X, beta[brn,], s, z, gamma[brn-1], eta[brn-1], rho2[brn-1])
    
    suggest <- 0.1
    if (sum(accept_gamma) > 10 & runif(1) < 0.95) {
      suggest <- scf_gamma[brn-1] * sd(logit(gamma[1:(brn-1)], gamINT))
    }
    gamma[brn] <- RWMH_gamma(y, X, beta[brn,], s, z, sigma, gamma[brn-1], suggest)
    accept_gamma[brn] <- gamma[brn] != gamma[brn-1]
    scf_gamma[brn] <- scf_gamma[brn-1] + (3.3 * accept_gamma[brn] - 1) * DD / sqrt(brn) * (suggest != 0.1)
    
    suggest <- 0.1
    if (sum(accept_eta) > 10 & runif(1) < 0.95) {
      suggest <- scf_eta[brn-1] * sd(log(eta[1:(brn-1)]))
    }
    eta[brn] <- RWMH_eta(sigma, rho2[brn-1], a, b, eta[brn-1], suggest)
    accept_eta[brn] <- eta[brn] != eta[brn-1]
    scf_eta[brn] <- scf_eta[brn-1] + (3.3 * accept_eta[brn] - 1) * DD / sqrt(brn) * (suggest != 0.1)
    
    rho2[brn] <- Gibbs_rho2(beta[brn,], sigma, eta[brn], v)
    v <- Gibbs_v(beta[brn,], rho2[brn], lambda, method)
    lambda <- Gibbs_lambda(v, r, delta, method)
    
    if (brn %% 1000 == 0) { print(paste0("==== burn ", brn, " done! ====")) }
  }
  
  beta_brn <- beta[1:nburn,]
  gamma_brn <- gamma[1:nburn]
  eta_brn <- eta[1:nburn]
  rho2_brn <- rho2[1:nburn]
  
  beta[1, ] <- beta[nburn, ]
  gamma[1] <- gamma[nburn]
  eta_brn[1] <- eta[nburn]
  rho2_brn[1] <- rho2[nburn]
  accept_gamma <- accept_eta <- rep(0, niter)
  scf_gamma[1] <- scf_gamma[nburn]
  scf_eta[1] <- scf_eta[nburn]
  
  for (itr in 2:niter) {
    beta[itr,] <- Gibbs_beta(y, X, s, z, sigma, gamma[itr-1], rho2[itr-1], v)
    s <- Gibbs_s(y, X, beta[itr,], z, sigma, gamma[itr-1])
    z <- Gibbs_z(y, X, beta[itr,], s, sigma, gamma[itr-1])
    sigma <- Gibbs_sigma(y, X, beta[itr,], s, z, gamma[itr-1], eta[itr-1], rho2[itr-1])
    
    suggest <- 0.1
    if (sum(accept_gamma) > 10 & runif(1) < 0.95) {
      suggest <- scf_gamma[itr-1] * sd(logit(gamma[1:(itr-1)], gamINT))
    }
    gamma[itr] <- RWMH_gamma(y, X, beta[itr,], s, z, sigma, gamma[itr-1], suggest)
    accept_gamma[itr] <- gamma[itr] != gamma[itr-1]
    scf_gamma[itr] <- scf_gamma[itr-1] + (3.3 * accept_gamma[itr] - 1) * DD / sqrt(itr) * (suggest != 0.1)
    
    suggest <- 0.1
    if (sum(accept_eta) > 10 & runif(1) < 0.95) {
      suggest <- scf_eta[itr-1] * sd(log(eta[1:(itr-1)]))
    }
    eta[itr] <- RWMH_eta(sigma, rho2[itr-1], a, b, eta[itr-1], suggest)
    accept_eta[itr] <- eta[itr] != eta[itr-1]
    scf_eta[itr] <- scf_eta[itr-1] + (3.3 * accept_eta[itr] - 1) * DD / sqrt(itr) * (suggest != 0.1)
    
    rho2[itr] <- Gibbs_rho2(beta[itr,], sigma, eta[itr], v)
    v <- Gibbs_v(beta[itr,], rho2[itr], lambda, method)
    lambda <- Gibbs_lambda(v, r, delta, method)
    
    if (itr %% 5000 == 0) { print(paste0("==== iter ", itr, " done! ====")) } 
  }
  
  if (!model.fit) {
    sample.ind <- 1:(niter / thin) * thin
    return(list(beta = beta[sample.ind,], gamma = gamma[sample.ind], eta = eta[sample.ind], rho2 = rho2[sample.ind],
                beta_brn = beta_brn, gamma_brn = gamma_brn, eta_brn = eta_brn, rho2_brn = rho2_brn))
  } else {
    sample.ind <- 1:(niter / thin) * thin
    result_BQR_GAH <- cbind(beta[sample.ind,], gamma[sample.ind], eta[sample.ind], rho2[sample.ind])
    med_result_BQR_GAH <- apply(result_BQR_GAH, 2, median)
    avg_result_BQR_GAH <- colMeans(result_BQR_GAH)
    logLs_GAH <- cbind(apply(result_BQR_GAH, 1, function(res) { logGAH(y, X, res[1:p], res[p+1], res[p+2], res[p+3], tau0) }),
                       logGAH(y, X, avg_result_BQR_GAH[1:p], avg_result_BQR_GAH[p+1], avg_result_BQR_GAH[p+2], avg_result_BQR_GAH[p+3], tau0))
    DIC <- 2 * sum(logLs_GAH[, length(sample.ind) + 1]) - 4 * sum(rowMeans(logLs_GAH[, 1:length(sample.ind)]))
    LPML <- -sum(log(rowMeans(exp(-logLs_GAH[, 1:length(sample.ind)]))))
    return(list(beta = beta[sample.ind,], gamma = gamma[sample.ind], eta = eta[sample.ind], rho2 = rho2[sample.ind],
                beta_brn = beta_brn, gamma_brn = gamma_brn, eta_brn = eta_brn, rho2_brn = rho2_brn,
                beta_hat = med_result_BQR_GAH[1:p], gamma_hat = med_result_BQR_GAH[p+1], eta_hat = med_result_BQR_GAH[p+2], rho2_hat = med_result_BQR_GAH[p+3],
                DIC = DIC, LPML = LPML))
  }
}