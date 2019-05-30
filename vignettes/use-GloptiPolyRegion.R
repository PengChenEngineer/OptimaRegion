## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup---------------------------------------------------------------
library(OptimaRegion)

## ------------------------------------------------------------------------
# Box and Draper Benmark Problem (1987, CH9)
fun_3D_quad <- function(x){
  dim(x) <- c(3, 1)
  b <- matrix(c(-1.5, 2.13, -1.81), 3, 1)
  B <- matrix(c(9.38, 7.13, 3.27,
                7.13, 12.54, 2.73,
                3.27, 2.73, 10.42),
              3, 3, TRUE) / 2
  57.31 - t(b) %*% x - t(x) %*% B %*% x
}
bounds_3D_quad <- rbind(rep(-2, 3), rep(2, 3))
noise_fun_3D_quad <- 10

## ------------------------------------------------------------------------
sim_design <- function(n, m = 1, bounds, method = "LHD"){
  design_matrix <- matrix(NA, nrow = n, ncol = ncol(bounds))
  if(method == "random"){
    for(i in 1:ncol(bounds)){
      design_matrix[, i] <- runif(n, min = bounds[1, i], max = bounds[2, i])
    }
  }else if(method == "LHD"){
    design_matrix <- lhs::geneticLHS(n = n, k = ncol(bounds), criterium = "Maximin")
    design_matrix <- sweep(design_matrix, 2, bounds[2, ] - bounds[1, ], "*")
    design_matrix <- sweep(design_matrix, 2, bounds[1, ], "+")
  }
  # replicates
  design_matrix[rep(1:n, rep(m, n)), ]
}
sim_response <- function(design_matrix, fun, sigma_noise){
  N <- nrow(design_matrix) # num_locations * num_replicates
  # f <- matrix(NA, nrow = N, ncol = 1)
  y <- matrix(NA, nrow = N, ncol = 1)
  set.seed(as.numeric(Sys.time()))
  for(i in 1:N) y[i] <- fun(design_matrix[i, ]) + rnorm(1, 0, sigma_noise)
  # f + MASS::mvrnorm(n = 1, mu = rep(0, N), Sigma = sigma_noise * diag(N))
  y
}

## ---- fig.width = 7, fig.height = 7, fig.align = "center", warning = FALSE----
library(magrittr)
set.seed(123)
design_matrix <- sim_design(n = 200, m = 1, bounds_3D_quad)
response <- sim_response(design_matrix, fun_3D_quad, noise_fun_3D_quad * 0)
system.time(
  CRO_3D_quad <- GloptiPolyRegion(X = design_matrix,
                                  y = response,
                                  degree = 2,
                                  lb = bounds_3D_quad[1, ],
                                  ub = bounds_3D_quad[2, ],
                                  B = 200,
                                  verbose = FALSE)
)
CRO_3D_quad$boost_optimum

