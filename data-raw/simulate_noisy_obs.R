# Define subroutines to generate design matrix and noisy responses
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

# Define Box and Draper (1987, CH9) quadratic function in 3 vars
fun_3D_quad <- function(x){
  dim(x) <- c(3, 1)
  b <- matrix(c(-1.5, 2.13, -1.81), 3, 1)
  B <- matrix(c(9.38, 7.13, 3.27,
                7.13, 12.54, 2.73,
                3.27, 2.73, 10.42),
              3, 3, TRUE) / 2
  57.31 - t(b) %*% x - t(x) %*% B %*% x
}
true_optima_fun_3D_quad <- c(0.46, -0.46, 0.15)
bounds_3D_quad <- rbind(rep(-2, 3), rep(2, 3))
noise_fun_3D_quad <- 10

# Define a cubic function in 5 vars
fun_5D_cubic <- function(x){
  10 - (x[1] - 1.5)^2 - (x[2] - 2)^2 - (x[3] - 2.5)^2 - (x[4] - 3)^2 - (x[5] - 3.5)^2 + 
    .1 * x[1]^3 - .1 * x[2]^3 - .1 * x[3]^3 - .1 * x[4]^3 - .1 * x[5]^3 + 
    (x[2] - x[3]) * x[4]
}
true_optima_fun_5D_cubic <- c(2.27925, 2.43625, 1.01800, 2.65325, 2.53550)
bounds_5D_cubic <- rbind(rep(0, 5), rep(5, 5))
noise_fun_5D_cubic <- 6

# Simulat noisy obs from Box and Draper quadratic function in 3 vars
set.seed(123)
design_matrix <- sim_design(n = 200, m = 1, bounds_3D_quad)
response <- sim_response(design_matrix, fun_3D_quad, noise_fun_3D_quad * 1)
quad_3D <- list(design_matrix = design_matrix, response = response)
# usethis::use_data(quad_3D, overwrite = TRUE)

# Simulat noisy obs from the cubic function in 5 vars
set.seed(123)
design_matrix <- sim_design(n = 300, m = 15, bounds_5D_cubic)
response <- sim_response(design_matrix, fun_5D_cubic, noise_fun_5D_cubic * 1)
cubic_5D <- list(design_matrix = design_matrix, response = response)
# usethis::use_data(cubic_5D, overwrite = TRUE)