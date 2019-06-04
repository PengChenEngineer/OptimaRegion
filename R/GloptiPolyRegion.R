#' Compute confidence regions for polynomial models' optima
#' 
#' @param X numeric matrix of shape (N, k); N is the sample size; k is the 
#'          number of variables, which can be 2, 3, 4 and 5; X specifies the 
#'          design matrix
#' @param y numeric vector of shape (N, 1); y specifies the responses
#' @param degree integer scalor; degree specifies the order of the polynomial 
#'               model, which can be 2 or 3
#' @param lb numeric vector of shape (1, k); lb specifies the lower bounds for
#'           the k variables
#' @param ub numeric vector of shape (1, k); ub specifies the upper bounds for 
#'           the k variables 
#' @param B integer scalor; B specifies the number of bootstrap operations
#' @param alpha numeric scalor between 0 and 1; alpha specifies the nominal 
#'              confidence level, 1 - alpha, of the confidence region
#' @param maximization boolean scalor; if specifies whether the algorithm 
#'                     computes the confidence region for the maxima or minima
#' @param verbose boolean scalor; it specifies whether to display running status
#' @return a figure displaying the confidence region of the true optimum,
#'         projected onto each pairwise-variable planes and a list consisting
#'         of 2 components: 
#'         \describe{
#'           \item{boot_optima}{numeric matrix of shape ((1 - alpha)B, k);
#'                                it contains the (1 - alpha)B bootstrap optima}
#'           \item{boost_optimum}{numeric vector of shape (1, k); the boosted 
#'                                  optimum; computed by taking the column average
#'                                  of boot_optima}
#'         }
#' @inheritSection OptRegionQuad Author(s)
#' @examples 
#' \dontrun{
#' # Example 1: run GloptiPolyRegion on a quadratic, 3 vars example
#' out <- GloptiPolyRegion(X = quad_3D[, 1:3], y = quad_3D[, 4], degree = 2, 
#'                         lb = c(-2, -2, -2), ub = c(2, 2, 2), B = 2000, alpha = 0.1, 
#'                         maximization = TRUE, verbose = TRUE)
#' # check result
#' str(out)
#' 
#' # define subroutines to draw 3D confidence reigon
#' library(rgl)
#' plot_3D_CR_demo <- function(X){
#'   plot3d(X, col = "green",
#'          type = "p", size = 5, alpha = 0.01,
#'          xlab = "x1", ylab = "x2", zlab = "x3",
#'          xlim = c(-2, 2), ylim = c(-2, 2), zlim = c(-2, 2))
#'   crownhull(X, col = "green", alpha = 0.2)
#'   points3d(X, add = TRUE, col = "green",
#'            size = 2, alpha = 0.5)
#'   X_ave <- apply(X, 2, mean)
#'   points3d(X_ave[1],
#'            X_ave[2],
#'            X_ave[3],
#'            add = TRUE, col = "red",
#'            size = 10)
#' }
#' crownhull <- function(xyz, plotit = TRUE, col = "green", alpha = 0.8){
#' if(is.list(xyz) && !is.data.frame(xyz))
#'   p <- as.matrix(do.call("rbind", xyz))
#' else
#'   p <- as.matrix(xyz)
#'  ch <- geometry::convhulln(p, "FA")
#' if(plotit){
#'   ch2 <- t(geometry::convhulln(p, "Qt"))
#'   triangles3d(p[ch2,1], p[ch2,2], p[ch2,3], col = col, alpha = alpha, add = TRUE)
#' }
#' return(list(crownvolume = ch$vol, crownsurface = ch$area))
#' }
#'  
#' # draw 3D confidence region based on out
#' plot_3D_CR_demo(out$boot_optima)
#' 
#' # Example 2: run GloptiPolyRegion on a cubic, 5 vars example
#' out <- GloptiPolyRegion(X = cubic_5D$design_matrix, y = cubic_5D$response,  
#'                         degree = 3, lb = rep(0, 5), ub = rep(5, 5), B = 2000,  
#'                         alpha = 0.05, maximization = TRUE, verbose = TRUE)
#' # check result
#' str(out)
#' }
#' @export                          
GloptiPolyRegion <- function(X, y, degree, lb, ub, B = 200, alpha = 0.05,
                            maximization = TRUE, verbose = TRUE){
  # Check polynomial order -- -----------------------------------------------
  if(degree < 2 || degree > 3){
    stop("This function accepts only quadratic or cubic polynomials!")
  }
  # Check number of variables -----------------------------------------------
  if(ncol(X) < 2 || ncol(X) > 5){
    stop("This function accepts only 2 - 5 variables!")
  }
  # Simplify function arguments ---------------------------------------------
  plot_CR = TRUE # draw pairwise projected CR's 
  scale = TRUE # scale X to [-1, 1]
  # Original fit ------------------------------------------------------------
  if(verbose) print("Fit original model ...")
  if(scale) X <- encode_orthogonal(X, lb, ub) # scale X to [-1, 1]
  model_original <- fit_polym(X, y, degree) # fit polynomial model
  theta_hat <- model_original$coefficients; cov_theta_hat <- vcov(model_original)
  p <- length(theta_hat); n <- length(y)
  # Bootstrap ---------------------------------------------------------------
  if(verbose) print("Bootstrap ...")
  theta_hat_star <- matrix(NA, nrow = B, ncol = p)
  theta_hat_star_studentized <- matrix(NA, nrow = B, ncol = p)
  boot_optima <- matrix(NA, nrow = B, ncol = ncol(X))
  fit <- fitted(model_original); names(fit) <- NULL
  res <- resid(model_original) - mean(resid(model_original))
  # accumulate B "optimizable" bootstrap instances
  counter <- 0
  while(counter < B){
    model_star <- fit_polym(X, fit + res[sample(1:n, replace = TRUE)], degree)
    opt <- try(
      # optimize a polynomial model fitted through the polym formula using 
      # the GloptiPolyR solver
      Gloptipolym(coefficients = model_star$coefficients,
                  k = ncol(X), degree = degree,
                  lb = rep(-1, ncol(X)), ub = rep(1, ncol(X)),
                  maximization = maximization)
    )
    if(!inherits(opt, "try-error")){
      counter <- counter + 1
      if(verbose) print(paste(counter, "th bootstrapping ..."))
      theta_hat_star[counter, ] <- model_star$coefficients
      e <- eigen(vcov(model_star))
      S_sqrt_inv <- solve(e$vectors %*% diag(sqrt(e$values)) %*% t(e$vectors))
      theta_hat_star_studentized[counter, ] <- sqrt(n) * S_sqrt_inv %*%
        (theta_hat_star[counter, ] - theta_hat)
      boot_optima[counter, ] <- opt$solution
    }
  }
  # Trim -------------------------------------------------------------------
  if(verbose) print("Trimming ...")
  d <- vector(length = B)
  d <- DepthProc::depthTukey(theta_hat_star_studentized,
                             theta_hat_star_studentized, ndir = 3000)
  order_d <- order(d); ind_alpha <- alpha * B + 1
  indices <- order_d[ind_alpha:B]
  # Optimization ------------------------------------------------------------
  # extract already optimized results that remains after trimming
  if(verbose) print("Optimizing over bootstrapped models that remains...")
  boot_optima <- boot_optima[indices, ]
  if(scale) boot_optima <- decode_orthogonal(boot_optima, lb, ub)
  boost_optimum <- apply(boot_optima, 2, mean)
  # plotting ----------------------------------------------------------------
  if(plot_CR){
    if(verbose) print("Ploting the confidence region ... ")
    draw_2D_CRs(boot_optima, boost_optimum, lb, ub)
  }
  # return ------------------------------------------------------------------
  list(boot_optima = boot_optima, boost_optimum = boost_optimum)
}