#' Confidence Regions for Response Surface Optima
#'
#' This function bootstraps a Confidence Region (CR) on the resposne surface optima.
#' One can choose from several popular Response Surface Models (RSM) to capture the true
#' reponse pattern, such as polynomial model and Thin Plate Spline (TPS) model
#' \insertCite{DelCastilloCR}{OptimaRegion}.
#' 
#'@param X A numeric matrix of shape (N, k), where N represents the number of data points
#'         and k represents the number of regressors. 
#'         It contains the N design points from the experimental data.
#'         The value that k can take depends on the RSM argument.
#'@param y A numeric vector of shape (N, 1), where N represents the number of data points.
#'         It contains the N response values from the experimental data.
#'@param alpha A numeric scalor that specifies the confidence level, 1 - alpha, of the CR.
#'             Its value must be set between 0 and 1.
#'             Default is 0.05.
#'@param num_boot An integer scalor that specifies how many times the bootstrap operation
#'                will be repeated to construct the CR.
#'                Default is 200.
#'@param RSM A character string to choose the type of the RSM, 
#'            which must be either "quad", "tps", or "poly".
#'            "quad" fits the experimental data to a quadratic polynomial model
#'            with 2 regressors;
#'            "tps" fits the experimental data to a TPS model with 2 regressors;
#'            "poly" fits the experimental data to a quadratic or cubic polynomial model 
#'            (specified by the degree argument) with 2 - 5 regressors.
#'            For all model types, 
#'            the optima of the response surface will be constrained by the bounds of the
#'            regressors (specified by the constr_lb and constr_ub arguments).
#'            For "quad" and "tps" models, 
#'            the optima of the response surface can be further constrained to lie inside 
#'            a triangle defined by the original (0,0) and two other vertices
#'            (specified by the constr_triangle, constr_vertex_1, and constr_vertex_2 
#'            arguments).
#'@param lambda A numeric scalor that specifies the value of the penalty parameter 
#'              if RSM is "tps". 
#'              Default is 0.04.
#'@param degree A integer scalor that specifies the degree of the polynomial model
#'              if RSM is "poly". 
#'              It can be set to 2 or 3.
#'              Default is 2.
#'@param maximization A boolean scalor.
#'                    If TRUE, the function returns a CR on the maxima of the response surface;
#'                    if FALSE, the function returns a CR on the minima of the response surface.
#'                    Default is TRUE.
#'@param constr_lb A numeric vector of shape (1, k) that specifies the lower bound constraint
#'                 for each of the k regressors.
#'@param constr_ub A numeric vector of shape (1, k) that specifies the upper bound constraint
#'                 for each of the k regressors.
#'@param constr_triangle A boolean scalor.
#'                       If TRUE, the optima of the RSM will also be constrained to lie inside
#'                       a triangle defined by the original (0,0) and two other vertices.
#'                       Note only the "quad" model and the "tps" model have this option.
#'                       Dafault is FALSE.
#'@param constr_vertex_1 A numeric vector of shape (1, 2) that specifies one of the other two 
#'                       vertices if constr_triangle is TRUE. 
#'                       (NOTE: vertices numbered clockwise)
#'@param constr_vertex_2 A numeric vector of shape (1, 2) that specifies one of the other two 
#'                       vertices if constr_triangle is TRUE. 
#'                       (NOTE: vertices numbered clockwise)
#'@param verbose A boolean scalor.
#'               If TRUE, function status will be printed.
#'               Default is FALSE.
#'
#' @references{
#'  \insertAllCited{}
#' }
#' 
#' @examples
#' \dontrun{
#' # Example 1: randomly generated 2-variable response surface data
#' X <- cbind(runif(100, -2, 2), runif(100, -2, 2))
#' y <- as.matrix(72 - 11.78 * X[, 1] + 0.74 * X[, 2] - 7.25 * X[, 1]^2 - 7.55 * X[, 2]^2 -
#'   4.85 * X[, 1] * X[, 2] + rnorm(100, 0, 8))
#' # Find a 95 percent confidence region for the maximum of a quadratic polynomial
#' # fitted) to these data
#' out <- OptRegion(
#'   X = X, y = y, B = 200, LB = c(-2, -2), UB = c(2, 2), RSM = "quad"
#' )
#' plot(out, xlab = "X1", ylab = "X2")
#'
#' # Example 2: a mixture-amount experiment in two components (Drug dataset) with
#' # non-normal data. Note triangular experimental region. Resulting 95%
#' # confidence region is pushed against the constraint and results in a
#' # "thin line"
#' out <- OptRegion(
#'   X = Drug[, 1:2], y = Drug[, 3], B = 500, LB = c(0, 0), UB = c(0.08, 11), RSM = "quad",
#'   triangularRegion = TRUE, vertex1 = c(0.02, 11), vertex2 = c(0.08, 1.8)
#' )
#' plot(out, xlab = "Component 1 (mg.)", ylab = "Component 2 (mg.)")
#'
#' # Example 3: randomly generated 2-variable response surface data
#' X <- cbind(runif(100, -2, 2), runif(100, -2, 2))
#' y <- as.matrix(72 - 11.78 * X[, 1] + 0.74 * X[, 2] - 7.25 * X[, 1]^2 -
#'   7.55 * X[, 2]^2 - 4.85 * X[, 1] * X[, 2] + rnorm(100, 0, 8))
#' # Find a 95 percent confidence region for the maximum of a Thin Plate Spline
#' # model fitted to these data
#' out <- OptRegion(X = X, y = y, B = 200, LB = c(-2, -2), UB = c(2, 2), RSM = "tps")
#' plot(out, xlab = "X1", ylab = "X2")
#'
#' # Example 4: a mixture-amount experiment in two components (Drug dataset) with
#' # non-normal data. Note triangular experimental region. Resulting 95p confidence
#' # region of the maxima of a TPS model has area > 0. Contrast with region for
#' # quadratic polynomial model. Note: 500 bootstrap iterations may take a few minutes.
#' out <- OptRegion(
#'   X = Drug[, 1:2], y = Drug[, 3], B = 500, lambda = 0.05,
#'   LB = c(0, 0), UB = c(0.08, 11), RSM = "tps",
#'   triangularRegion = TRUE, vertex1 = c(0.02, 11), vertex2 = c(0.08, 1.8)
#' )
#' plot(out, xlab = "Component 1 (mg.)", ylab = "Component 2 (mg.)")
#'
#' # Example 5: run GloptiPolyRegion on a quadratic, 3 vars example
#' out <- OptRegion(
#'   X = quad_3D[, 1:3], y = quad_3D[, 4], B = 500, alpha = 0.1,
#'   LB = c(-2, -2, -2), UB = c(2, 2, 2), RSM = "poly", degree = 2,
#'   maximization = TRUE, verbose = TRUE
#' )
#' plot(out, c("x1", "x2", "x3"))
#'
#' # Example 6: run GloptiPolyRegion on a cubic, 5 vars example
#' out <- OptRegion(
#'   X = cubic_5D$design_matrix, y = cubic_5D$response, B = 200, alpha = 0.05,
#'   LB = rep(0, 5), UB = rep(5, 5), RSM = "poly", degree = 3,
#'   maximization = TRUE, verbose = TRUE
#' )
#' plot(out, c("x1", "x2", "x3", "x4", "x5"))
#' }
#' @export
OptRegion <- function(X, y, alpha = 0.05, num_boot = 200, 
                      RSM, lambda = 0.04, degree = 2,
                      maximization = TRUE, constr_lb, constr_ub,
                      constr_triangle = FALSE, constr_vertex_1 = NULL, constr_vertex_2 = NULL,
                      verbose = FALSE) {
  constr_vertex_1 <- t(constr_vertex_1)
  constr_vertex_2 <- t(constr_vertex_2)
  if (RSM == "quad") {
    res <- OptRegionQuad(
      X = X, y = y, alpha = alpha, nosim = num_boot,
      maximization = maximization, LB = constr_lb, UB = constr_ub,
      triangularRegion = constr_triangle, vertex1 = constr_vertex_1, vertex2 = constr_vertex_2
    )
  } else if (RSM == "tps") {
    res <- OptRegionTps(
      X = X, y = y, alpha = alpha, nosim = num_boot,
      lambda = lambda,
      maximization = maximization, LB = constr_lb, UB = constr_ub,
      triangularRegion = constr_triangle, vertex1 = constr_vertex_1, vertex2 = constr_vertex_2
    )
  } else if (RSM == "poly") {
    res <- GloptiPolyRegion(
      X = X, y = y, alpha = alpha, B = num_boot, 
      degree = degree, 
      maximization = maximization, lb = constr_lb, ub = constr_ub, 
      verbose = verbose
    )
  } else {
    stop("Function does not support this RSM type.")
  }
  # return
  res
}
