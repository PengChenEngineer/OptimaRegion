#' Encode using the orthogonal convention
#' 
#' @param Xi numeric matrix of shape (N, k); Xi specifies the points in
#'           original units
#' @param lb numeric vector of shape (1, k); lb specifies the lower bounds for
#'           the k variables in original units
#' @param ub numeric vector of shape (1, k); ub specifies the upper bounds for 
#'           the k variables in original units
#' @return numeric matrix of shape (N, k); it specifies the points in the
#'         (-1, 1) coded units
encode_orthogonal <- function(Xi, lb, ub){
  Xi <- t(Xi) # put data in columns for vectorized implementation
  M <- (lb + ub) / 2; dim(M) <- c(nrow(Xi), 1); res <- sweep(Xi, 1, M, "-")
  R <- ub - lb; S <- diag(R) / 2; res <- solve(S, res)
  # return
  t(res)
}

#' Decode using the orthogonal convention
#' 
#' @param X numeric matrix of shape (N, k); X specifies the points in
#'          original units
#' @inheritParams encode_orthogonal
#' @return numeric matrix of shape (N, k); res specifies the points in
#'         original units    
decode_orthogonal <- function(X, lb, ub){
  X <- t(X) 
  M <- (lb + ub) / 2; dim(M) <- c(nrow(X), 1);
  R <- ub - lb; S <- diag(R) / 2
  res <- S %*% X; res <- sweep(res, 1, M, "+")
  #return
  t(res)
}