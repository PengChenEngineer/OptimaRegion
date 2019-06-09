# OptimaRegion

Compute the confidence region on the optima of a response surface model

Contributors: Enrique del Castillo, Peng Chen, Adam Meyers, John Hunt, and James Rapkin

## Note
- rcmdcheck::rcmdcheck(build_args = c('--compact-vignettes=gs+qpdf'))
- system("R CMD Rd2pdf . --title=Package OptimaRegion --output=./manual.pdf --force --no-clean --internals")

## To do 
- [ ] vignettes vs. qpdf warning
- [ ] Make arguments' names consistent between old and new functions
- [ ] Change xlab and ylab for old functions?

### 3D View
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
