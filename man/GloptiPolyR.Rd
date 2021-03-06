% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GloptiPolyR.R
\name{GloptiPolyR}
\alias{GloptiPolyR}
\title{Global optimization of up to cubic polynomial functions (up to 5 variables)}
\usage{
GloptiPolyR(P)
}
\arguments{
\item{P}{A list of list; Each sub-list has 2 components:
1. a multi-dimensional array corresponding to a objective or constraint function
2. an attribute of the objective or constraint function}
}
\value{
Returns the optimal solution and its corresponding objective value
}
\description{
Optimize a quadratic or cubic polynomial functionin 2 ~ 5 variables
with bound constraints
\insertCite{DelCastilloCR}{OptimaRegion}.
}
\details{
GloptipolyR is an R implementation of the “Gloptipoly” algorithm
\insertCite{lasserre2001global}{OptimaRegion}
}
\section{Author(s)}{

Enrique del Castillo \email{exd13@psu.edu},
Peng Chen \email{pfc5098@psu.edu},
Adam Meyers \email{akm5733@psu.edu},
John Hunt \email{J.Hunt@exeter.ac.uk} and
James Rapkin \email{jr297@exeter.ac.uk}.
}

\examples{
# Optimize the following quadratic function in 3 variables
# f(x) = -1.5 x_1 + 2.13 x_2 - 1.81 x_3 + 7.13 x_1 x_2 +
#        3.27 x_1 x_3  + 2.73 x_2 x_3 +
#        4.69 x_1^2 + 6.27 x_2^2 + 5.21 x_3^2.

# The input for GloptiPolyR is a list of 7 sub-lists,
# each of which corresponds to the objective function or a constraint
# function, respectively. See del Castillo et al. (2019) for details.
P <- list()
p_f <- list()
p_g_1 <- list()
p_g_2 <- list()
p_g_3 <- list()
p_g_4 <- list()
p_g_5 <- list()
p_g_6 <- list()

p_f$c <- array(0, dim = c(3, 3, 3))
p_f$c[2, 1, 1] <- -1.5
p_f$c[1, 2, 1] <- 2.13
p_f$c[1, 1, 2] <- -1.81
p_f$c[2, 2, 1] <- 7.13
p_f$c[2, 1, 2] <- 3.27
p_f$c[1, 2, 2] <- 2.73
p_f$c[3, 1, 1] <- 4.69
p_f$c[1, 3, 1] <- 6.27
p_f$c[1, 1, 3] <- 5.21

p_g_1$c <- array(0, dim = c(3, 3, 3))
p_g_1$c[1, 1, 1] <- 2
p_g_1$c[2, 1, 1] <- 1

p_g_2$c <- array(0, dim = c(3, 3, 3))
p_g_2$c[1, 1, 1] <- -2
p_g_2$c[2, 1, 1] <- 1

p_g_3$c <- array(0, dim = c(3, 3, 3))
p_g_3$c[1, 1, 1] <- 2
p_g_3$c[1, 2, 1] <- 1

p_g_4$c <- array(0, dim = c(3, 3, 3))
p_g_4$c[1, 1, 1] <- -2
p_g_4$c[1, 2, 1] <- 1

p_g_5$c <- array(0, dim = c(3, 3, 3))
p_g_5$c[1, 1, 1] <- 2
p_g_5$c[1, 1, 2] <- 1

p_g_6$c <- array(0, dim = c(3, 3, 3))
p_g_6$c[1, 1, 1] <- -2
p_g_6$c[1, 1, 2] <- 1

# Set the attribute for the objective function as either ``min'' or ``max''.
p_f$t <- "min"

# Set the attributes for the constraint functions as either ``>='' or ``<=''.
p_g_1$t <- ">="
p_g_2$t <- "<="
p_g_3$t <- ">="
p_g_4$t <- "<="
p_g_5$t <- ">="
p_g_6$t <- "<="

# Now we put together the input P and use it to call GloptiPolyR
P <- list(p_f, p_g_1, p_g_2, p_g_3, p_g_4, p_g_5, p_g_6)
GloptiPolyR(P)
}
\references{
{
 \insertAllCited{}
}
}
