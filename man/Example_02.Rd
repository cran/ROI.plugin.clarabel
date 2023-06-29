\name{Example-2}
\title{SOCP 2}
\alias{ROI.plugin.clarabel_Example_2}
\description{
The following example is also known as \code{Problem 10} from the 
Hock-Schittkowski-Collection \code{Hock and Schittkowski (1981)}.

\deqn{minimize \ \ x - y \\}
\deqn{subject \ to \ \ -3 x^2 + 2 x y  + 1 \geq 0}
}

\references{
W. Hock, K. Schittkowski (1981): Test Examples for Nonlinear Programming Codes, Lecture Notes in Economics and Mathematical Systems, Vol. 187, Springer
}

\examples{

Sys.setenv("ROI_LOAD_PLUGINS" = FALSE)
library(ROI)
library(ROI.plugin.clarabel)

obj <- L_objective(c(1, -1))
L <- chol(rbind(c(3, -1), c(-1, 1)))
con <- C_constraint(L = rbind(0, -L), cones = K_soc(3), rhs = c(1, 0, 0))
op <- OP(objective = obj, constraints = con, 
         bounds = V_bound(li = 1:2, lb = rep(-Inf, 2)))
x <- ROI_solve(op, solver="clarabel")
x
## Optimal solution found.
## The objective value is: -1.000000e+00
solution(x)
## [1] -4.622464e-16  1.000000e+00


}
