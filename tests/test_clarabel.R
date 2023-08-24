if (isTRUE(Sys.getenv("_R_ROI_NO_CHECK_SOLVERS_") == "")) {
    Sys.setenv("ROI_LOAD_PLUGINS" = FALSE)
}
library(ROI)
library(ROI.plugin.clarabel)

check <- function(domain, condition, level=1, message="", call=sys.call(-1L)) {
    if ( isTRUE(condition) ) return(invisible(NULL))
    msg <- sprintf("in %s", domain)
    if ( all(nchar(message) > 0) ) msg <- sprintf("%s\n\t%s", msg, message)
    stop(msg)
    return(invisible(NULL))
}


## SOCP - Example - 1
## min:  1 x1 + 1 x2 + 1 x3
## s.t.     x1 == sqrt(2)
##          x1 >= ||(x2, x3)||
##
## c(sqrt(2), -1, -1)
test_cp_01 <- function(solver) {
    obj <- c(1, 1, 1)
    A <- rbind(c(1, 0, 0))
    b <- c(sqrt(2))
    G <- diag(x=-1, 3)
    h <- rep(0, 3)

    bound <- V_bound(li = 1:3, lb = rep(-Inf, 3))

    lc <- C_constraint(L = rbind(A, G), 
                       cones = c(K_zero(1), K_soc(3)), 
                       rhs = c(b, h))
    x <- OP(objective = obj, constraints = lc, types = rep("C", 3),
            bounds =  bound, maximum = FALSE)

    opt <- ROI_solve(x, solver = solver)
    
    check("CP-01@01", equal(sum(abs(opt$solution - c(sqrt(2), -1, -1))), 0))
    check("CP-01@02", equal(opt$objval, (sqrt(2) - 2)))
}

## SOCP - Example - 2
## min:  0 x1 - 2 x2 - 2 x3 + 0 x4 - 2 x5 - 2 x6
## s.t.     x1 == sqrt(2)
##          x4 == sqrt(2)
##          x1 >= ||(x2, x3)||
##          x4 >= ||(x5, x6)||
##
## c(sqrt(2), 1, 1, sqrt(2), 1, 1)
test_cp_02 <- function(solver) {
    obj <- c(0, -2, -2, 0, -2, -2)
    A <- rbind(c(1, 0, 0, 0, 0, 0),
               c(0, 0, 0, 1, 0, 0))
    b <- c(sqrt(2), sqrt(2))
    G <- diag(x=-1, 6)
    h <- rep(0, 6)

    lc <- C_constraint(L = rbind(A, G), 
                       cones = c(K_zero(2), K_soc(c(3, 3))), 
                       rhs = c(b, h))
    x <- OP(objective = obj, constraints = lc)

    opt <- ROI_solve(x, solver=solver)
    check("CP-02@01", equal(sum(abs(opt$solution - c(sqrt(2), 1, 1, sqrt(2), 1, 1))), 0))
}

## EXPP - Example - 1
## min:  x + y + z
## s.t.
## y e^(x/y) <= z
## y > 0
## x := 1
## y := 2
## c(1, 2, 2*exp(1/2))
test_cp_03 <- function(solver) {
    obj <- c(1, 1, 1)
    A <- rbind(c(1, 0, 0),
               c(0, 1, 0))
    b <- c(1, 2)
    G <- -diag(3)
    h <- rep(0, 3)

    lc <- C_constraint(L = rbind(A, G), 
                       cones = c(K_zero(2), K_expp(1)), 
                       rhs = c(b, h))
    x <- OP(objective = obj, constraints = lc)

    opt <- ROI_solve(x, solver = solver)
    check("CP-03@01", max(abs(opt$solution - c(1, 2, 2*exp(1/2)))) < 1e-3)
}

## EXPP - Example - 2
## max:  x + y + z
## s.t.
## y e^(x/y) <= z
## y > 0
## y == 2
## z == 2 * exp(1/2)
## c(1, 2, 2*exp(1/2))
test_cp_04 <- function(solver) {
    obj <- c(1, 1, 1)
    A <- rbind(c(0, 1, 0),
               c(0, 0, 1))
    b <- c(2, 2*exp(1/2))
    G <- diag(x=-1, 3)
    h <- rep(0, 3)

    lc <- C_constraint(L = rbind(A, G), 
                       cones = c(K_zero(2), K_expp(1)), 
                       rhs = c(b, h))
    x <- OP(objective = obj, constraints = lc, maximum = TRUE)

    opt <- ROI_solve(x, solver=solver)
    check("CP-04@01", equal(opt$solution , c(1, 2, 2*exp(1/2))))
}

## EXPP - Example - 3
## max:  x + y + z
## s.t.
## y e^(x/y) <= z
## y > 0
## y == 1
## z == exp(1)
## c(1, 1, exp(1))
test_cp_05 <- function(solver) {
    obj <- c(1, 1, 1)
    A <- rbind(c(0, 1, 0),
               c(0, 0, 1))
    b <- c(1, exp(1))
    G <- diag(x=-1, 3)
    h <- rep(0, 3)
    ## cones <- list("free"=c(1, 2), "expp"=list(3:5))
    ## bound <- as.C_bound(cones)

    lc <- C_constraint(L = rbind(A, G), 
                       cones = c(K_zero(2), K_expp(1)), 
                       rhs = c(b, h))
    x <- OP(objective = obj, constraints = lc, 
            types = rep("C", 3), maximum = TRUE)

    opt <- ROI_solve(x, solver = solver)
    check("CP-05@01", equal(opt$solution , c(1, 1, exp(1))))
}

## POWP - Example - 1
## max:  x + y + z
## s.t.
##      x^a * y ^ (1-a) >= |z|
##      x == 4
##      y == 4
##      a == 1/2
##
## c(4, 4, 4)
test_cp_07 <- function(solver) {
    obj <- c(1, 1, 1)
    A <- rbind(c(1, 0, 0),
               c(0, 1, 0))
    b <- c(4, 4)
    G <- diag(x=-1, 3)
    h <- rep(0, 3)

    cc <- C_constraint(L = rbind(A, G), 
                       cones = c(K_zero(2), K_powp(0.5)), 
                       rhs = c(b, h))
    x <- OP(objective = obj, constraints = cc, 
            types = rep("C", 3), maximum = TRUE)

    opt <- ROI_solve(x, solver=solver)
    check("CP-07@01", equal(opt$solution, c(4, 4, 4)))
    check("CP-07@02", equal(opt$objval, 12 ))
}


## QP - Example - 1
##
## from the quadprog package
## (c) S original by Berwin A. Turlach R port by Andreas Weingessel
## GPL-3
##
## min: -(0 5 0) %*% x + 1/2 x^T x
## under the constraints:      A^T x >= b
## with b = (-8,2,0)^T
## and      (-4  2  0)
##      A = (-3  1 -2)
##          ( 0  0  1)
## we can use solve.QP as follows:
##
## library(quadprog)
## D <- diag(1, 3)
## d <- c(0, 5, 0)
## A <- cbind(c(-4, -3, 0),
##            c( 2,  1, 0),
##            c( 0, -2, 1))
## b <- c(-8, 2, 0)
##
## sol <- solve.QP(D, d, A, bvec=b)
## deparse(sol$solution)
## deparse(sol$value)
test_qp_01 <- function(solver) {
    A <- cbind(c(-4, -3, 0),
               c( 2,  1, 0),
               c( 0, -2, 1))
    x <- OP(Q_objective(diag(3), L =  c(0, -5, 0)),
            L_constraint(L = t(A),
                         dir = rep(">=", 3),
                         rhs = c(-8, 2, 0)))

    opt <- ROI_solve(x, solver = solver)
    solution <- c(0.476190476190476, 1.04761904761905, 2.0952380952381)
    check("QP-01@01", equal(opt$solution, solution) )
    check("QP-01@02", equal(opt$objval, -2.38095238095238) )
}

## This Test detects non-conform objective functions.
## minimize 0.5 x^2 - 2 x + y
## s.t. x <= 3
## Type 1:   0.5 x'Qx + c'Lx => c(2, 0)  objval=-2
## Type 2:       x'Qx + c'Lx => c(3, 0)  objval=-3.75
test_qp_02 <- function(solver) {
    zero <- .Machine$double.eps * 100
    qo <- Q_objective(Q=rbind(c(1, 0), c(0, zero)), L=c(-2, 1))
    lc1 <- L_constraint(L=matrix(c(1, 0), nrow=1), dir="<=", rhs=3)
    lc2 <- L_constraint(L=matrix(c(1, 0), nrow=1), dir=">=", rhs=0)
    x <- OP(qo, c(lc1, lc2))

    opt <- ROI_solve(x, solver=solver)
    solution <- c(2, 0)
    check("QP-02@01", equal(opt$solution, solution) )
    check("QP-02@02", equal(opt$objval, -2) )
}

## as qp_01 but maximize
test_qp_03 <- function(solver) {
    A <- cbind(c(-4, -3, 0),
               c( 2,  1, 0),
               c( 0, -2, 1))
    x <- OP(Q_objective(-diag(3), L = -c(0, -5, 0)),
            L_constraint(L = t(A),
                         dir = rep(">=", 3),
                         rhs = c(-8, 2, 0)),
            maximum = TRUE)

    opt <- ROI_solve(x, solver=solver)
    solution <- c(0.476190476190476, 1.04761904761905, 2.0952380952381)
    check("QP-01@01", equal(opt$solution, solution) )
    check("QP-01@02", equal(opt$objval, 2.38095238095238) )
}

## SDP - Example - 1
## for the example definition see ROI.plugin.scs inst/doc
## or http://cvxopt.org/userguide/coneprog.html
test_cp_09 <- function(solver) {
    ## this function or something similar should go into ROI
    obj <- c(1, -1, 1)
    A1 <- matrix(c(-7, -11, -11,  3), 2)
    A2 <- matrix(c( 7, -18, -18,  8), 2)
    A3 <- matrix(c(-2,  -8,  -8,  1), 2)
    a  <- matrix(c(33,  -9,  -9, 26), 2)
    B1 <- matrix(c(-21, -11,  0, -11,  10,   8,  0,    8, 5), 3)
    B2 <- matrix(c(  0,  10,  16, 10, -10, -10,  16, -10, 3), 3)
    B3 <- matrix(c( -5,   2, -17,  2,  -6,   8, -17,   8, 6), 3)
    b  <- matrix(c( 14,   9,  40,  9,  91,  10,  40,  10,15), 3)

    ## PSD matrices have to be vectorizedb <- ordered_cones$cone == scs_cones["psd"]
    G1 <- vech(A1, A2, A3, lower=FALSE, scale=TRUE)
    h1 <- vech(a, lower=FALSE, scale=TRUE)
    G2 <- vech(B1, B2, B3, lower=FALSE, scale=TRUE)
    h2 <- vech(b, lower=FALSE, scale=TRUE)
    h <- c(h1, h2)
    bounds <- V_bound(li=1:3, lb=rep(-Inf, 3)) 

    x <- OP(objective = obj,
            constraints = C_constraint(L = rbind(G1, G2), 
                                       cones = K_psd(c(3, 6)), 
                                       rhs = h),
            types = rep("C", length(obj)),
            bounds =  bounds,
            maximum = FALSE)

    opt <- ROI_solve(x, solver = solver)
    
    ## NOTE: The solutions I compare with are from cvxopt where I used the default settings,
    ##       therefore it is possible that scs just provides a solution with a smaler eps
    sol <- c(-0.36775, 1.89833, -0.88746)
    check("CP-09@01", equal(solution(opt, "objval"), drop(obj %*% sol)))
    
    ## solution from cvxopt
    ## [-3.68e-01 1.90e+00 -8.88e-01]
    ## or c(-0.367666090041563, 1.89832827158511, -0.887550426343585)
    check("CP-09@02", isTRUE(sum(abs(solution(opt) - sol)) < 1e-3))

    ## [ 3.96e-03 -4.34e-03]
    ## [-4.34e-03  4.75e-03]
    ## c(0.00396107103000518, -0.00433836779348354, -0.00433836779348354,  0.00475162592559036) 
    sol_psd_1 <- c( 0.00396107103000518, -0.00433836779348354, 
                   -0.00433836779348354,  0.00475162592559036)

    opt_sol_psd_1 <- as.numeric(as.matrix(solution(opt, "psd")[[1]]))
    check("CP-09@03", isTRUE(sum(abs(opt_sol_psd_1 - sol_psd_1)) < 1e-5))
    
    ## [ 5.58e-02 -2.41e-03  2.42e-02]
    ## [-2.41e-03  1.04e-04 -1.05e-03]
    ## [ 2.42e-02 -1.05e-03  1.05e-02]
    ## c(0.0558011514407859, -0.00240909203896524, 0.0242146296992217,  -0.00240909203896524, 
    ##   0.000104021271556218, -0.00104543254168053,  0.0242146296992217, -0.00104543254168053, 
    ##   0.0105078600239678) 
    sol_psd_2 <- c( 0.0558011514407859, -0.00240909203896524,   0.0242146296992217,  
                   -0.00240909203896524, 0.000104021271556218, -0.00104543254168053,  
                    0.0242146296992217, -0.00104543254168053,   0.0105078600239678)
    opt_sol_psd_2 <- as.numeric(as.matrix(solution(opt, "psd")[[2]]))
    check("CP-09@04", isTRUE(max(abs(opt_sol_psd_2 - sol_psd_2)) < 1e-4))
}



if ( !any("clarabel" %in% names(ROI_registered_solvers())) ) {
    ## This should never happen.
    cat("ROI.plugin.clarabel could not be found among the registered solvers.\n")
} else {
    print("Start Testing!")
    local({test_cp_01("clarabel")})
    local({test_cp_02("clarabel")})
    local({test_cp_03("clarabel")})
    local({test_cp_04("clarabel")})
    local({test_cp_05("clarabel")})
    local({test_cp_07("clarabel")})
    local({test_qp_01("clarabel")})
    local({test_qp_02("clarabel")})
    local({test_qp_03("clarabel")})
    local({test_cp_09("clarabel")})
}
