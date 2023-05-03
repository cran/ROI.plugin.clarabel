## ROI plugin: clarabel
## based on clarabel interface

cone_dims <- function(x, ctype) {
    wcol <- which(colnames(x) == ctype)
    if ( length(wcol) == 0 ) return(NULL)
    as.integer(table(x$v[x$j == wcol]))
}

cone_counts <- function(x, ctype) {
    wcol <- which(colnames(x) == ctype)
    if ( length(wcol) == 0 ) return(0)
    sum(x$v[x$j == wcol])
}


to_dense_vector <- function(x, len) {
    y <- rep.int(0L, len)
    if ( is.null(x$ind) ) return(y)
    y[x$ind] <- x$val
    y
}

calc_zero_dims <- function(x) {
    y <- sum(x$cone == clarabel_cones["zero"])
    if ( y == 0L ) return(NULL)
    y
}

calc_lin_dims <- function(x) {
    y <- sum(x$cone == clarabel_cones["nonneg"])
    if ( y == 0L ) return(NULL)
    y
}

calc_expp_dims <- function(x) {
    y <- x$id[x$cone == clarabel_cones['expp']]
    if ( !length(y) ) return(NULL)
    length(unique(y))
}

calc_soc_dims <- function(x) {
    y <- x$id[x$cone == clarabel_cones['soc']]
    if ( !length(y) ) return(NULL)
    as.integer(table(y))
}

calc_pow_dims <- function(x) {
    powp <- powd <- NULL
    ids <- unique(x$id[x$cone == clarabel_cones['powp']])
    if ( length(ids) )
        powp <- sapply(as.character(ids), function(id) x$params[[id]]['a'], USE.NAMES = FALSE)

    # ids <- unique(x$id[x$cone == clarabel_cones['powd']])
    # if ( length(ids) )
    #     powd <- sapply(as.character(ids), function(id) -x$params[[id]]['a'], USE.NAMES = FALSE)

    unname(c(powp, powd))
}

calc_dims <- function(cones) {
    dims <- list()
    dims$z <- calc_zero_dims(cones)
    dims$l <- calc_lin_dims(cones)
    dims$q <- calc_soc_dims(cones)
    dims$ep <- calc_expp_dims(cones)
    #dims$s <- calc_psd_dims(cones)
    dims$p <- calc_pow_dims(cones)

    o_class <- unlist(lapply(dims, class), FALSE, FALSE)
    o_lens <- lapply(dims, length)
    cl <- rep(o_class, o_lens)
    dims <- as.list(unlist(dims))
    for (i in seq(dims)) {
        class(dims[[i]]) <- cl[[i]]
    }
    dims
}

clarabel_cones <-  c("zero" = 1L, "nonneg" = 2L, "soc" = 3L, "expp" = 5L, "powp" = 7L)


convert_mat <- function(x) {
    ind <- order(x$j, x$i)
    list(matbeg = c(0L, cumsum(tabulate(x$j[ind], x$ncol))),
         matind = x$i[ind] - 1L,
         values = x$v[ind])
}


solve_OP <- function(x, control = list()) {

    constr <- as.C_constraint(constraints(x))

    ## check if "clarabel" supports the provided cone types
    stopifnot(all(constr$cones$cone %in% clarabel_cones))

    q <- as.vector(terms(objective(x))[["L"]])
    P <- terms(objective(x))[["Q"]]

    if ( !is.null(P) ) P <- as.matrix(P)
    if ( maximum(x) ) {
        q <- -q
        if ( !is.null(P) ) P <- -P
    }

    AL <- AU <- NULL
    AL.rhs <- AU.rhs <- double()

    ## lower bounds
    lower_bounds <- to_dense_vector(bounds(x)$lower, length(objective(x)))
    not_is_clarabel_default <- !is.infinite(lower_bounds)
    if ( any(not_is_clarabel_default) ) {
        li <- which(not_is_clarabel_default)
        AL <- simple_triplet_matrix(i = seq_along(li), j = li, 
                                    v = rep.int(-1, length(li)), 
                                    nrow = length(li), ncol = length(q))
        AL.rhs <- -lower_bounds[not_is_clarabel_default]
    }

    ## upper bounds
    ui <- bounds(x)$upper$ind
    ub <- bounds(x)$upper$val
    if ( length(ui) ) {
        AU <- simple_triplet_matrix(i = seq_along(ui), j = ui,
                                    v = rep.int(1, length(ui)), 
                                    nrow = length(ui), ncol = length(q))
        AU.rhs <- ub
    }

    A <- rbind(constr$L, AL, AU)
    A.rhs <- c(constr$rhs, AL.rhs, AU.rhs)
    cones <- c(constr$cones, K_lin(length(AL.rhs)), K_lin(length(AU.rhs)))

    if ( nrow(constr) > 0 ) {
        i <- with(cones, order(cone, id))
        ordered_cones <- list(cone = cones$cone[i], id = cones$id[i])
        A <- A[i,]
        A.rhs <- A.rhs[i]
        dims <- calc_dims(cones)
    } else {
        dims <- list()
    }    

    if ( is.null(control$verbose) ) control$verbose <- FALSE

    solver_call <- list(clarabel, A = A, b = A.rhs, q = q, P = P, cone = dims, control = control, strict_cone_order=FALSE)
    mode(solver_call) <- "call"
    if ( isTRUE(control$dry_run) )
        return(solver_call) 

    out <- eval(solver_call)
    out$len_objective <- length(objective(x))
    out$len_dual_objective <- nrow(constraints(x))

    objval <- tryCatch(objective(x)(out$x), error = function(e) as.numeric(NA))
    ROI_plugin_canonicalize_solution(solution = out$x,
                                     optimum  = objval,
                                     status   = out[["status"]],
                                     solver   = "clarabel",
                                     message  = out)
}
