## ROI plugin: clarabel
## based on clarabel interface
make_clarabel_signatures <- function()
    ROI_plugin_make_signature( objective = c("L", "Q"),
                               constraints = c("X", "L", "C"),
                               types = c("C"),
                               bounds = c("X", "C", "V", "CV"),
                               cones = c("X", "zero", "nonneg", "soc", "psd", "expp", "powp"),
                               maximum = c(TRUE, FALSE) )


## SOLVER CONTROLS
.add_controls <- function(solver) {
    ## Main algorithm settings
    ROI_plugin_register_solver_control( solver, "max_iter", "max_iter" )
    ROI_plugin_register_solver_control( solver, "time_limit", "max_time" )
    ROI_plugin_register_solver_control( solver, "verbose", "verbose" )
    ROI_plugin_register_solver_control( solver, "max_step_fraction", "X")
    ## Full accuracy settings
    ROI_plugin_register_solver_control( solver, "tol_gap_abs", "abs_tol" )
    ROI_plugin_register_solver_control( solver, "tol_gap_rel", "rel_tol" )
    ROI_plugin_register_solver_control( solver, "tol_feas", "tol" )
    ROI_plugin_register_solver_control( solver, "tol_infeas_abs", "X" )
    ROI_plugin_register_solver_control( solver, "tol_infeas_rel", "X" )
    ROI_plugin_register_solver_control( solver, "tol_ktratio", "X" )
    ## Reduced accuracy settings
    ROI_plugin_register_solver_control( solver, "reduced_tol_gap_abs", "X" )
    ROI_plugin_register_solver_control( solver, "reduced_tol_gap_rel", "X" )
    ROI_plugin_register_solver_control( solver, "reduced_tol_feas", "X" )
    ROI_plugin_register_solver_control( solver, "reduced_tol_infeas_abs", "X" )
    ROI_plugin_register_solver_control( solver, "reduced_tol_infeas_rel", "X" )
    ROI_plugin_register_solver_control( solver, "reduced_tol_ktratio", "X" )
    ## data equilibration settings
    ROI_plugin_register_solver_control( solver, "equilibrate_enable", "X" )
    ROI_plugin_register_solver_control( solver, "equilibrate_max_iter", "X" )
    ROI_plugin_register_solver_control( solver, "equilibrate_min_scaling", "X" )
    ROI_plugin_register_solver_control( solver, "equilibrate_max_scaling", "X" )
    ## Step size settings
    ROI_plugin_register_solver_control( solver, "linesearch_backtrack_step", "X" )
    ROI_plugin_register_solver_control( solver, "min_switch_step_length", "X" )
    ROI_plugin_register_solver_control( solver, "min_terminate_step_length", "X" )
    ## Linear solver settings
    ROI_plugin_register_solver_control( solver, "direct_kkt_solver", "X" )
    ROI_plugin_register_solver_control( solver, "direct_solve_method", "method" )
    ## static regularization parameters
    ROI_plugin_register_solver_control( solver, "static_regularization_enable", "X" )
    ROI_plugin_register_solver_control( solver, "static_regularization_constant", "X" )
    ROI_plugin_register_solver_control( solver, "static_regularization_proportional", "X" )
    ## dynamic regularization parameters
    ROI_plugin_register_solver_control( solver, "dynamic_regularization_enable", "X" )
    ROI_plugin_register_solver_control( solver, "dynamic_regularization_eps", "X" )
    ROI_plugin_register_solver_control( solver, "dynamic_regularization_delta", "X" )
    ROI_plugin_register_solver_control( solver, "iterative_refinement_enable", "X" )
    ROI_plugin_register_solver_control( solver, "iterative_refinement_reltol", "X" )
    ROI_plugin_register_solver_control( solver, "iterative_refinement_abstol", "X" )
    ROI_plugin_register_solver_control( solver, "iterative_refinement_max_iter", "X" )
    ROI_plugin_register_solver_control( solver, "iterative_refinement_stop_ratio", "X" )
    ROI_plugin_register_solver_control( solver, "presolve_enable", "presolve" )
    
    invisible( TRUE )
}

.onLoad <- function( libname, pkgname ) {
    ## Solver plugin name (based on package name)
    if( ! pkgname %in% ROI_registered_solvers() ){
        ## Register solver methods here.
        ## One can assign several signatures a single solver method
        solver <- ROI_plugin_get_solver_name( pkgname )
        ROI_plugin_register_solver_method( 
            signatures = make_clarabel_signatures(),
            solver = solver,
            method = getFunction( "solve_OP", where = getNamespace(pkgname)) )
        ## Finally, for status code canonicalization add status codes to data base
        .add_status_codes()
        .add_controls( solver )
    }
}

