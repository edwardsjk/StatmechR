##' A function to fit a generic epi model given an epi model and a 
##' curve generating function. Currently uses optim with the defualt
##' parameters
##'
##' @param ecs a list of epidemic curves we are trying to fit.
##' @param epi_mdl_func the epidemic model we are trying to fit. Should take a set of parameters and produce a curve of the given length
##' @param epi_mdl_pars a data frame containing the starting parameters for the epidemic model
##' @param error_func the error function we are trying to minimze
##' @param ... additional parameters to the objective function
##'
##' @return a data frame with the epidemic model parameters and the covergence result.
fit_epi_mdl <- function (ecs, epi_mdl_func, epi_mdl_pars, error_func, ...) {
    
    ##' takes in the parameter set and the epidemic curve
    obj_fxn <- function(par, ec) {
        pred_curve <- epi_mdl_func(par, length(ec))
        err <- error_func(pred_curve, ec, ...)
        return(error)
    }

    converged <- rep(0, legnth(ecs))

    for (i in 1:length(ecs)) {
        ec <- ecs[[i]]

        tmp <- optim(epi_mdl_pars[i,],obj_fxn, ec=ec)

        coverged[i]<- tmp$convergence
        epi_mdl_pars[i,] <- tmp$par

    }

    epi_mdl_pars$converged <- converged

    return(epi_mdl_pars)
}

