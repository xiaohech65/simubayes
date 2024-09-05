## -----------------------------------------------------------------------------
##
##                 STAN
##
## -----------------------------------------------------------------------------

#' Call STAN models
#'
#' Call STAN models. Called by \code{psrwe_powerp}.
#'
#' @param lst_data List of study data to be passed to STAN
#' @param stan_mdl STAN model name
#' @param chains STAN parameter. Number of Markov chainsm
#' @param iter STAN parameter. Number of iterations
#' @param warmup STAN parameter. Number of burnin.
#' @param control STAN parameter. See \code{rstan::stan} for details.
#' @param ... other options to call STAN sampling such as \code{thin},
#'     \code{algorithm}. See \code{rstan::sampling} for details.#'
#'
#' @return Result from STAN sampling
#'
#' @export
#'
stb_stan <- function(lst_data,
                     stan_mdl = c("nb_mix", "nb_pp", "b_pp", "n_pp"),
                     chains = 4, iter = 2000, warmup = 1000, cores = 4,
                     ...) {
    
    stan_mdl <- match.arg(stan_mdl)
    mod <- instantiate::stan_package_model(
      name = stan_mdl,
      package = "simubayes2"
    )
    stan_rst <- mod$sample(data = lst_data,
                           chains = chains,
                           parallel_chains = cores,
                           iter_warmup = warmup,
                           iter_sampling = iter - warmup,
                           adapt_delta = 0.95,
                           ...)

    stan_rst
}


## -----------------------------------------------------------------------------
##
##                 SHINY
##
## -----------------------------------------------------------------------------

#' Run Web-Based application
#'
#' Call Shiny to run \code{statidea} as a web-based application.
#'
#' @details
#'
#' A web browser will be brought up for users to access the GUI
#'
#'
#' @export
#'
stb_shiny <- function(appname = "mstage", pkgname = "simubayes2") {
    req_pkgs        <- c("shiny", "shinythemes", "DT",
                         "shinybusy")

    chk_uninstalled <- sapply(req_pkgs,
                              function(x) {
                                  !requireNamespace(x,
                                                    quietly = TRUE)
                              })

    chk_inx         <- which(chk_uninstalled)

    if (0 < length(chk_inx)) {

        msg <- paste("For the Shiny app to work, please install ",
                     ifelse(1 < length(chk_inx), "packages ", "package "),
                     paste(req_pkgs[chk_inx], collapse = ", "),
                     " by \n install.packages(",
                     paste(paste("'",
                                 req_pkgs[chk_inx],
                                 "'",
                                 sep = ""), collapse = ", "),
                     ") \n  ",
                     sep = "")

        stop(msg, call. = FALSE)
    }

    app_dir <- system.file(appname, package = pkgname)
    if (app_dir == "") {
        stop("Could not find Shiny directory. Try re-installing the package.",
             call. = FALSE)
    }

    shiny::runApp(app_dir, display.mode = "normal")
}
