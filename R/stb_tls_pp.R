#' Parametric approximation of Bayesian power prior for multiple binomial historical control data
#'
#'
stb_tl_bayes_bin_pp_control <- function(N,
                                        y,
                                        prior1,
                                        prior2,
                                        power.param,
                                        stan_mdl = "b_pp",
                                        ...){
  ## call stan
  lst_data <- list(N0 = as.array(N),
                   y0 = as.array(y),
                   k = length(N),
                   c = prior1,
                   d = prior2,
                   a_0 = power.param)

  rst <- stb_stan(lst_data, stan_mdl = stan_mdl, ...)
  rst <- rst$draws("theta")
  rst <- posterior::as_draws_rvars(rst)
  rst <- posterior::draws_of(rst$theta)
  PP.c <- RBesT::automixfit(as.vector(rst), type='beta')
  PP.c
}


#' Parametric approximation of Bayesian power prior for multiple normal historical control data
#'
#'
stb_tl_bayes_norm_pp_control <- function(y.m,
                                         y.se,
                                         prior1,
                                         prior2,
                                         power.param,
                                         stan_mdl = "n_pp",
                                         ...){

  lst_data <- list(y0 = as.array(y.m),
                   sigma = as.array(y.se),
                   k = length(y.m),
                   m = prior1,
                   s = prior2,
                   a_0 = power.param)

  rst <- stb_stan(lst_data, stan_mdl = stan_mdl, ...)
  rst <- rst$draws("mu")
  rst <- posterior::as_draws_rvars(rst)
  rst <- posterior::draws_of(rst$mu)
  PP.c <- RBesT::automixfit(as.vector(rst), type='norm')
  
  PP.c

}


#' Parametric approximation of Bayesian power prior
#'
#' Parametric approximation of Bayesian power prior for multiple historical
#' difference between treatment and control
#'
#'
stb_tl_bayes_pp_dif <- function(y.dif,
                                y.dif.se,
                                prior1,
                                prior2,
                                power.param,
                                stan_mdl = "n_pp",
                                ...) {
  ## call stan
  lst_data <- list(y0 = as.array(y.dif),
                   sigma = as.array(y.dif.se),
                   k = length(y.dif),
                   m = prior1,
                   s = prior2,
                   a_0 = power.param)

  rst <- stb_stan(lst_data, stan_mdl = stan_mdl, ...)
  rst <- rst$draws("mu")
  rst <- posterior::as_draws_rvars(rst)
  rst <- posterior::draws_of(rst$mu)
  PP <- RBesT::automixfit(as.vector(rst), type='norm')

  PP
}

stb_tl_bayes_pp_or <- function(logOR,
                               logOR.se,
                               prior1,
                               prior2,
                               power.param,
                               stan_mdl = 'n_pp',
                               ...) {
    ## call stan
    lst_data <- list(y0 = as.array(logOR),
                     sigma = as.array(logOR.se),
                     k = length(logOR),
                     m = prior1,
                     s = prior2,
                     a_0 = power.param)

    rst <- stb_stan(lst_data, stan_mdl = stan_mdl, ...)
    rst <- rst$draws("mu")
    rst <- posterior::as_draws_rvars(rst)
    rst <- posterior::draws_of(rst$mu)
    
    PP <-RBesT::automixfit(as.vector(rst), type='norm')
    PP
}


#' Parametric approximation of Bayesian power prior
#'
#' @export
#'
stb_tl_bayes_para_pp <- function(dtype = c('binary', 'normal'),
                                borrow = c('control', 'diff', 'logOR'),
                                ...){

    if(borrow == 'control'){
        if (dtype=='binary'){
            PP = stb_tl_bayes_bin_pp_control(...)
        } else if(dtype == 'normal') {
            PP = stb_tl_bayes_norm_pp_control(...)
        }
    } else if(borrow == 'diff'){
        PP = stb_tl_bayes_pp_dif(...)
    } else if(borrow == 'logOR'){
        PP = stb_tl_bayes_pp_or(...)
    }

    print(plot(PP))
    return(PP)
}
