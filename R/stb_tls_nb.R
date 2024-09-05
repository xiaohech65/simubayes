## -----------------------------------------------------------------------------
##
## This file contains R tool functions for Bayesian analysis.
##
## All functions start with "stb_tl_bayes_".
##
##
##
## -----------------------------------------------------------------------------

#' Bayesian informative prior for negative binomial data
#'
#' Draw posterior samples of log risk ratio for negative binomial study.
#'
#' @param cur_data
#' @param ext_prior
#' @param pri_mu0
#'
#'
#' @param ... Parameters for STAN sampling
#'
#'
#' @export
#'
stb_tb_bayes_nb_borrow <- function(cur_data,
                                   ext_prior,
                                   pri_mu0,
                                   pri_k,
                                   pri_log_rr_weak,
                                   weight,
                                   model = c("nb_mix", "nb_pp"),
                                   seed = NULL,
                                   ...) {

    if (!is.null(seed)) {
        old_seed <- set.seed(seed)
    }

    model <- match.arg(model)

    ## orgnize data
    ## y0  <-

    ## call stan
    lst_data <- list()
    ## lst_data <- list(N0 = n0,
    ##                  N1 = n1,
    ##                  )
    rst <- stb_stan(lst_data, stan_mdl = model, ...)

    ## reset random seed
    if (!is.null(seed)) {
        set.seed(old_seed)
    }

    ## return
    rst <- extract(rst, par = c('theta', 'mu0', "k"))
    rst
}


#' Simulate Negative Binomial with Event time
#'
#' Simulate negative binomial by event time
#'
#' @param n sample size
#' @param mu mean of the negative binomial
#' @param k  dispersion parameter
#' @param n_event number of event time
#'
#' @export
#'
stb_tl_rc_simu_events <- function(n, mu, k = 1, n_event = 5, seed = NULL) {

    if (!is.null(seed))
        old_seed <- set.seed(seed)

    f_single <- function(x) {
        lambda  <- rgamma(1, k, k) * mu
        ts      <- rexp(n_event, lambda)
        cum_ts  <- cumsum(ts)
        t_start <- c(0,      cum_ts)
        t_end   <- c(cum_ts, Inf)

        data.frame(id     = x,
                   inx    = seq_len(n_event + 1),
                   lambda = lambda,
                   start  = t_start,
                   end    = t_end,
                   time   = t_end - t_start)
    }

    rst <- lapply(1:n, function(x) f_single(x))

    if (!is.null(seed))
        set.seed(old_seed)

    rbindlist(rst)
}

#' Censor Recurrent Event Data
#'
#'
#' @param dat_rc recurrent events
#' @param t follow up time
#'
#' @export
#'
stb_tl_rc_events_censor <- function(dat_rc, t = 1) {
    dat_rc %>%
        filter(day_start < t) %>%
        mutate(censor  = if_else(day_end <= t, 0,       1),
               day_end = if_else(day_end <= t, day_end, t),
               time    = day_end - day_start)
}


#' Convert Recurrent Events to NB
#'
#'
#' @param dat_rc recurrent events data frame
#' @param t follow up time
#'
#' @export
#'
stb_tl_rc_events_to_nb <- function(dat_rc, t = 1, var_id = "sid") {

    dat_rc$sid      <- dat_rc[[var_id]]
    dat            <- stb_tl_rc_events_censor(dat_rc, t = t)
    dat$day_censor <- t

    ## minus 1 since there is always a censored row
    dat_count <- dat %>%
        group_by(sid) %>%
        summarize(y = n() - 1) %>%
        mutate(y = if_else(is.na(y), 0, y))

    rst <- dat_count %>%
        left_join(dat_rc %>%
                  select(sid, lambda, day_enroll) %>%
                  distinct(),
                  by = "sid")

    rst[[var_id]] <- rst$sid

    list(data_event = dat,
         data_nb    = rst)
}


#' Generate censoring data
#'
#'
rcurrent_censor <- function(data_full) {
    data_full %>%
        filter(day_start <= day_eos) %>%
        mutate(censor  = if_else(day_end <= day_eos, 0, 1),
               day_end = if_else(day_end <= day_eos, day_end, day_eos),
               time    = day_end - day_start)
}
