## -----------------------------------------------------------------------------
##
##  DESCRIPTION:
##      Bayesian Design for Binary Endpoint
##
##  DATE:
##      July, 2024
## -----------------------------------------------------------------------------

#' Describe the design
#'
#'
bayes_mstage_bin_describe <- function(x, ...) {
    cat("Type:\n")
    cat("    Bayesian multi-stage design for binary outcomes with historical borrowing.\n\n")

    cat("Design Parameters:\n")
    cat("    sample_size:     stage 1 control sample size (default 50)\n")
    cat("    ratio_by_arm:    randomization proportion for control and treatment \n")
    cat("                     (default c(1, 1))\n")
    cat("    ratio_by_stage:  randomization proportion for each stage \n")
    cat("                     (default c(1, 1, 1))\n")
    cat("    p_by_arm:        response rate for control and treatment \n")
    cat("                     (default c(0.4, 0.7))\n")
    cat("    prior_by_arm:    list of beta priors by arm or list of normal prior for treatment difference\n")
    cat("                     (default list(c(1,1), c(1,1)) \n")
    cat("    n_post:          number of posterior samples \n")
    cat("                     (default 10000)\n")
    cat("    test_direction:  difference or logOR greater than or less than threshold \n")
    cat("                     (default 'greater than', can be 'less than')\n")
    cat("    decision_h0:     Threshold under null hypothesis \n")
    cat("                     (default 0.2)\n")
    cat("    decision_suc:    threshold of success in each stage \n")
    cat("                     (default c(0.9,0.9,0.9)\n")
    cat("    decision_fut:    threshold of futility in each stage \n")
    cat("                     (default c(0.5,0.4,0.3)\n")
    cat("    borrow:          borrow historical control or historical log OR or difference between treatment and control \n")
    cat("                     (default 'control', can be 'diff' or 'logOR') \n")
    cat("    test_statistic:  test statistic when borrowing historical control \n")
    cat("                     (default 'diff', can be 'logOR')\n")
}


internal_bayes_mstage_bin_dpara <- function() {
    list(sample_size     = 50,
         ratio_by_arm    = c(1, 1),
         ratio_by_stage  = c(1, 1, 1),
         p_by_arm        = c(0.4, 0.7),
         prior_by_arm    = list(c(1,1), c(1,1)),
         n_post          = 10000,
         test_direction  = 'greater than',
         decision_h0     = 0.2,
         decision_suc    = c(0.9,0.9,0.9),
         decision_fut    = c(0.5,0.4,0.3),
         borrow          = 'control',
         test_statistic  = 'diff'
         )
}

#' @export
bayes_mstage_gen_data_bin = function(lst_design, seed=NULL){

    if (!is.null(seed))
        old_seed = set.seed(seed)

                                        # if(sum(lst_design$ratio_by_arm)!=1 | sum(lst_design$ratio_by_stage)!=1){
                                        #   stop('Proportions does not sum to 1 !')
                                        # }

    if( (length(lst_design$ratio_by_stage) != length(lst_design$decision_suc)) |
        (length(lst_design$ratio_by_stage) != length(lst_design$decision_fut)) |
        (length(lst_design$decision_suc) != length(lst_design$decision_fut))
       ){
        stop('Number of stages and number of success/futility criteria do not correspond! ')
    }

    ssc1 = lst_design$sample_size
    rt_arm = lst_design$ratio_by_arm[2] / lst_design$ratio_by_arm[1]
    rt_stage = lst_design$ratio_by_stage / lst_design$ratio_by_stage[1]

    ssc = ceiling(ssc1*rt_stage)
    sst = ceiling(ssc*rt_arm)
    rps_sim_c = sapply(ssc, function(ssc){rbinom(1, ssc, p=lst_design$p_by_arm[1])})
    rps_sim_t = sapply(sst, function(sst){rbinom(1, sst, p=lst_design$p_by_arm[2])})

    n_stage = length(rt_stage)

    rst = data.frame(stage = rep(seq_len(n_stage), 2),
                     arm = c(rep('control', n_stage), rep('treatment', n_stage)),
                     n = c(ssc, sst),
                     y = c(rps_sim_c, rps_sim_t),
                     p = c(rep(lst_design$p_by_arm[1], n_stage), rep(lst_design$p_by_arm[2], n_stage))
                     )

    ## reset
    if (!is.null(seed))
        set.seed(old_seed)

    return(rst)
}

#' @export
get_postsample_bin_m = function(lst_design, ...){

    gendat = bayes_mstage_gen_data_bin(lst_design = lst_design, ...)

    ssc = (gendat$n)[seq_len(max(gendat$stage))]
    sst = (gendat$n)[(max(gendat$stage)+1):length(gendat$stage)]

    rpsc = (gendat$y)[seq_len(max(gendat$stage))]
    rpst = (gendat$y)[(max(gendat$stage)+1):length(gendat$stage)]

    lst_design = convert_prior_format_bin(lst_design)

    if (lst_design$borrow == 'diff'){

        get_postd = function(pc_hat, pt_hat, prior, ssc, sst){
            postd = postmix(prior, m=pt_hat-pc_hat,
                            se=max(sqrt(pc_hat*(1-pc_hat)/ssc + pt_hat*(1-pt_hat)/sst),0.0001))
            return(postd)
        }

        pc_hat = (gendat$y/gendat$n)[seq_len(max(gendat$stage))]
        pt_hat = (gendat$y/gendat$n)[(max(gendat$stage)+1):length(gendat$stage)]

        if (sum(pc_hat==0 | pt_hat==0) > 0){
            warning('Problematic asymptotic approximation due to small sample size or response rate !!!')
        }

        sample_lst = list()
        postd = lst_design$prior_by_arm[[1]]

        for(m in seq_len(max(gendat$stage))){

            postd = get_postd(pc_hat[m], pt_hat[m],
                              postd,
                              ssc[m], sst[m])

            k = sample.int(ncol(postd), lst_design$n_post, replace = TRUE, prob = postd[1,])
            sample_postd = rnorm(lst_design$n_post, postd[2, k], postd[3, k])
            sample_lst[m] = list(sample_postd)
        }
    }

    if (lst_design$borrow == 'logOR'){

        get_postd = function(a, b, c, d, prior){

            postd = postmix(prior, m=log(a*d/(b*c)),
                            se=sqrt(1/a+1/b+1/c+1/d))

            return(postd)
        }

        b = gendat$y[seq_len(max(gendat$stage))]
        a = gendat$y[(max(gendat$stage)+1):length(gendat$stage)]
        d = gendat$n[seq_len(max(gendat$stage))] - gendat$y[seq_len(max(gendat$stage))]
        c = gendat$n[(max(gendat$stage)+1):length(gendat$stage)] - gendat$y[(max(gendat$stage)+1):length(gendat$stage)]

        if(sum(a==0 | b==0| c==0 | d==0) > 0){
            warning('Problematic asymptotic approximation due to small sample size or response rate !!!')
        }

        sample_lst = list()
        postd = lst_design$prior_by_arm[[1]]

        for(m in seq_len(max(gendat$stage))){

            postd = get_postd(max(a[m], 0.0001), max(b[m], 0.0001),
                              max(c[m], 0.0001), max(d[m], 0.0001),
                              postd)

            k = sample.int(ncol(postd), lst_design$n_post, replace = TRUE, prob = postd[1,])
            sample_postd = rnorm(lst_design$n_post, postd[2, k], postd[3, k])
            sample_lst[m] = list(sample_postd)
        }
    }

    if (lst_design$borrow == 'control'){

        sample_lst = list()
        postd_c = lst_design$prior_by_arm[[1]]
        postd_t = lst_design$prior_by_arm[[2]]
        for(m in seq_len(max(gendat$stage))){
                                        #control
            postd_c = postmix(postd_c, n = ssc[m], r=rpsc[m])

            k = sample.int(ncol(postd_c),
                           lst_design$n_post,
                           replace = TRUE,
                           prob = postd_c[1,])

            sample_post_c = rbeta(lst_design$n_post,
                                  postd_c[2, k],
                                  postd_c[3, k])
                                        #treatment
            postd_t = postmix(postd_t,
                              n = sst[m],
                              r=rpst[m])

            k = sample.int(ncol(postd_t),
                           lst_design$n_post,
                           replace = TRUE,
                           prob = postd_t[1,])

            sample_post_t = rbeta(lst_design$n_post,
                                  postd_t[2, k],
                                  postd_t[3, k])

            ##diff
            if (lst_design$test_statistic == 'diff'){
                sample_postd = sample_post_t - sample_post_c
                sample_lst[m] = list(sample_postd)
            }

            if ((lst_design$test_statistic == 'logOR')){
                sample_postd   = log(sample_post_t/ (1-sample_post_t) / (sample_post_c/(1-sample_post_c)))
                sample_lst[m] = list(sample_postd)
            }
        }
    }

    return(sample_lst)
}

#' @export
get_postprobs_bin_m = function(lst_design, ...) {

    sample_lst  = get_postsample_bin_m(lst_design, ...)
    decision_df = get_trial_decision_from_smplst(lst_design, sample_lst)

    return(list(decision_df, sample_lst))
}
