## -----------------------------------------------------------------------------
##
##  DESCRIPTION:
##      Bayesian Design for Normal Endpoint
##
##  DATE:
##      July, 2024
## -----------------------------------------------------------------------------

#' Describe the design
#'
#'
bayes_mstage_norm_describe <- function(x, ...) {
  cat("Type:\n")
  cat("    Bayesian multi-stage design for normal outcomes with historical borrowing.\n\n")
  
  cat("Design Parameters:\n")
  cat("    sample_size:     stage 1 control sample size (default 50) \n")
  cat("    ratio_by_arm:    randomization proportion for control and treatment \n")
  cat("                     (default c(1, 1))\n")
  cat("    ratio_by_stage:  randomization proportion for each stage \n")
  cat("                     (default c(1, 1, 1))\n")
  cat("    m_s_control:     mean and sd for control arm\n")
  cat("                     (default c(10, 3))\n")
  cat("    m_s_treat:       mean and sd for treatent arm\n")
  cat("                     (default c(15, 3))\n")
  cat("    prior_by_arm:    list of normal priors by arm or list of normal prior for treatment difference\n")
  cat("                     (default list(c(0,1), c(0,1)) \n")
  cat("    n_post:          number of posterior samples \n")
  cat("                     (default 10000)\n")
  cat("    test_direction:  difference greater than or less than threshold \n")
  cat("                     (default 'greater than', can be 'less than')\n")
  cat("    decision_h0:     Threshold under null hypothesis \n")
  cat("                     (default 0.2)\n")
  cat("    decision_suc:    Threshold of success in each stage \n")
  cat("                     (default c(0.9,0.9,0.9))\n")
  cat("    decision_fut:    Threshold of futility in each stage \n")
  cat("                     (default c(0.5,0.4,0.3))\n")
  cat("    borrow:          borrow historical control only or difference between treatment and control \n")
  cat("                     (default 'control', can be 'diff')\n")
}


internal_bayes_mstage_norm_dpara <- function() {
  list(sample_size     = 50,
       ratio_by_arm    = c(1, 1),
       ratio_by_stage  = c(1, 1, 1), 
       m_s_control     = c(10, 3),
       m_s_treat       = c(15, 3),
       prior_by_arm    = list(c(0,1), c(0,1)),
       n_post          = 10000,
       test_direction  = 'greater than',
       decision_h0     = 4,
       decision_suc    = c(0.9,0.9,0.9),
       decision_fut    = c(0.5,0.4,0.3),
       borrow          = 'control'
  )
}


#' @export
bayes_mstage_gen_data_norm <- function(lst_design, seed=NULL){
  
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
  rps_sim_c = lapply(ssc, function(ssc){rnorm(ssc, lst_design$m_s_control[1], lst_design$m_s_control[2])})
  rps_sim_t = lapply(sst, function(sst){rnorm(sst, lst_design$m_s_treat[1], lst_design$m_s_treat[2])})
  
  n_stage = length(rt_stage)
  
  rst = data.frame(stage = rep(seq_len(n_stage), 2), 
                   arm = c(rep('control', n_stage), rep('treatment', n_stage)),
                   n = c(ssc, sst),
                   y.m = c(sapply(1:n_stage, function(i){mean(rps_sim_c[[i]])}), 
                           sapply(1:n_stage, function(i){mean(rps_sim_t[[i]])})),
                   y.sd = c(sapply(1:n_stage, function(i){sd(rps_sim_c[[i]])}), 
                            sapply(1:n_stage, function(i){sd(rps_sim_t[[i]])}))
  )
  
  ## reset
  if (!is.null(seed))
    set.seed(old_seed)
  
  return(rst)
}



#' @export
get_postsample_norm = function(lst_design, ...){
  
  gendat = bayes_mstage_gen_data_norm(lst_design = lst_design, ...)
  
  lst_design = convert_prior_format_norm(lst_design)
  
  ssc = (gendat$n)[seq_len(max(gendat$stage))]
  sst = (gendat$n)[(max(gendat$stage)+1):length(gendat$stage)]
  
  m.c = (gendat$y.m)[seq_len(max(gendat$stage))]
  m.t = (gendat$y.m)[(max(gendat$stage)+1):length(gendat$stage)]
  
  sd.c = (gendat$y.sd)[seq_len(max(gendat$stage))]
  sd.t = (gendat$y.sd)[(max(gendat$stage)+1):length(gendat$stage)]
  
  if (lst_design$borrow == 'diff'){
    
    get_postd = function(m.c, sd.c, m.t, sd.t, prior, ssc, sst){
      postd = postmix(prior, m=m.t-m.c, se=sqrt(sd.c**2/ssc + sd.t**2/sst))
      return(postd)
    }
    
    sample_lst = list()
    postd = lst_design$prior_by_arm[[1]]
    for(m in seq_len(max(gendat$stage))){
      postd = get_postd(m.c[m], sd.c[m], m.t[m], sd.t[m], postd, ssc[m], sst[m])
      k = sample.int(ncol(postd), lst_design$n_post, replace = TRUE, prob = postd[1,])
      sample_postd = rnorm(lst_design$n_post, postd[2, k], postd[3, k])
      sample_lst = append(sample_lst, list(sample_postd))
    }
  }
  
  if (lst_design$borrow == 'control'){
    
    sample_lst = list()
    postd_c = lst_design$prior_by_arm[[1]]
    postd_t = lst_design$prior_by_arm[[2]]
    for(m in seq_len(max(gendat$stage))){
      #control
      postd_c = postmix(postd_c, m=m.c[m], se=sd.c[m]/sqrt(ssc[m]))
      k = sample.int(ncol(postd_c), lst_design$n_post, replace = TRUE, prob = postd_c[1,])
      sample_post_c = rnorm(lst_design$n_post, postd_c[2, k], postd_c[3, k])
      #treatment
      postd_t = postmix(postd_t, m=m.t[m], se=sd.t[m]/sqrt(sst[m]))
      k = sample.int(ncol(postd_t), lst_design$n_post, replace = TRUE, prob = postd_t[1,])
      sample_post_t = rnorm(lst_design$n_post, postd_t[2, k], postd_t[3, k])
      #diff
      sample_postd = sample_post_t - sample_post_c
      sample_lst = append(sample_lst, list(sample_postd))
    }
  }
  return(sample_lst)
}


#' @export
get_postprobs_norm = function(lst_design, ...){
  
  sample_lst = get_postsample_norm(lst_design, ...)
  
  decision_df = get_trial_decision_from_smplst(lst_design, sample_lst)
  
  return(list(decision_df, sample_lst))
}

