postmix = function(...){
  suppressMessages(RBesT::postmix(...))
}

#' @export
get_trial_decision_from_smplst = function(lst_design, sample_lst){

  decision_df = NULL
  indeterm = TRUE

  for(i in 1:length(sample_lst)){

    if(!indeterm){
      temp_row = c(success, futility, indeterm)
      names(temp_row) = paste0(c('success_stage', 'futility_stage', 'indeterm_stage'), i)
      decision_df = c(decision_df, temp_row)
    }

    if(indeterm){

      sample_postd = sample_lst[[i]]
      if(lst_design$test_direction == 'greater than'){
        postprob = sum(sample_postd > lst_design$decision_h0)/lst_design$n_post
      }
      if(lst_design$test_direction == 'less than'){
        postprob = sum(sample_postd < lst_design$decision_h0)/lst_design$n_post
      }

      success = postprob > lst_design$decision_suc[i]
      futility = postprob < lst_design$decision_fut[i]
      indeterm = !(postprob > lst_design$decision_suc[i]) & !(postprob < lst_design$decision_fut[i])

      temp_row = c(success, futility, indeterm)
      names(temp_row) = paste0(c('success_stage', 'futility_stage', 'indeterm_stage'), i)
      decision_df = c(decision_df, temp_row)
    }

  }
  return(decision_df)
}



#' @export

stb_tl_get_oc = function(x, ss_seq, n_rep, seed, n_core,...){
  oc_tab = c()
  for(ss in ss_seq){
    x@design_para$sample_size = ss
    y = stb_create_trial(x=x, seed=seed, ...)
    probs = unlist(stb_create_simustudy(x=x, n_rep=n_rep, seed=seed, n_core=n_core, ...)@rst_key)

    ssc = (y@data$n)[seq_len(max(y@data$stage))]
    sst = (y@data$n)[(max(y@data$stage)+1):length(y@data$stage)]
    ss_total = sum(ssc) + sum(sst)
    temp_row = c(ss_total, ssc, sst, probs)

    names(temp_row) = c('total sample size ',
                        sapply(1:max(y@data$stage),
                               function(i){paste('sample size control stage', i)}),
                        sapply(1:max(y@data$stage),
                               function(i){paste('sample size treatment stage', i)}),
                        sapply(1:max(y@data$stage),
                               function(i){paste(c('success rate stage', 'futility rate stage', 'indeterminate rate stage'), i)})
    )

    oc_tab = rbind(oc_tab, temp_row)
  }

  oc_tab = t(oc_tab)
  colnames(oc_tab) = sapply(1:ncol(oc_tab), function(i){paste0('scenario ', i)})

  plt = stb_tl_plot_oc(x, oc_tab)
  print(plt)

  return(oc_tab)
}


stb_tl_plot_oc = function(x, oc_tab){
  plt_lst = list()
  for(i in 1:length(x@design_para$ratio_by_stage)){
    pltdat = data.frame(sample_size = oc_tab[1,],
                        success_rate = oc_tab[paste('success rate stage', i),],
                        futility_rate = oc_tab[paste('futility rate stage', i),],
                        indeterminate_rate = oc_tab[paste('indeterminate rate stage', i),])
    plt_stage = ggplot2::ggplot(pltdat, aes(x=sample_size))+
      ggplot2::geom_line(aes(y=success_rate, color='success'), linewidth=2) +
      ggplot2::geom_line(aes(y=futility_rate, color='futility'), linewidth =2) +
      ggplot2::geom_line(aes(y=indeterminate_rate, color='indeterminate'), linewidth=2) +
      ggplot2::labs(title = paste('stage',i), x = 'total sample size', y = 'cumulative probability', color = 'decision') +
      ggplot2::scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0,1))

    plt_lst[[i]] = plt_stage
    #print(plt_stage)
  }

  plt = do.call(ggpubr::ggarrange, c(plt_lst, nrow=1, common.legend=TRUE, legend = 'top'))
  ggpubr::annotate_figure(plt, top = ggpubr::text_grob(paste0('sample size proportion per arm = ',
                                                               paste(round(x@design_para$ratio_by_arm,2), collapse = ' : '),
                                                               '\n',
                                                               'sample size proportion per stage = ',
                                                               paste(round(x@design_para$ratio_by_stage,2), collapse = ' : ')))
  )
}

#' @export
stb_tl_get_bayes_update_plt_bin = function(lst_design,
                                           data,
                                           rst, seed=NULL, ...){

    lst_design = convert_prior_format_bin(lst_design)

    if (lst_design$borrow =='control'){
        prior_c = lst_design$prior_by_arm[[1]]
        prior_t = lst_design$prior_by_arm[[2]]
        set.seed(seed)
        k = sample.int(ncol(prior_c), lst_design$n_post ,
                       replace = TRUE, prob = prior_c[1,])
        sample_prior_c = rbeta(lst_design$n_post, prior_c[2, k], prior_c[3, k])

        k = sample.int(ncol(prior_t), lst_design$n_post , replace = TRUE, prob = prior_t[1,])
        sample_prior_t = rbeta(lst_design$n_post, prior_t[2, k], prior_t[3, k])
                                        #diff_prior
        sample_prior = sample_prior_t - sample_prior_c
    }

    if(lst_design$borrow =='diff' | lst_design$borrow =='logOR' ){

        prior = lst_design$prior_by_arm[[1]]
        set.seed(seed)
        k = sample.int(ncol(prior), lst_design$n_post, replace = TRUE, prob = prior[1,])
        sample_prior = rnorm(lst_design$n_post, prior[2, k], prior[3, k])
    }

    sample_list = rst[[2]]


    pltdat = data.frame(theta = c(sample_prior, unlist(sample_list)),
                        stage=factor(c(rep('prior', lst_design$n_post),
                                       rep(1:length(sample_list), each = lst_design$n_post))))

    ggplot(pltdat, aes(x=theta, color=stage, fill=stage)) +
        geom_density(linewidth=1, alpha=0.3) +
        ggtitle(paste0('density plot for posterior samples, ', 'total sample size = ', sum(data$n)))
}



#' @export
stb_tl_get_bayes_update_plt_norm = function(lst_design, data, rst, seed=NULL, ...){

  lst_design = convert_prior_format_norm(lst_design)

  if(lst_design$borrow =='control'){
    prior_c =lst_design$prior_by_arm[[1]]
    prior_t = lst_design$prior_by_arm[[2]]
    set.seed(seed)
    k = sample.int(ncol(prior_c), lst_design$n_post , replace = TRUE, prob = prior_c[1,])
    sample_prior_c = rnorm(lst_design$n_post, prior_c[2, k], prior_c[3, k])

    k = sample.int(ncol(prior_t), lst_design$n_post , replace = TRUE, prob = prior_t[1,])
    sample_prior_t = rnorm(lst_design$n_post, prior_t[2, k], prior_t[3, k])

    sample_prior = sample_prior_t - sample_prior_c
  }

  if(lst_design$borrow =='diff'){

    prior = lst_design$prior_by_arm[[1]]
    set.seed(seed)
    k = sample.int(ncol(prior), lst_design$n_post, replace = TRUE, prob = prior[1,])
    sample_prior = rnorm(lst_design$n_post, prior[2, k], prior[3, k])
  }

  sample_list = rst[[2]]

  pltdat = data.frame(theta = c(sample_prior, unlist(sample_list)),
                      stage=factor(c(rep('prior', lst_design$n_post),
                                     rep(1:length(sample_list), each = lst_design$n_post))))

  ggplot(pltdat, aes(x=theta, color=stage, fill=stage)) +
    geom_density(linewidth=1, alpha=0.3) +
    ggtitle(paste0('density plot for posterior samples, ', 'total sample size = ', sum(data$n)))
  #print(plt)
}


convert_prior_format_norm = function(lst_design){
  if(lst_design$borrow == 'diff'){
    if(class(lst_design$prior_by_arm[[1]])[1] =='numeric'){
    lst_design$prior_by_arm[[1]] = RBesT::mixnorm(c(1,lst_design$prior_by_arm[[1]]))
    }
  }
  if(lst_design$borrow == 'control'){
    lst_design$prior_by_arm[[2]] = RBesT::mixnorm(c(1,lst_design$prior_by_arm[[2]]))
    if(class(lst_design$prior_by_arm[[1]])[1] =='numeric'){
      lst_design$prior_by_arm[[1]] = RBesT::mixnorm(c(1,lst_design$prior_by_arm[[1]]))
    }
  }
  return(lst_design)
}


convert_prior_format_bin = function(lst_design){
  if(lst_design$borrow == 'diff'){
    if(class(lst_design$prior_by_arm[[1]])[1] =='numeric'){
      lst_design$prior_by_arm[[1]] = RBesT::mixnorm(c(1,lst_design$prior_by_arm[[1]]))
    }
  }
  if(lst_design$borrow == 'control'){
    lst_design$prior_by_arm[[2]] = RBesT::mixbeta(c(1,lst_design$prior_by_arm[[2]]))
    if(class(lst_design$prior_by_arm[[1]])[1] =='numeric'){
      lst_design$prior_by_arm[[1]] = RBesT::mixbeta(c(1,lst_design$prior_by_arm[[1]]))
    }
  }
  return(lst_design)
}
