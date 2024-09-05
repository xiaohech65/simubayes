cl = parallel::makeCluster(detectCores()-2)
doParallel::registerDoParallel(cl)

gen_raw = function(lst) {
  rst = NULL
  for (i in seq_len(length(lst))) {
    temp_row = lst[[i]][[1]]
    rst = rbind(rst, temp_row)
  }
  return(list(rst))
}

gen_summary = function(lst) {
  df = as.data.frame(lst[[1]])
  len = nrow(df)
  return(list(colSums(df, na.rm=TRUE)/len))
}

stb_tl_get_bayes_update_shiny = function(lst_design, dtype, ss_seq){
  if(dtype == 'binary'){
    lst_design$sample_size = min(ss_seq)
    gendat = bayes_mstage_gen_data_bin(lst_design)
    rst = get_postprobs_bin_m(lst_design)
    plt1 = stb_tl_get_bayes_update_plt_bin(lst_design, gendat, rst, seed=1)

    lst_design$sample_size = max(ss_seq)
    gendat = bayes_mstage_gen_data_bin(lst_design)
    rst = get_postprobs_bin_m(lst_design)
    plt2 = stb_tl_get_bayes_update_plt_bin(lst_design, gendat, rst, seed=1)
  }
  if(dtype == 'normal'){
    lst_design$sample_size = min(ss_seq)
    gendat = bayes_mstage_gen_data_norm(lst_design)
    rst = get_postprobs_norm(lst_design)
    plt1 = stb_tl_get_bayes_update_plt_norm(lst_design, gendat, rst, seed=1)

    lst_design$sample_size = max(ss_seq)
    gendat = bayes_mstage_gen_data_norm(lst_design)
    rst = get_postprobs_norm(lst_design)
    plt2 = stb_tl_get_bayes_update_plt_norm(lst_design, gendat, rst, seed=1)
  }
  plt_lst = list(plt1, plt2)
  do.call(ggpubr::ggarrange, c(plt_lst, nrow=1, common.legend=TRUE, legend = 'top'))
}

stb_tl_get_oc_shiny = function(lst_design, dtype, ss_seq, n_rep){
  oc_tab = c()
  for(ss in ss_seq){
    lst_design$sample_size = ss
    if(dtype == 'binary'){
      gendat = bayes_mstage_gen_data_bin(lst_design)
      rst = foreach(i=1:n_rep,
                    .export = ls(globalenv()),
                    .packages=c('RBesT', 'simubayes2'))  %dopar%{
                      get_postprobs_bin_m(lst_design)
                    }
    }
    if(dtype == 'normal'){
      gendat = bayes_mstage_gen_data_norm(lst_design)
      rst = foreach(i=1:n_rep,
                   .export = ls(globalenv()),
                   .packages=c('RBesT', 'simubayes2')) %dopar%{
                     get_postprobs_norm(lst_design)
                   }
    }
    rst = gen_raw(rst)
    probs = gen_summary(rst)[[1]]

    ssc = (gendat$n)[seq_len(max(gendat$stage))]
    sst = (gendat$n)[(max(gendat$stage)+1):length(gendat$stage)]
    ss_total = sum(ssc) + sum(sst)
    temp_row = c(ss_total, ssc, sst, probs)

    names(temp_row) = c('total sample size ',
                        sapply(1:max(gendat$stage),
                               function(i){paste('sample size control stage', i)}),
                        sapply(1:max(gendat$stage),
                               function(i){paste('sample size treatment stage', i)}),
                        sapply(1:max(gendat$stage),
                               function(i){paste(c('success rate stage', 'futility rate stage', 'indeterminate rate stage'), i)})
    )

    oc_tab = rbind(oc_tab, temp_row)
  }

  oc_tab = t(oc_tab)
  colnames(oc_tab) = sapply(1:ncol(oc_tab), function(i){paste0('scenario ', i)})

  return(oc_tab)
}


stb_tl_plot_oc_shiny = function(lst_deisign, oc_tab){
  plt_lst = list()
  for(i in 1:length(lst_deisign$ratio_by_stage)){
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
  }

  plt = do.call(ggpubr::ggarrange, c(plt_lst, nrow=1, common.legend=TRUE, legend = 'top'))
  ggpubr::annotate_figure(plt, top = ggpubr::text_grob(paste0('sample size proportion per arm = ',
                                                              paste(round(lst_deisign$ratio_by_arm,2), collapse = ' : '),
                                                              '\n',
                                                              'sample size proportion per stage = ',
                                                              paste(round(lst_deisign$ratio_by_stage,2), collapse = ' : ')))
  )
}


stb_tl_plot_type1_shiny = function(lst_deisign, oc_tab){
  plt_lst = list()
  for(i in 1:length(lst_deisign$ratio_by_stage)){
    pltdat = data.frame(sample_size = oc_tab[1,],
                        type1_error = oc_tab[paste('success rate stage', i),])
    plt_stage = ggplot2::ggplot(pltdat, aes(x=sample_size))+
      ggplot2::geom_line(aes(y=type1_error, color='type 1 error rate'), linewidth=2) +
      ggplot2::labs(title = paste('stage',i), x = 'total sample size', y = '', color = '') +
      ggplot2::scale_y_continuous(breaks = seq(0, 0.5, 0.05), limits = c(0, 0.5))

    plt_lst[[i]] = plt_stage
  }

  plt = do.call(ggpubr::ggarrange, c(plt_lst, nrow=1, common.legend=TRUE, legend = 'top'))
  ggpubr::annotate_figure(plt, top = ggpubr::text_grob(paste0('sample size proportion per arm = ',
                                                              paste(round(lst_deisign$ratio_by_arm,2), collapse = ' : '),
                                                              '\n',
                                                              'sample size proportion per stage = ',
                                                              paste(round(lst_deisign$ratio_by_stage,2), collapse = ' : ')))
  )
}
