
#' Parametric approximation of robust meta-analytic predictive prior for multiple historical control data
#'
#' 
rmap_control_bin = function(y,
                            N,
                            v = 1,
                            m = 0,
                            s = 1,
                            weight_noninfo = 0.5,
                            prior1_noninfo,
                            prior2_noninfo){
  
  map_mcmc_c = RBesT::gMAP(cbind(y,N-y)~1, 
                           data=data.frame(y,N), 
                           family='binomial', 
                           tau.dist="HalfNormal",
                           tau.prior=v, 
                           beta.prior=cbind(m,s))
  
  map_para_c = RBesT::automixfit(map_mcmc_c)
  print(plot(map_para_c))
  noninfo = RBesT::mixbeta(c(1, prior1_noninfo, prior2_noninfo))
  map_robust_c= RBesT::mixcombine(map_para_c, noninfo, weight = c(1-weight_noninfo, weight_noninfo))
  #map_robust_t = noninfo
  return(map_robust_c)
}


rmap_control_norm = function(y.m,
                             y.se,
                             v = 1,
                             m = 0,
                             s = 1,
                             weight_noninfo = 0.5,
                             prior1_noninfo,
                             prior2_noninfo){

  map_mcmc_c = RBesT::gMAP(cbind(y.m, y.se)~1, 
                           data=data.frame(y.m, y.se), 
                           family='gaussian', 
                           tau.dist="HalfNormal",
                           tau.prior=v, 
                           beta.prior=cbind(m,s))
  
  map_para_c = RBesT::automixfit(map_mcmc_c)
  print(plot(map_para_c))
  noninfo = RBesT::mixnorm(c(1, prior1_noninfo, prior2_noninfo))
  map_robust_c= RBesT::mixcombine(map_para_c, noninfo, weight = c(1-weight_noninfo, weight_noninfo))
  #map_robust_t = noninfo
  return(map_robust_c)
}

#' Parametric approximation of robust meta-analytic predictive prior for multiple historical differences between treatment and control
#'
#'
rmap_diff = function(y.dif,
                     y.dif.se,
                     v = 1,
                     m = 0,
                     s = 1,
                     weight_noninfo = 0.5,
                     prior1_noninfo,
                     prior2_noninfo){
  
  map_mcmc_diff = RBesT::gMAP(cbind(y.dif, y.dif.se)~1, 
                              data=data.frame(y.dif, y.dif.se), 
                              family='gaussian', 
                              tau.dist="HalfNormal", 
                              tau.prior=v, 
                              beta.prior=cbind(m,s))
  
  map_para_diff = RBesT::automixfit(map_mcmc_diff)
  print(plot(map_para_diff))
  noninfo = RBesT:: mixnorm(c(1, prior1_noninfo, prior2_noninfo))
  map_robust_diff= RBesT::mixcombine(map_para_diff, noninfo, weight = c(1-weight_noninfo, weight_noninfo))
  return(map_robust_diff)
}


rmap_logOR = function(logOR,
                      logOR.se,
                      v = 1,
                      m = 0,
                      s = 1,
                      weight_noninfo = 0.5,
                      prior1_noninfo,
                      prior2_noninfo){
  
  map_mcmc_diff = RBesT::gMAP(cbind(logOR, logOR.se)~1, 
                              data=data.frame(logOR, logOR.se), 
                              family='gaussian', 
                              tau.dist="HalfNormal", 
                              tau.prior=v, 
                              beta.prior=cbind(m,s))
  
  map_para_diff = RBesT::automixfit(map_mcmc_diff)
  print(plot(map_para_diff))
  noninfo = RBesT:: mixnorm(c(1, prior1_noninfo, prior2_noninfo))
  map_robust_diff= RBesT::mixcombine(map_para_diff, noninfo, weight = c(1-weight_noninfo, weight_noninfo))
  return(map_robust_diff)
}


#' Parametric approximation of robust meta-analytic predictive prior
#'
#' @export
#' 
stb_tl_bayes_para_rmap = function(dtype = c('binary', 'normal'), 
                                  borrow = c('control', 'diff', 'logOR'), 
                                  cores = 1, ...){
  
  borrow = match.arg(borrow)
  dtype = match.arg(dtype)
  
  if(borrow == 'control'){
    if(dtype == 'binary'){
      rMAP = rmap_control_bin(...)
    }
    if(dtype == 'normal'){
      rMAP = rmap_control_norm(...)
    }
  }
  if(borrow == 'diff'){
    rMAP = rmap_diff(...)
  }
  if(borrow == 'logOR'){
    rMAP = rmap_logOR(...)
  }
  return(rMAP)
}