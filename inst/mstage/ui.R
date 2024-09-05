library(shinybusy)
library(doParallel)
library(foreach)
library(simubayes2)
library(ggplot2)
library(ggpubr)

borrow_bin = c('control arm', 'treatment difference between arms', 'log OR')
borrow_norm = c('control arm', 'treatment difference between arms')
endpoint = c('binary','normal')
direction_input = c('greater than', 'less than')
prior = c('Power Prior', 'Meta-Analytic Predictive Prior')

parameter_tabs_bin = conditionalPanel(
  condition = "input.outcome=='binary'",
  splitLayout(
    numericInput('pi.c', withMathJax("$$ \\pi_c$$"), value=0.3, min=0, max=1),
    numericInput('pi.t', withMathJax("$$ \\pi_t$$"), value=0.6, min=0, max=1)
  )
)

parameter_tabs_norm = conditionalPanel(
  condition = "input.outcome=='normal'",
  splitLayout(
    numericInput("mu_c", withMathJax("$$\\mu_c $$"), value = 4),
    numericInput("sd_c", withMathJax("$$\\sigma_c $$"), value = 2, min=0.01)
  ),
  splitLayout(
    numericInput("mu_t", withMathJax("$$\\mu_t $$"), value = 6),
    numericInput("sd_t", withMathJax("$$\\sigma_t $$"), value = 2, min=0.01)
  )
)



MAP_tabs = conditionalPanel(
  condition = "input.prior=='Meta-Analytic Predictive Prior'",
  conditionalPanel(
    condition = "input.outcome == 'binary'",
    conditionalPanel(
      condition = "input.borrow=='control arm'",
      withMathJax("$$ y_i|\\pi_i \\sim \\text{Binomial}(N_i, \\pi_i), i=1,...,H$$"),
      withMathJax("$$\\pi_i \\sim \\text{LogitNormal}(\\mu, \\tau)$$"),
      withMathJax("$$ \\mu \\sim \\text{Normal}(m, s)$$"),
      withMathJax("$$ \\tau \\sim \\text{HalfNormal}(v)$$")
    ),
    conditionalPanel(
      condition = "input.borrow=='treatment difference between arms' | input.borrow == 'log OR'",
      withMathJax("$$ d_i|\\theta_{d_i} \\sim \\text{Normal}(\\theta_{d_i}, \\sigma_{d_i}), i=1,...,H$$"),
      withMathJax("$$\\theta_{d_i} \\sim \\text{Normal}(\\mu, \\tau)$$"),
      withMathJax("$$ \\mu \\sim \\text{Normal}(m, s)$$"),
      withMathJax("$$ \\tau \\sim \\text{HalfNormal}(v)$$")
    )
  ),
  conditionalPanel(
    condition = "input.outcome == 'normal'",
    conditionalPanel(
      condition = "input.borrow=='control arm'",
      withMathJax("$$ \\bar{y}_i|\\theta_i \\sim \\text{Normal}(\\theta_i, \\sigma_i), i=1,...,H$$"),
      withMathJax("$$\\theta_i \\sim \\text{Normal}(\\mu, \\tau)$$"),
      withMathJax("$$ \\mu \\sim \\text{Normal}(m, s)$$"),
      withMathJax("$$ \\tau \\sim \\text{HalfNormal}(v)$$")
    ),
    conditionalPanel(
      condition = "input.borrow=='treatment difference between arms'",
      withMathJax("$$ d_i|\\theta_{d_i} \\sim \\text{Normal}(\\theta_{d_i}, \\sigma_{d_i}), i=1,...,H$$"),
      withMathJax("$$\\theta_{d_i} \\sim \\text{Normal}(\\mu, \\tau)$$"),
      withMathJax("$$ \\mu \\sim \\text{Normal}(m, s)$$"),
      withMathJax("$$ \\tau \\sim \\text{HalfNormal}(v)$$")
    )
  ),
  conditionalPanel(
    condition = "input.borrow=='control arm'",

    conditionalPanel(
      condition="input.outcome=='normal'",
      textInput('y.m', 'historical control effect (y.m) (comma delimited)', value=5),
      textInput('y.se', 'standard error of historical control effect (y.se) (comma delimited)', value=1),
      numericInput('v','half-normal tau prior (v)', value=1, min=0.01),
      numericInput('m','normal mu prior (m)', value=0),
      numericInput('s','normal mu prior (s)', value=1, min=0.01),
      sliderInput('w','weight of non-informative prior', value=0.5, min=0, max=1),
      numericInput('p1.noninfo.n', 'normal noninformative prior to mix with (mean)', value=0),
      numericInput('p2.noninfo.n', 'normal noninformative prior to mix with (sd)', value=1, min=0.01),
      numericInput('t.a.n', 'normal prior for treatment arm (mu)', value = 0, min=0.01),
      numericInput('t.b.n', 'normal prior for treatment arm (tau)', value = 1, min=0.01)
    ),

    conditionalPanel(
      condition = "input.outcome=='binary'",
      textInput('y', 'historical control responses (y) (comma delimited)', value=20),
      textInput('N', 'historical control sample size (N) (comma delimited)', value=100),
      numericInput('v','half-normal tau prior (v)', value=1, min=0.01),
      numericInput('m','normal mu prior (m)', value=0),
      numericInput('s','normal mu prior (s)', value=1, min=0.01),
      sliderInput('w','weight of non-informative prior', value=0.5, min=0, max=1),
      numericInput('p1.noninfo', 'beta noninformative prior to mix with (a)', value=1),
      numericInput('p2.noninfo', 'beta noninformative prior to mix with (b)', value=1),
      numericInput('t.a', 'beta prior for treatment arm (a)', value = 1, min=0.01),
      numericInput('t.b', 'beta prior for treatment arm (b)', value = 1, min=0.01),

    )
  ),

  conditionalPanel(
    condition = "input.borrow=='treatment difference between arms'",
    textInput('y.dif', 'treatment difference (y.dif) ', value=0),
    textInput('y.dif.se', 'standard error of treatment difference (y.dif.se) (comma delimited)', value=1),
    numericInput('v','half-normal tau prior (v)', value=1, min=0.01),
    numericInput('m','normal mu prior (m)', value=0),
    numericInput('s','normal mu prior (s)', value=1, min=0.01),
    sliderInput('w','weight of non-informative prior', value=0.5, min=0, max=1),
    numericInput('p1.noninfo.n', 'normal noninformative prior to mix with (mean)', value=0),
    numericInput('p2.noninfo.n', 'normal noninformative prior to mix with (sd)', value=1, min=0.01)
  ),
  conditionalPanel(
    condition = "input.borrow=='log OR'",
    textInput('logOR', 'log OR (comma delimited)', value=0),
    textInput('logOR.se', 'standard error of log OR (comma delimited)', value=1),
    numericInput('v','half-normal tau prior (v)', value=1, min=0.01),
    numericInput('m','normal mu prior (m)', value=0),
    numericInput('s','normal mu prior (s)', value=1, min=0.01),
    sliderInput('w','weight of non-informative prior', value=0.5, min=0, max=1),
    numericInput('p1.noninfo.n', 'normal noninformative prior to mix with (mean)', value=0),
    numericInput('p2.noninfo.n', 'normal noninformative prior to mix with (sd)', value=1, min=0.01)
  )
)


PP_tabs = conditionalPanel(
  condition = "input.prior=='Power Prior'",
  conditionalPanel(
    condition = "input.outcome == 'binary'",
    conditionalPanel(
      condition = "input.borrow=='control arm'",
      withMathJax("$$ y_i | \\pi \\sim \\text{Binomial}(N_i, \\pi), i=1,...,H $$"),
      withMathJax("$$ \\pi \\sim \\text{Beta}(a, b) $$"),
      withMathJax("$$\\pi | D_0 \\sim \\prod_{i=1}^{H} L(\\pi|y_i, N_i)^k f(\\pi)$$")
    ),
    conditionalPanel(
      condition = "input.borrow=='treatment difference between arms' | input.borrow == 'log OR'",
      withMathJax("$$ d_i | \\theta_d \\sim \\text{Normal}(\\theta_d, \\sigma_i), i=1,...,H $$"),
      withMathJax("$$ \\theta_d \\sim \\text{Normal}(m, s) $$"),
      withMathJax("$$\\theta_d | D_0 \\sim \\prod_{i=1}^H L(\\theta_d | d_i, \\sigma_i)^k f(\\theta_d) $$")
    ),
  ),
  conditionalPanel(
    condition = "input.outcome == 'normal'",
    conditionalPanel(
      condition = "input.borrow=='control arm'",
      withMathJax("$$ \\bar{y}_i | \\theta \\sim \\text{Normal}(\\theta, \\sigma_i), i=1,...,H $$"),
      withMathJax("$$ \\theta \\sim \\text{Normal}(m, s) $$"),
      withMathJax("$$\\theta | D_0 \\sim \\prod_{i=1}^H L(\\theta |\\bar{y}_i, \\sigma_i)^k f(\\theta) $$")
    ),
    conditionalPanel(
      condition = "input.borrow=='treatment difference between arms'",
      withMathJax("$$ d_i | \\theta_d \\sim \\text{Normal}(\\theta_d, \\sigma_i), i=1,...,H $$"),
      withMathJax("$$ \\theta_d \\sim \\text{Normal}(m, s) $$"),
      withMathJax("$$\\theta_d | D_0 \\sim \\prod_{i=1}^H L(\\theta_d |d_i, \\sigma_i)^k f(\\theta_d) $$")
    )
  ),
  conditionalPanel(
    condition = "input.borrow=='control arm'",
    sliderInput('power','power parameter (k) ', value=0.3, min=0, max=1),
    conditionalPanel(
      condition = "input.outcome=='binary'",
      textInput('y', 'historical control responses (y) (comma delimited)', value=20),
      textInput('N', 'historical control sample sizes(N) (comma delimited)', value=100),
      numericInput('a','beta prior for control arm (a)', value=1, min=0.01),
      numericInput('b','beta prior for control arm (b)', value=1, min=0.01),
      numericInput('t.a', 'beta prior for treatment arm (a)', value = 1, min=0.01),
      numericInput('t.b', 'beta prior for treatment arm (b)', value = 1, min=0.01)
    ),
    conditionalPanel(
      condition = "input.outcome=='normal'",
      textInput('y.m', 'control mean (y.m) (comma delimited)', value=4),
      textInput('y.se', 'standard error of control mean (y.se) (comma delimited)', value=1),
      numericInput('a.n','normal prior for control arm (m)', value=0),
      numericInput('b.n','normal prior for control arm (s)', value=1, min=0.01),
      numericInput('t.a.n', 'normal prior for treatment arm (m)', value = 0, min=0.01),
      numericInput('t.b.n', 'normal prior for treatment arm (s)', value = 1, min=0.01)
    )
  ),
  conditionalPanel(
    condition = "input.borrow=='treatment difference between arms'",
    textInput('y.dif', 'treatment difference (y.dif) (comma delimited)', value=0),
    textInput('y.dif.se', 'standard error of treatment difference (y.dif.se) (comma delimited)', value=1),
    sliderInput('power','power parameter (k)', value=0.3, min=0, max=1),
    numericInput('a.n','normal prior (m)', value=0),
    numericInput('b.n','normal prior (s)', value=1, min=0.01)
  ),
  conditionalPanel(
    condition = "input.borrow=='log OR'",
    textInput('logOR', 'log OR (comma delimited)', value=0),
    textInput('logOR.se', 'standard error of log OR (comma delimited)', value=1),
    sliderInput('power','power parameter (k)', value=0.3, min=0, max=1),
    numericInput('a.n','normal prior (m)', value=0),
    numericInput('b.n','normal prior (s)', value=1, min=0.01)
  )
)

side_panel1 =  sidebarPanel(
  radioButtons('outcome','type of your clinical outcome',endpoint),
  conditionalPanel(
    condition = "input.outcome == 'binary'",
    radioButtons('borrow', 'borrow historical information from', borrow_bin)
  ),
  conditionalPanel(
    condition = "input.outcome == 'normal'",
    radioButtons('borrow', 'borrow historical information from', borrow_norm)
  ),
  radioButtons('prior', 'select type of prior', prior),
  MAP_tabs,
  PP_tabs
)


side_panel2 =   sidebarPanel(
  radioButtons('direct','test direction w.r.t threshold delta:',direction_input),
  conditionalPanel(
    condition = "input.outcome == 'binary' && input.borrow == 'control arm'",
    radioButtons('teststat', 'test statistic', c('diff', 'logOR'))
  ),
  conditionalPanel(
    condition = "input.direct == 'greater than' && input.outcome =='binary'",
    conditionalPanel(
      condition = "input.borrow == 'treatment difference between arms'",
      withMathJax("$$\\text{success criterion: }$$"),
      withMathJax("$$P(\\theta_d=\\pi_t-\\pi_c > \\delta | D_0, D) > c_s $$"),
      withMathJax("$$\\text{futility criterion: }$$"),
      withMathJax("$$P(\\theta_d=\\pi_t-\\pi_c > \\delta | D_0, D) < c_f $$")
    ),
    conditionalPanel(
      condition = "input.borrow == 'log OR'",
      withMathJax("$$\\text{success criterion: }$$"),
      withMathJax("$$P(\\theta_d = \\log \\frac{\\pi_t/(1-\\pi_t)}{\\pi_c/(1-\\pi_c)} > \\delta | D_0, D) > c_s $$"),
      withMathJax("$$\\text{futility criterion: }$$"),
      withMathJax("$$P(\\theta_d = \\log \\frac{\\pi_t/(1-\\pi_t)}{\\pi_c/(1-\\pi_c)} > \\delta | D_0, D) < c_f $$")
    ),
    conditionalPanel(
      condition = "input.borrow == 'control arm' && input.teststat == 'diff'",
      withMathJax("$$\\text{success criterion: }$$"),
      withMathJax("$$P(\\theta_d=\\pi_t-\\pi_c > \\delta | D_0, D) > c_s $$"),
      withMathJax("$$\\text{futility criterion: }$$"),
      withMathJax("$$P(\\theta_d=\\pi_t-\\pi_c > \\delta | D_0, D) < c_f $$")
    ),
    conditionalPanel(
      condition = "input.borrow == 'control arm' && input.teststat == 'logOR'",
      withMathJax("$$\\text{success criterion: }$$"),
      withMathJax("$$P(\\theta_d = \\log \\frac{\\pi_t/(1-\\pi_t)}{\\pi_c/(1-\\pi_c)} > \\delta | D_0, D) > c_s $$"),
      withMathJax("$$\\text{futility criterion: }$$"),
      withMathJax("$$P(\\theta_d = \\log \\frac{\\pi_t/(1-\\pi_t)}{\\pi_c/(1-\\pi_c)} > \\delta | D_0, D) < c_f $$")
    )
  ),
  conditionalPanel(
    condition = "input.direct == 'less than' && input.outcome =='binary'",
    conditionalPanel(
      condition = "input.borrow == 'treatment difference between arms'",
      withMathJax("$$\\text{success criterion: }$$"),
      withMathJax("$$P(\\theta_d=\\pi_t-\\pi_c<\\delta | D_0, D) > c_s $$"),
      withMathJax("$$\\text{futility criterion: }$$"),
      withMathJax("$$P(\\theta_d=\\pi_t-\\pi_c<\\delta | D_0, D) < c_f $$")
    ),
    conditionalPanel(
      condition = "input.borrow == 'log OR'",
      withMathJax("$$\\text{success criterion: }$$"),
      withMathJax("$$P(\\theta_d = \\log \\frac{\\pi_t/(1-\\pi_t)}{\\pi_c/(1-\\pi_c)}<\\delta | D_0, D) > c_s $$"),
      withMathJax("$$\\text{futility criterion: }$$"),
      withMathJax("$$P(\\theta_d = \\log \\frac{\\pi_t/(1-\\pi_t)}{\\pi_c/(1-\\pi_c)}<\\delta | D_0, D) < c_f $$")
    ),
    conditionalPanel(
      condition = "input.borrow == 'control arm' && input.teststat == 'diff'",
      withMathJax("$$\\text{success criterion: }$$"),
      withMathJax("$$P(\\theta_d=\\pi_t-\\pi_c<\\delta | D_0, D) > c_s $$"),
      withMathJax("$$\\text{futility criterion: }$$"),
      withMathJax("$$P(\\theta_d=\\pi_t-\\pi_c<\\delta | D_0, D) < c_f $$")
    ),
    conditionalPanel(
      condition = "input.borrow == 'control arm' && input.teststat == 'logOR'",
      withMathJax("$$\\text{success criterion: }$$"),
      withMathJax("$$P(\\theta_d = \\log \\frac{\\pi_t/(1-\\pi_t)}{\\pi_c/(1-\\pi_c)}<\\delta | D_0, D) > c_s $$"),
      withMathJax("$$\\text{futility criterion: }$$"),
      withMathJax("$$P(\\theta_d = \\log \\frac{\\pi_t/(1-\\pi_t)}{\\pi_c/(1-\\pi_c)}<\\delta | D_0, D) < c_f $$")
    )
  ),
  conditionalPanel(
    condition = "input.direct == 'greater than' && input.outcome =='normal'",
    withMathJax("$$\\text{success criterion: }$$"),
    withMathJax("$$P(\\theta_d = \\theta_t-\\theta_c>\\delta| D_0, D) > c_s $$"),
    withMathJax("$$\\text{futility criterion: }$$"),
    withMathJax("$$P(\\theta_d = \\theta_t-\\theta_c>\\delta| D_0, D) < c_f $$")
  ),
  conditionalPanel(
    condition = "input.direct == 'less than' && input.outcome =='normal'",
    withMathJax("$$\\text{success criterion: }$$"),
    withMathJax("$$P(\\theta_d = \\theta_c-\\theta_t<\\delta| D_0, D) > c_s $$"),
    withMathJax("$$\\text{futility criterion: }$$"),
    withMathJax("$$P(\\theta_d = \\theta_c-\\theta_t>\\delta| D_0, D) < c_f $$")
  ),
  conditionalPanel(
    condition = "input.outcome=='normal'",
    withMathJax("$$ \\text{simulate trial data D:} $$"),
    withMathJax("$$ y_c \\sim \\text{Normal}(\\mu_c, \\sigma_c) $$"),
    withMathJax("$$ y_t \\sim \\text{Normal}(\\mu_t, \\sigma_t) $$"),
    conditionalPanel(
      condition = "input.borrow=='treatment difference between arms'",
      withMathJax("$$ d \\sim \\text{Normal}(\\hat{\\mu_t} - \\hat{\\mu_c}, \\sqrt{\\hat{\\sigma_t}^2/N_t + \\hat{\\sigma_c}^2/N_c}) $$")
    )
  ),
  conditionalPanel(
    condition = "input.outcome=='binary'",
    withMathJax("$$ \\text{simulate trial data D:} $$"),
    withMathJax("$$ y_c \\sim \\text{Binomial}(N_c, \\pi_c) $$"),
    withMathJax("$$ y_t \\sim \\text{Binomial}(N_t, \\pi_t) $$"),
    conditionalPanel(
      condition = "input.borrow=='treatment difference between arms'",
      withMathJax("$$ d \\sim \\text{Normal}(\\hat{\\pi_t} - \\hat{\\pi_c}, \\sqrt{\\hat{\\pi_t} (1-\\hat{\\pi_t})/N_t + \\hat{\\pi_c} (1-\\hat{\\pi_c})/N_c}) $$")
    )
  ),
  parameter_tabs_bin,
  parameter_tabs_norm,
  numericInput('target', withMathJax("$$\\text{ decision threshold } \\delta $$"), value=0),
  #sliderInput('success', withMathJax("$$\\text{success criterion } c_s $$"), value=0.9, min=0, max=1),
  #sliderInput('futility', withMathJax("$$\\text{futility criterion } c_f $$"), value=0.5, min=0, max=1)
  textInput('success', 'success criteria per stage (comma delimited)', value = '0.9, 0.9, 0.9'),
  textInput('futility', 'futility criteria per stage (comma delimited)', value = '0, 0, 0'),
)


side_panel3 =   sidebarPanel(
  numericInput('ntrials','number of simulated trials', value=2000, min=1),
  strong('control sample size (stage 1)'),
  splitLayout(
    numericInput('ssc1', 'from', value = 30, min=2),
    numericInput('ssc2', 'to', value = 60, min=2),
    numericInput('SSi', 'interval', value=10, min=1)
  ),

  textInput('rt1', 'allocation ratio per arm (comma delimited)', value = '1, 1'),
  textInput('rt2', 'cumulative allocation ratio per stage (comma delimited)', value = '0.5, 0.75, 1'),
  actionButton('start', 'Generate Operating Characteristics', class='btn-success')
)

input_tab = tabPanel(
  title = 'settings',
  side_panel1,
  side_panel2,
  side_panel3
)


output_tab = tabPanel(
  title = 'results',
  add_busy_bar(color='red'),
  plotOutput('OC_plot', width = "60%"),
  tableOutput('OC_tab'),
  plotOutput('OC_plot_h0', width = '60%'),
  tableOutput('type1'),
  plotOutput('bayes1', width = '60%')
)

ui <- navbarPage('Bayesian Group Sequential Design with Historical Information Borrowing',
                 input_tab,
                 output_tab
)
