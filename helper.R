show_selex <- function(model = c('double-normal', 'logistic', 'len-double-normal', 'len-logistic'),
                       initial_pars = c(4,-2,0,0,-5,-3),
                       ages = 1:10, lenbins = seq(from = 2, to = 130, by = 2)) {
  require(ggplot2)
  
  # model = 'double-normal'; initial_pars = c(4,-2,0,0,-5,-3); ages = 1:10
  # model = 'logistic'; initial_pars = c(1.5,0.3); ages = 1:10
  # model = 'len-double-normal'; initial_pars = c(50,-1,4,4,-5,-2); lenbins = seq(from = 2, to = 130, by = 2)
  if(!model %in% c('double-normal', 'logistic', 'len-double-normal', 'len-logistic')) stop("this function only works for the age and length-based double-normal and logistic slx opts")
  if(model %in% c('double-normal', 'len-double-normal') & length(initial_pars) != 6) stop("initial parameter list for double-normal or len-double-normal should be a vector with length of 6")
  if(model %in% c('logistic', 'len-logistic') & length(initial_pars) != 2) stop("initial parameter list for logistic should be a vector with length of 2")
  
  n_ages = length(ages)
  amin = min(ages)
  tmp = vector(length = n_ages)
  
  n_lengths = length(lenbins)
  lmin = min(lenbins)
  ltmp = vector(length = n_lengths)
  
  if(model == 'logistic') {
    a50 = initial_pars[1]
    k = initial_pars[2]
    tmp = 1/(1+exp(-(ages-a50)/k))
  }
  
  if(model == 'len-logistic') {
    l50 = initial_pars[1]
    k = initial_pars[2]
    ltmp = 1/(1+exp(-(lenbins-l50)/k))
  }
  
  if(model == 'double-normal') {
    p_1 = initial_pars[1]
    p_2 = initial_pars[2]
    p_3 = initial_pars[3]
    p_4 = initial_pars[4]
    p_5 = 1/(1+exp(-initial_pars[5]))
    p_6 = 1/(1+exp(-initial_pars[6]))
    binwidth = 1
    gammax = p_1 + binwidth + (0.99*n_ages - p_1 - binwidth)/(1 + exp(-p_2))
    for(age in 1:n_ages) {
      alpha = p_5 + (1 - p_5)*(exp(-(age - p_1)^2/exp(p_3)) - exp(-(amin - p_1)^2/exp(p_3)))/(1-exp(-(amin - p_1)^2/exp(p_3)));
      beta = 1 + (p_6 - 1)*(exp(-(age - gammax)^2/exp(p_4)) - 1)/(exp(-(n_ages - gammax)^2/exp(p_4)) - 1);
      j_1 = 1/(1 + exp(-20*(age - p_1)/(1  + abs(age - p_1))))
      j_2 = 1/(1 + exp(-20*(age - gammax)/(1  + abs(age - gammax))))
      tmp[age] = alpha * (1 - j_1) + j_1*((1 - j_2) + j_2*beta)
    }
  }
  
  if(model == 'len-double-normal') {
    p_1 = initial_pars[1]
    p_2 = initial_pars[2]
    p_3 = initial_pars[3]
    p_4 = initial_pars[4]
    p_5 = 1/(1+exp(-initial_pars[5]))
    p_6 = 1/(1+exp(-initial_pars[6]))
    binwidth = 1
    gammax = p_1 + binwidth + (0.99*n_lengths - p_1 - binwidth)/(1 + exp(-p_2))
    for(lenbin in 1:n_lengths) {
      alpha = p_5 + (1 - p_5)*(exp(-(lenbin - p_1)^2/exp(p_3)) - exp(-(lmin - p_1)^2/exp(p_3)))/(1-exp(-(lmin - p_1)^2/exp(p_3)));
      beta = 1 + (p_6 - 1)*(exp(-(lenbin - gammax)^2/exp(p_4)) - 1)/(exp(-(n_lengths - gammax)^2/exp(p_4)) - 1);
      j_1 = 1/(1 + exp(-20*(lenbin - p_1)/(1  + abs(lenbin - p_1))))
      j_2 = 1/(1 + exp(-20*(lenbin - gammax)/(1  + abs(lenbin - gammax))))
      ltmp[lenbin] = alpha * (1 - j_1) + j_1*((1 - j_2) + j_2*beta)
    }
  }
  
  graphics.off()
  if(grepl('len', model)) {
    data.frame(length = lenbins, selex = ltmp) %>% 
      ggplot(aes(x = length, y = selex)) +
      geom_point() +
      geom_line()
  } else {
  data.frame(age = ages, selex = tmp) %>% 
    ggplot(aes(x = age, y = selex)) +
    geom_point() +
    geom_line()
    }
}