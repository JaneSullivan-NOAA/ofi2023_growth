show_selex <- function(model = c('double-normal', 'logistic'),
                       initial_pars = c(4,-2,0,0,-5,-3),
                       ages = 1:10) {
  require(ggplot2)
  
  # model = 'double-normal'; initial_pars = c(4,-2,0,0,-5,-3); ages = 1:10
  # model = 'logistic'; initial_pars = c(1.5,0.3); ages = 1:10
  if(!model %in% c('double-normal', 'logistic')) stop("this function only works for the double-normal and logistic slx opts")
  if(model == 'double-normal' & length(initial_pars) != 6) stop("initial parameter list for double-normal should be a vector with length of 6")
  if(model == 'logistic' & length(initial_pars) != 2) stop("initial parameter list for logistic should be a vector with length of 2")
  
  amin = min(ages)
  tmp = vector(length = n_ages)
  
  if(model == 'logistic') {
    a50 = initial_pars[1]
    k = initial_pars[2]
    tmp = 1/(1+exp(-(ages-a50)/k))
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
  graphics.off()
  data.frame(age = ages, selex = tmp) %>% 
    ggplot(aes(x = age, y = selex)) +
    geom_point() +
    geom_line()
}