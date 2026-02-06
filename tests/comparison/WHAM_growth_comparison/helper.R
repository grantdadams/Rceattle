# https://giancarlomcorrea.netlify.app/labs/OFI_WK_2023/examples/helper.R

show_selex <- function(model = c('double-normal', 'logistic', 'len-double-normal', 'len-logistic'),
                       initial_pars = c(4,-2,0,0,-5,-3),
                       ages = 1:10, lenbins = seq(from = 2, to = 130, by = 2)) {
  require(ggplot2)

  # model = 'double-normal'; initial_pars = c(4,-2,0,0,-5,-3); ages = 1:10
  # model = 'logistic'; initial_pars = c(1.5,0.3); ages = 1:10
  # model = 'len-double-normal'; initial_pars = c(50,-1,4,4,-5,-2); lenbins = seq(from = 2, to = 130, by = 2); binwidth = 2
  if(!model %in% c('double-normal', 'logistic', 'len-double-normal', 'len-logistic')) stop("this function only works for the age and length-based double-normal and logistic slx opts")
  if(model %in% c('double-normal', 'len-double-normal') & length(initial_pars) != 6) stop("initial parameter list for double-normal or len-double-normal should be a vector with length of 6")
  if(model %in% c('logistic', 'len-logistic') & length(initial_pars) != 2) stop("initial parameter list for logistic should be a vector with length of 2")

  binwidth = ages[2]-ages[1]
  lbinwidth = lenbins[2]-lenbins[1]

  n_ages = length(ages)
  amin = min(ages)
  tmp = vector(length = n_ages)

  n_lengths = length(lenbins)
  new_lenbins = lenbins + 0.5 * lbinwidth
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
    ltmp = 1/(1+exp(-(new_lenbins-l50)/k))
  }

  if(model == 'double-normal') {
    p_1 = initial_pars[1]
    p_2 = initial_pars[2]
    p_3 = initial_pars[3]
    p_4 = initial_pars[4]
    p_5 = 1/(1+exp(-initial_pars[5]))
    p_6 = 1/(1+exp(-initial_pars[6]))
    gammax = p_1 + binwidth + (0.99*n_ages - p_1 - binwidth)/(1 + exp(-p_2))
    for(age in 1:n_ages) {
      alpha = p_5 + (1 - p_5)*(exp(-(age - p_1)^2/exp(p_3)) - exp(-(amin - p_1)^2/exp(p_3)))/(1-exp(-(amin - p_1)^2/exp(p_3)));
      beta = 1 + (p_6 - 1)*(exp(-(age - gammax)^2/exp(p_4)) - 1)/(exp(-(n_ages - gammax)^2/exp(p_4)) - 1);
      j_1 = 1/(1 + exp(-20*(age - p_1)/(1  + abs(age - p_1))))
      j_2 = 1/(1 + exp(-20*(age - gammax)/(1  + abs(age - gammax))))
      tmp[age] = alpha * (1 - j_1) + j_1*((1 - j_2) + j_2*beta)
    }
    cat("\n",
        "Initial values for age-based double normal selectivity:\n",
        "p1 = first age at which selex = 1\n",
        "p2 = controls the width of the peak (small to negative number gives narrow peak, larger number gives wider peak)\n",
        "p3 = slope of the ascending limb (small to negative number gives knife edge/steep slope, larger number gives gradual slope)\n",
        "p4 = slope of the descending limb (same as p3)\n",
        "p5 = selex at the minimum length (parameterized in logit-space, p5=-10 is selex=0, p5=0 is selex=0.5, p5=10 is selex=1)\n",
        "p6 = selex at the max age (same as p5)\n\n",
        "gamma (function of p2) is the oldest age where selex = 1, gamma =", round(gammax,2), " "
    )
  }

  if(model == 'len-double-normal') {
    p_1 = initial_pars[1]
    p_2 = initial_pars[2]
    p_3 = initial_pars[3]
    p_4 = initial_pars[4]
    p_5 = 1/(1+exp(-initial_pars[5]))
    p_6 = 1/(1+exp(-initial_pars[6]))
    lmax = max(new_lenbins)
    lmin = min(new_lenbins)
    gammax = p_1 + lbinwidth + (0.99*lmax - p_1 - lbinwidth)/(1 + exp(-p_2))
    for(l in 1:n_lengths) {
      alpha = p_5 + (1 - p_5)*(exp(-(new_lenbins[l] - p_1)^2/exp(p_3)) - exp(-(lmin - p_1)^2/exp(p_3)))/(1-exp(-(lmin - p_1)^2/exp(p_3)));
      beta = 1 + (p_6 - 1)*(exp(-(new_lenbins[l] - gammax)^2/exp(p_4)) - 1)/(exp(-(lmax - gammax)^2/exp(p_4)) - 1);
      j_1 = 1/(1 + exp(-20*(new_lenbins[l] - p_1)/(1  + abs(new_lenbins[l] - p_1))))
      j_2 = 1/(1 + exp(-20*(new_lenbins[l] - gammax)/(1  + abs(new_lenbins[l] - gammax))))
      ltmp[l] = alpha * (1 - j_1) + j_1*((1 - j_2) + j_2*beta)
    }
    cat("\n",
        "Initial values for length-based double normal selectivity:\n",
        "p1 = first length at which selex = 1\n",
        "p2 = controls the width of the peak (small to negative number gives narrow peak, larger number gives wider peak)\n",
        "p3 = slope of the ascending limb (small to negative number gives knife edge/steep slope, larger number gives gradual slope)\n",
        "p4 = slope of the descending limb (same as p3)\n",
        "p5 = selex at the minimum length (parameterized in logit-space, p5=-10 is selex=0, p5=0 is selex=0.5, p5=10 is selex=1)\n",
        "p6 = selex at the max length (same as p5)\n\n",
        "gamma (function of p2) is the largest length where selex = 1, gamma =", round(gammax,0), " "
    )
  }

  graphics.off()
  if(grepl('len', model)) {
    data.frame(length = lenbins, selex = ltmp) %>%
      ggplot(aes(x = length, y = selex)) +
      geom_point() +
      geom_line() +
      expand_limits(y = 0)
  } else {
    data.frame(age = ages, selex = tmp) %>%
      ggplot(aes(x = age, y = selex)) +
      geom_point() +
      geom_line() +
      expand_limits(y = 0)
  }
}
