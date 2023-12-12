# Case 1: compare empirical WAA and nonparametric WAA 

# Model info: one fishery that occurs in July, one survey that occurs in
# January, age comps for both fishery and survey

# remotes::install_github(repo = 'GiancarloMCorrea/wham', ref='growth', INSTALL_opts = c("--no-docs", "--no-multiarch", "--no-demo"))
library(readxl)
library(wham)
library(ggplot2)
library(dplyr)
source("helper.R") # not needed, includes function to help with selectivity inputs
dir.create("case1a_empiricalWAA")
dir.create("case1b_nonparametricWAA")
runmodels = FALSE # use RDS model object files instead of running them

#  Read data ----
catch_df = readxl::read_xlsx(path = 'case1.xlsx', sheet = 1)
index_df = readxl::read_xlsx(path = 'case1.xlsx', sheet = 2)
fsh_comp_df = readxl::read_xlsx(path = 'case1.xlsx', sheet = 3)
index_comp_df = readxl::read_xlsx(path = 'case1.xlsx', sheet = 4)
maturity_df = readxl::read_xlsx(path = 'case1.xlsx', sheet = 5)
waa_jan1_df = readxl::read_xlsx(path = 'case1.xlsx', sheet = 6)
waa_jul1_df = readxl::read_xlsx(path = 'case1.xlsx', sheet = 7)

# Make input data:
input_data = list()
input_data$ages = 1:10 
input_data$years = 1976:2020 
input_data$n_fleets = 1  
input_data$n_indices = 1 
n_years = length(input_data$years)
n_ages = length(input_data$ages)
# Agg catch:
input_data$agg_catch = matrix(catch_df$catch, ncol = input_data$n_fleets, nrow = n_years) # Obs
input_data$catch_cv = matrix(catch_df$cv, ncol = input_data$n_fleets, nrow = n_years) # Obs error
# Age comps (fishery):
input_data$catch_paa = array(as.matrix(fsh_comp_df[,3:12]), dim = c(input_data$n_fleets, n_years, n_ages)) # Obs
input_data$catch_Neff = matrix(fsh_comp_df$Nsamp, ncol = input_data$n_fleets, nrow = n_years) # Obs error
# Agg index:
input_data$agg_indices = matrix(index_df$index, ncol = input_data$n_indices, nrow = n_years) # Obs
input_data$index_cv = matrix(index_df$cv, ncol = input_data$n_indices, nrow = n_years) # Obs error
# Additional info:
input_data$units_indices = matrix(1L, nrow = n_years, ncol = input_data$n_indices) # 0 = numbers, 1 = biomass
input_data$fracyr_indices = matrix(index_df$fr_yr, ncol = input_data$n_indices, nrow = n_years) # fraction of the year when survey occurs
# Age comps (index):
input_data$index_paa = array(as.matrix(index_comp_df[,3:12]), dim = c(input_data$n_indices, n_years, n_ages)) # Obs
input_data$index_Neff = matrix(index_comp_df$Nsamp, ncol = input_data$n_indices, nrow = n_years) # Obs error
# Selex pointers:
input_data$selblock_pointer_fleets = matrix(1L, ncol = input_data$n_fleets, nrow = n_years)
input_data$selblock_pointer_indices = matrix(2L, ncol = input_data$n_indices, nrow = n_years)
# weight-at-age information:
input_data$waa = array(0, dim = c(2, n_years, n_ages))
input_data$waa[1,,] = as.matrix(waa_jan1_df[, 2:11]) # jan1 (survey, spawning)
input_data$waa[2,,] = as.matrix(waa_jul1_df[, 2:11]) # jul1 (fishery)
input_data$waa_pointer_fleets = 2
input_data$waa_pointer_indices = 1
input_data$waa_pointer_totcatch = 2
input_data$waa_pointer_ssb = 1
input_data$waa_pointer_jan1 = 1
# More information:
input_data$maturity = as.matrix(maturity_df[,2:11]) # maturity
input_data$fracyr_SSB = matrix(0, ncol = 1, nrow = n_years) # spawning fraction (0 = spawn at beginning of year)
input_data$Fbar_ages = 1:10 # ages to include in mean F calculation
input_data$bias_correct_process = 1 # do process bias correction, 0 = no, 1 = yes
input_data$bias_correct_observation = 1 # do obs bias correction, 0 = no, 1 = yes 

# Empirical WAA model ----

# Make input WHAM object (empirical WAA approach)
my_input1a = prepare_wham_input(model_name = 'Case1_empiricalWAA',
                                basic_info = input_data,
                                NAA_re = list(N1_model = 1, # estimating a mean recruitment with yearly recruitment as random effects
                                              N1_pars = c(1e+05, 0), # initial numbers in the first age class, and equilib F rate generating the rest of the NAA in the first year
                                              recruit_model = 2, # random about mean
                                              recruit_pars = 1e+05, # mean rec
                                              sigma = 'rec', # rand eff on rec devs, all other ages deterministic
                                              cor = 'iid'), 
                                M = list(model = 'constant', initial_means = 0.35),
                                selectivity = list(model = c('double-normal', 'logistic'),
                                                   initial_pars = list(c(4,-2,0,0,-5,-3), c(1.5,0.3)),
                                                   n_selblocks = 2),
                                catchability = list(initial_q = 1))

show_selex(model = "double-normal", initial_pars = c(4,-2,0,0,-5,-3), ages = input_data$ages)
show_selex(model = "logistic", initial_pars = c(1.5,0.3), ages = input_data$ages)

# Fix some parameters:
my_input1a$map$logit_q = factor(NA)
my_input1a$map$log_N1_pars = factor(c(1,NA))

# Run model:
if(isTRUE(runmodels)) {
  my_model1a = wham::fit_wham(MakeADFun.silent = TRUE, input = my_input1a, do.retro = FALSE, do.osa = FALSE)
  proj1a = wham::project_wham(model = my_model1a, proj.opts = list(n.yrs = 3))
  saveRDS(my_model1a, file = "case1a_empiricalWAA/my_model1a.RDS")
  saveRDS(proj1a, file = "case1a_empiricalWAA/proj1a.RDS")
  
} else {
  my_model1a <- readRDS(file = "case1a_empiricalWAA/my_model1a.RDS")
  proj1a <- readRDS(file = "case1a_empiricalWAA/proj1a.RDS")
}

# Explore outputs:
my_model1a$opt 
my_model1a$sdrep
my_model1a$rep
names(my_model1a)

# Make figures:
wham::plot_wham_output(mod = my_model1a, dir.main = "case1a_empiricalWAA", out.type = 'pdf')

# Nonparametric WAA model ----

# Add waa_cv and activate their use:
input_data$waa_cv = array(0.1, dim = dim(input_data$waa))
input_data$use_catch_waa = matrix(1L, ncol = input_data$n_fleets, nrow = n_years)
input_data$use_index_waa = matrix(1L, ncol = input_data$n_indices, nrow = n_years)

my_input1b = wham::prepare_wham_input(model_name = 'Case_1b',
                                basic_info = input_data,
                                NAA_re = list(N1_model = 1, N1_pars = c(1e+05, 0),
                                              recruit_model = 2, recruit_pars = 1e+05,
                                              sigma = 'rec', cor = 'iid'),
                                M = list(model = 'constant', initial_means = 0.35),
                                selectivity = list(model = c('double-normal', 'logistic'),
                                                   initial_pars = list(c(4,-2,0,0,-5,-3), c(1.5,0.3)),
                                                   n_selblocks = 2),
                                catchability = list(initial_q = 1),
                                WAA = list(re = '2dar1', WAA_vals = colMeans(input_data$waa[1,,])))


# Fix some parameters:
my_input1b$map$logit_q = factor(NA)
my_input1b$map$log_N1_pars = factor(c(1,NA))

# Run models:
if(isTRUE(runmodels)) {
  my_model1b = wham::fit_wham(MakeADFun.silent = TRUE, input = my_input1b, do.retro = FALSE, do.osa = FALSE)
  proj1b = wham::project_wham(model = my_model1b, proj.opts = list(n.yrs = 3))
  saveRDS(my_model1b, file = "case1b_nonparametricWAA/my_model1b.RDS")
  saveRDS(proj1b, file = "case1b_nonparametricWAA/proj1b.RDS")
} else {
  my_model1b <- readRDS(file = "case1b_nonparametricWAA/my_model1b.RDS")
  proj1b <- readRDS(file = "case1b_nonparametricWAA/proj1b.RDS")
}

# Make plots:
wham::plot_wham_output(mod = my_model1b, dir.main = "case1b_nonparametricWAA", out.type = 'pdf')

# Compare projections ----

x = summary(proj1a$sdrep)
x[rownames(x) == "pred_waa",]

unique(rownames(x)) # list of estimated parameters and derived quantities with SE
x = x[rownames(x) == "log_SSB",] # SSB estimates with SE
tmp = tibble::as_tibble(cbind(proj1a$years_full, exp(cbind(x, x[,1] + qnorm(0.975)*cbind(-x[,2],x[,2])))/1000)) # calculate 95% CI
colnames(tmp) <- c("year", "SSB","SSB_se","SSB_lower","SSB_upper")
ssb <- tmp %>% dplyr::mutate(model = "empirical_waa")

x = summary(proj1b$sdrep)
pred_waa <- x[rownames(x) == "pred_waa",]
unique(rownames(x)) # list of estimated parameters and derived quantities with SE
x = x[rownames(x) == "log_SSB",] # SSB estimates with SE
tmp = tibble::as_tibble(cbind(proj1b$years_full, exp(cbind(x, x[,1] + qnorm(0.975)*cbind(-x[,2],x[,2])))/1000)) # calculate 95% CI
colnames(tmp) <- c("year", "SSB","SSB_se","SSB_lower","SSB_upper")
ssb <- ssb %>% 
  bind_rows(tmp %>% dplyr::mutate(model = "nonparametric_waa")) 

ggplot(ssb, aes(x = year, y = SSB, col = model, fill = model, shape = model)) +
  geom_ribbon(aes(ymin = SSB_lower, ymax = SSB_upper), alpha = 0.3, col = NA) +
  geom_line() +
  geom_point() +
  geom_vline(xintercept = max(proj1b$years), lty = 2) +
  expand_limits(y = 0) +
  theme_bw(15)
ggsave("case1_ssb.png", device = "png")

# Compare WAA -----
tmp <- as.data.frame(cbind(proj1a$years, my_input1a$data$waa[1,,]))
colnames(tmp) <- c("year", 1:10)
waa <- tmp %>% 
  dplyr::mutate(model = "empirical_waa") %>% 
  tidyr::pivot_longer(cols = 2:11, names_to = 'age', values_to = "weight")

# Default empirical WAA approach uses the most recent 5 yr mean WAA for
# projections
waa <- dplyr::bind_rows(
  waa,
  tidyr::expand_grid(
    year = max(waa$year)+1:3,
    model = "empirical_waa",
    waa %>% 
      dplyr::filter(year > max(year)-5) %>% 
      dplyr::group_by(age) %>% 
      dplyr::summarise(weight = mean(weight)))
)
          
tmp <- as.data.frame(cbind(proj1a$years_full, array(unname(pred_waa[,1]), dim = c(2,length(proj1b$years_full),10))[1,,])) 
colnames(tmp) <- c("year", 1:10)
waa <- waa %>% 
  dplyr::bind_rows(
    tmp %>% 
      dplyr::mutate(model = "nonparametric_waa") %>% 
      tidyr::pivot_longer(cols = 2:11, names_to = 'age', values_to = "weight")
  ) %>% 
  dplyr::mutate(age = factor(age, levels = as.character(1:10)))

ggplot(waa, aes(x = year, y = weight, col = model, shape = model)) +
  geom_point() +
  geom_line() +
  facet_wrap(~age, scales = "free_y", dir = "v", ncol = 2) +
  theme_bw(15) +
  theme(legend.position = "top")
ggsave("case1_waa.png", device = "png")



              