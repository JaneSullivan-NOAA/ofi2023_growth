# Case 3: fit to conditional age-at-length (CAAL) data 

# remotes::install_github(repo = 'GiancarloMCorrea/wham', ref='growth', INSTALL_opts = c("--no-docs", "--no-multiarch", "--no-demo"))
library(readxl)
library(wham)
library(ggplot2)
library(dplyr)
source("helper.R")
dir.create("case3_CAAL")
runmodels = FALSE

# Case 3 ------------------------------------------------------------------
catch_df = readxl::read_xlsx(path = 'case3.xlsx', sheet = 1)
index_df = readxl::read_xlsx(path = 'case3.xlsx', sheet = 2)
fsh_lcomp_df = readxl::read_xlsx(path = 'case3.xlsx', sheet = 3)
index_lcomp_df = readxl::read_xlsx(path = 'case3.xlsx', sheet = 4)
caal_df = readxl::read_xlsx(path = 'case3.xlsx', sheet = 5)
maturity_df = readxl::read_xlsx(path = 'case3.xlsx', sheet = 6)

# Make input data:
input_data = list()
input_data$ages = 1:10 # ages
input_data$years = 1976:2020 # years
input_data$lengths = seq(from = 2, to = 130, by = 2) # lengths
input_data$n_fleets = 1 # number of fleets
input_data$n_indices = 1 # number of indices
n_years = length(input_data$years)
n_ages = length(input_data$ages)
n_lengths = length(input_data$lengths)
# Agg catch:
input_data$agg_catch = matrix(catch_df$catch, ncol = input_data$n_fleets, nrow = n_years) # Obs
input_data$catch_cv = matrix(catch_df$cv, ncol = input_data$n_fleets, nrow = n_years) # Obs error
# Length comps (fishery)
input_data$catch_pal = array(as.matrix(fsh_lcomp_df[,3:67]), dim = c(input_data$n_fleets, n_years, n_lengths)) # Obs
input_data$catch_NeffL = matrix(fsh_lcomp_df$Nsamp, ncol = input_data$n_fleets, nrow = n_years) # Obs error
input_data$use_catch_pal = matrix(1, nrow = n_years, ncol = input_data$n_fleets) # 1 = fit, 0 = don't fit
# Agg index:
input_data$agg_indices = matrix(index_df$index, ncol = input_data$n_indices, nrow = n_years) # Obs
input_data$index_cv = matrix(index_df$cv, ncol = input_data$n_indices, nrow = n_years) # Obs error
# Additional information:
input_data$units_indices = matrix(1L, nrow = n_years, ncol = input_data$n_indices) # 0 = numbers, 1 = biomass
input_data$fracyr_indices = matrix(index_df$fr_yr, ncol = input_data$n_indices, nrow = n_years) # fraction of the year when survey occurs
# Length comps (index):
input_data$index_pal = array(as.matrix(index_lcomp_df[,3:67]), dim = c(input_data$n_indices, n_years, n_lengths)) # Obs
input_data$index_NeffL = matrix(index_lcomp_df$Nsamp, ncol = input_data$n_indices, nrow = n_years) # Obs error
input_data$use_index_pal = matrix(1, nrow = n_years, ncol = input_data$n_indices) # 1 = fit, 0 = don't fit
# CAAL index:
input_data$index_caal = array(0, dim = c(input_data$n_indices, n_years, n_lengths, n_ages))
for(i in seq_along(input_data$years)){
  tmp = caal_df %>% filter(year == input_data$years[i])
  input_data$index_caal[1,i,,] = as.matrix(tmp[,4:13])
}
# CAAL Neff index:
input_data$index_caal_Neff = array(0, dim = c(n_years, input_data$n_indices, n_lengths))
for(i in seq_along(input_data$years)){
  tmp = caal_df %>% filter(year == input_data$years[i])
  input_data$index_caal_Neff[i,1,] = as.matrix(tmp[,3])
}
# CAAL index use/not use (use only when Nsamp > 0):
input_data$use_index_caal = array(0, dim = c(n_years, input_data$n_indices, n_lengths))
input_data$use_index_caal[input_data$index_caal_Neff > 0] = 1
# Dont forget to turn off the use of paa (default):
input_data$use_catch_paa = matrix(0, nrow = n_years, ncol = input_data$n_fleets) # 1 = fit, 0 = don't fit
input_data$use_index_paa = matrix(0, nrow = n_years, ncol = input_data$n_indices) # 1 = fit, 0 = don't fit
# Selex pointers:
input_data$selblock_pointer_fleets = matrix(1L, ncol = input_data$n_fleets, nrow = n_years)
input_data$selblock_pointer_indices = matrix(2L, ncol = input_data$n_indices, nrow = n_years)
# weight-at-age pointers information:
input_data$waa_pointer_fleets = 2
input_data$waa_pointer_indices = 1
input_data$waa_pointer_totcatch = 2
input_data$waa_pointer_ssb = 1
input_data$waa_pointer_jan1 = 1
# More information:
input_data$maturity = as.matrix(maturity_df[,2:11]) # maturity
input_data$fracyr_SSB = matrix(0, ncol = 1, nrow = n_years) # spawning fraction (0 = spawn at beginning of year)
input_data$Fbar_ages = 1:10 # es to include in mean F calculation 
input_data$bias_correct_process = 1 # do process bias correction, 0 = no, 1 = yes
input_data$bias_correct_observation = 1 # do obs bias correction, 0 = no, 1 = yes 

# Create WHAM input: (LAA parametric approach)
input3 = wham::prepare_wham_input(model_name = "Case_3",
                                  basic_info = input_data, 
                                  NAA_re = list(N1_model = 1, N1_pars = c(1e+05, 0),
                                                recruit_model = 2, recruit_pars = 1e+05,
                                                sigma = 'rec', cor = 'iid'), # Recruitment parameters
                                  M = list(model = 'constant', initial_means = 0.35), # M parameter
                                  selectivity = list(model = c('len-double-normal', 'len-logistic'),
                                                     initial_pars = list(c(50,-1,4,4,-5,-2), c(15, 3)),
                                                     n_selblocks = 2, fix_pars = list(c(1), NULL)
                                                     ), # Selectivity parameter
                                  catchability = list(initial_q = 1), # Catchability parameter
                                  growth = list(model = 'vB_classic', 
                                                init_vals = c(0.2, 90, 10), est_pars = 1:3, 
                                                SD_vals = c(1, 9)),
                                  LW = list(init_vals = c(5.56e-06, 3.2))
) 

show_selex(model = "len-double-normal", initial_pars = c(50,-1,4,4,-5,-2))


# fix len-logistic selex
(2+130-2)/(1+exp(-(-0.08407463)))
(2+130-2)/(1+exp(-(0.35989867)))
input3$data$selpars_lower;
input3$data$selpars_upper

# Fix some parameters:
input3$map$log_N1_pars = factor(c(1, NA))
input3$map$logit_q = factor(NA)
# Optional: fix sigmaR and turn off NAA (recruitment) as random variable (to make it faster):
input3$par$log_NAA_sigma = log(0.6)
input3$map$log_NAA_sigma = factor(NA)
input3$random = NULL

# Run model:
if(isTRUE(runmodels)) {
  my_model3 = wham::fit_wham(MakeADFun.silent = TRUE, input = input3, do.retro = FALSE, do.osa = FALSE)
  saveRDS(my_model3, file = "case3_CAAL/my_model3.RDS")

} else {
  my_model3 <- readRDS(file = "case3_CAAL/my_model3.RDS")
}

my_model3$opt 
my_model3$sdrep
my_model3$rep
names(my_model3)

# Make plots:
plot_wham_output(my_model3, dir.main = "case3_CAAL", out.type = 'pdf')
