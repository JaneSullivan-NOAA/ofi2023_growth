require(readxl)
require(wham)
require(dplyr)

# Case 3 ------------------------------------------------------------------
sheet1 = readxl::read_xlsx(path = 'case3.xlsx', sheet = 1)
sheet2 = readxl::read_xlsx(path = 'case3.xlsx', sheet = 2)
sheet3 = readxl::read_xlsx(path = 'case3.xlsx', sheet = 3)
sheet4 = readxl::read_xlsx(path = 'case3.xlsx', sheet = 4)
sheet5 = readxl::read_xlsx(path = 'case3.xlsx', sheet = 5)
sheet6 = readxl::read_xlsx(path = 'case3.xlsx', sheet = 6)

# Make input data:
input_data = list()
input_data$ages = 1:10 # edades
input_data$years = 1976:2020 # years
input_data$lengths = seq(from = 2, to = 130, by = 2) # tallas
input_data$n_fleets = 1 # numero de pesquerias
input_data$n_indices = 1 # numero de indices
n_years = length(input_data$years)
n_ages = length(input_data$ages)
n_lengths = length(input_data$lengths)
# Agg catch:
input_data$agg_catch = matrix(sheet1$catch, ncol = input_data$n_fleets, nrow = n_years) # Obs
input_data$catch_cv = matrix(sheet1$cv, ncol = input_data$n_fleets, nrow = n_years) # Obs error
# Length comps (fishery)
input_data$catch_pal = array(as.matrix(sheet3[,3:67]), dim = c(input_data$n_fleets, n_years, n_lengths)) # Obs
input_data$catch_NeffL = matrix(sheet3$Nsamp, ncol = input_data$n_fleets, nrow = n_years) # Obs error
input_data$use_catch_pal = matrix(1, nrow = n_years, ncol = input_data$n_fleets) # 1 = usar, 0 = no usar
# Agg index:
input_data$agg_indices = matrix(sheet2$index, ncol = input_data$n_indices, nrow = n_years) # Obs
input_data$index_cv = matrix(sheet2$cv, ncol = input_data$n_indices, nrow = n_years) # Obs error
# Additional information:
input_data$units_indices = matrix(1L, nrow = n_years, ncol = input_data$n_indices) # 0 = numbers, 1 = biomass
input_data$fracyr_indices = matrix(sheet2$fr_yr, ncol = input_data$n_indices, nrow = n_years) # fraccion del year
# Length comps (index):
input_data$index_pal = array(as.matrix(sheet4[,3:67]), dim = c(input_data$n_indices, n_years, n_lengths)) # Obs
input_data$index_NeffL = matrix(sheet4$Nsamp, ncol = input_data$n_indices, nrow = n_years) # Obs error
input_data$use_index_pal = matrix(1, nrow = n_years, ncol = input_data$n_indices) # 1 = usar, 0 = no usar
# CAAL index:
input_data$index_caal = array(0, dim = c(input_data$n_indices, n_years, n_lengths, n_ages))
for(i in seq_along(input_data$years)){
  tmp = sheet5 %>% filter(year == input_data$years[i])
  input_data$index_caal[1,i,,] = as.matrix(tmp[,4:13])
}
# CAAL Neff index:
input_data$index_caal_Neff = array(0, dim = c(n_years, input_data$n_indices, n_lengths))
for(i in seq_along(input_data$years)){
  tmp = sheet5 %>% filter(year == input_data$years[i])
  input_data$index_caal_Neff[i,1,] = as.matrix(tmp[,3])
}
# CAAL index use/not use (use only when Nsamp > 0):
input_data$use_index_caal = array(0, dim = c(n_years, input_data$n_indices, n_lengths))
input_data$use_index_caal[input_data$index_caal_Neff > 0] = 1
# Dont forget to turn off the use of paa (default):
input_data$use_catch_paa = matrix(0, nrow = n_years, ncol = input_data$n_fleets) # 1 = use, 0 = no use
input_data$use_index_paa = matrix(0, nrow = n_years, ncol = input_data$n_indices) # 1 = use, 0 = no use
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
input_data$maturity = as.matrix(sheet6[,2:11]) # maturity
input_data$fracyr_SSB = matrix(0, ncol = 1, nrow = n_years) # fraction of year SSB
input_data$Fbar_ages = 1:10 # edades para calcular Fbar
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
                                                     n_selblocks = 2), # Selectivity parameter
                                  catchability = list(initial_q = 1), # Catchability parameter
                                  growth = list(model = 'vB_classic', 
                                                init_vals = c(0.2, 90, 10), est_pars = 1:3, 
                                                SD_vals = c(1, 9)),
                                  LW = list(init_vals = c(5.56e-06, 3.2))
) 

# Fix some parameters:
input3$map$log_N1_pars = factor(c(1, NA))
input3$map$logit_q = factor(NA)
# Optional: fix sigmaR and turn off NAA (recruitment) as random variable (to make it faster):
input3$par$log_NAA_sigma = log(0.6)
input3$map$log_NAA_sigma = factor(NA)
input3$random = NULL

# Run model:
model3 = fit_wham(input = input3, do.retro = FALSE, do.osa = FALSE)

# Make plots:
plot_wham_output(model3, out.type = 'pdf')