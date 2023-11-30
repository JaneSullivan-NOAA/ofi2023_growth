require(readxl)
require(wham)


# -------------------------------------------------------------------------
# Case 1 ------------------------------------------------------------------

#  Read data:
sheet1 = readxl::read_xlsx(path = 'case1.xlsx', sheet = 1)
sheet2 = readxl::read_xlsx(path = 'case1.xlsx', sheet = 2)
sheet3 = readxl::read_xlsx(path = 'case1.xlsx', sheet = 3)
sheet4 = readxl::read_xlsx(path = 'case1.xlsx', sheet = 4)
sheet5 = readxl::read_xlsx(path = 'case1.xlsx', sheet = 5)
sheet6 = readxl::read_xlsx(path = 'case1.xlsx', sheet = 6)
sheet7 = readxl::read_xlsx(path = 'case1.xlsx', sheet = 7)

# Make input data:
input_data = list()
input_data$ages = 1:10 # edades
input_data$years = 1976:2020 # years
input_data$n_fleets = 1 # numero de pesquerias
input_data$n_indices = 1 # numero de indices
n_years = length(input_data$years)
n_ages = length(input_data$ages)
# Agg catch:
input_data$agg_catch = matrix(sheet1$catch, ncol = input_data$n_fleets, nrow = n_years) # Obs
input_data$catch_cv = matrix(sheet1$cv, ncol = input_data$n_fleets, nrow = n_years) # Obs error
# Age comps (fishery):
input_data$catch_paa = array(as.matrix(sheet3[,3:12]), dim = c(input_data$n_fleets, n_years, n_ages)) # Obs
input_data$catch_Neff = matrix(sheet3$Nsamp, ncol = input_data$n_fleets, nrow = n_years) # Obs error
# Agg index:
input_data$agg_indices = matrix(sheet2$index, ncol = input_data$n_indices, nrow = n_years) # Obs
input_data$index_cv = matrix(sheet2$cv, ncol = input_data$n_indices, nrow = n_years) # Obs error
# Additional info:
input_data$units_indices = matrix(1L, nrow = n_years, ncol = input_data$n_indices) # 0 = numbers, 1 = biomass
input_data$fracyr_indices = matrix(sheet2$fr_yr, ncol = input_data$n_indices, nrow = n_years) # fraccion del year
# Age comps (index):
input_data$index_paa = array(as.matrix(sheet4[,3:12]), dim = c(input_data$n_indices, n_years, n_ages)) # Obs
input_data$index_Neff = matrix(sheet4$Nsamp, ncol = input_data$n_indices, nrow = n_years) # Obs error
# Selex pointers:
input_data$selblock_pointer_fleets = matrix(1L, ncol = input_data$n_fleets, nrow = n_years)
input_data$selblock_pointer_indices = matrix(2L, ncol = input_data$n_indices, nrow = n_years)
# weight-at-age information:
input_data$waa = array(0, dim = c(2, n_years, n_ages))
input_data$waa[1,,] = as.matrix(sheet6[,2:11]) # jan1
input_data$waa[2,,] = as.matrix(sheet7[,2:11]) # jul1
input_data$waa_pointer_fleets = 2
input_data$waa_pointer_indices = 1
input_data$waa_pointer_totcatch = 2
input_data$waa_pointer_ssb = 1
input_data$waa_pointer_jan1 = 1
# More information:
input_data$maturity = as.matrix(sheet5[,2:11]) # madurez sexual
input_data$fracyr_SSB = matrix(0, ncol = 1, nrow = n_years) # fraction of year SSB
input_data$Fbar_ages = 1:10 # edades para calcular Fbar
input_data$bias_correct_process = 1 # do process bias correction, 0 = no, 1 = yes
input_data$bias_correct_observation = 1 # do obs bias correction, 0 = no, 1 = yes 

# ----------
# Make input WHAM object (empirical approach)
my_input1a = prepare_wham_input(model_name = 'Case_1a',
                                basic_info = input_data,
                                NAA_re = list(N1_model = 1, N1_pars = c(1e+05, 0),
                                              recruit_model = 2, recruit_pars = 1e+05,
                                              sigma = 'rec', cor = 'iid'),
                                M = list(model = 'constant', initial_means = 0.35),
                                selectivity = list(model = c('double-normal', 'logistic'),
                                                   initial_pars = list(c(4,-2,0,0,-5,-3), c(1.5,0.3)),
                                                   n_selblocks = 2),
                                catchability = list(initial_q = 1))


# Fix some parameters:
my_input1a$map$logit_q = factor(NA)
my_input1a$map$log_N1_pars = factor(c(1,NA))

# Run model:
my_model1a = fit_wham(input = my_input1a, do.retro = FALSE, do.osa = FALSE)

# Explore outputs:
my_model1a$opt 
my_model1a$sdrep
my_model1a$rep

# Make figures:
wham::plot_wham_output(mod = my_model1a, out.type = 'pdf')

# ----------
# Make input WHAM object (WAA nonparametric)
# Add waa_cv and activate their use:
input_data$waa_cv = array(0.1, dim = dim(input_data$waa))
input_data$use_catch_waa = matrix(1L, ncol = input_data$n_fleets, nrow = n_years)
input_data$use_index_waa = matrix(1L, ncol = input_data$n_indices, nrow = n_years)

my_input1b = prepare_wham_input(model_name = 'Case_1b',
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
my_model1b = fit_wham(input = my_input1b, do.retro = FALSE, do.osa = FALSE)

# Make plots:
plot_wham_output(mod = model1b, out.type = 'pdf')

# Project:
proj1b = wham::project_wham(model = model1b, proj.opts = list(n.yrs = 3))