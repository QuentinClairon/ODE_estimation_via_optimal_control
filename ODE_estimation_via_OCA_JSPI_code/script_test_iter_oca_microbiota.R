###### Microbiota Simulation: Well-specified #########
pathnames <- list.files(pattern="[.]R$", path="function_model_and_sims//", full.names=TRUE);
sapply(pathnames, FUN=source);
pathnames <- list.files(pattern="[.]R$", path="function_util_estimation//", full.names=TRUE);
sapply(pathnames, FUN=source);

library('deSolve')
library('optimx')

######## Data generation
# observation interval & true parameter specification

# level of noise tested in the article
sigma_obs = 0.01
#sigma_obs = 0.02
#sigma_obs = 0.03

Times_obs <- seq(0, 23, by = 23/25)
nb_obs= length(Times_obs)

dim_syst = 7
dim_control = 7

M1 = c(-0.20516,	0.098398,	0.16739,	-0.16461,	-0.14341,	0.019881,	-0.51535,	-0.39162,	0.34635,	0.0088853,	-0.26894)
M2 = c(0.062123,	-0.10489,	-0.043011,	-0.15466,	-0.1872,	0.027031,	-0.45919,	-0.41388,	0.3013,	0.022081,	-0.19657)
M3 = c(0.14373,	-0.19203,	-0.10162,	-0.13971,	-0.16537,	0.013651,	-0.50414,	-0.7724,	0.29257, -0.005959,	-0.20645)
M4 = c(0.22403,	0.13813,	0.00045883,	-0.83125,	-0.2238,	0.22027,	-0.20529,	-1.0097,	0.66639,	-0.038986,	-0.40032)
M5 = c(-0.18016,	-0.051261,	-5.03E-05,	-0.054212,	-0.70858,	0.016198,	-0.50756,	0.55363,	0.15757,	0.22438,	0.10635)
M6 = c(-0.11159,	-0.03721,	-0.042591,	0.041044,	0.26134,	-0.42266,	-0.18536,	-0.43231,	0.1647,	-0.061038,	-0.26461)
M7 = c(-0.12669,	-0.18576,	-0.12222,	0.3809,	0.4003,	-0.16078,	-1.2124,	1.3897,	-0.37922,	0.19189,	-0.096352)
M8 = c(-0.071257,	0.00060448,	0.080355,	-0.4548,	-0.50349,	0.16899,	-0.56222,	-4.3508,	0.44315,	-0.22341,	-0.2074)
M9 = c(-0.037541,	-0.033333,	-0.049912,	-0.090424,	-0.10211,	0.03229,	-0.18179,	-0.30301,	-0.055765,	0.01436,	-0.0076697)
M10 = c(-0.04225,	-0.013105,	0.02398,	-0.11784,	-0.32893,	0.020748,	0.054767,	-2.0963,	0.11124,	-0.19213,	0.023816)
M11 = c(-0.3742,	0.27843,	0.24887,	-0.16829,	0.08399,	0.033691,	-0.23242,	-0.39513,	0.31454,	-0.038764,	-0.3841)

matrix_M = rbind(M1,M2,M3,M4,M5,M6,M7,M8,M9,M10,M11)

vect_g = c(0.36807, 0.31023, 0.3561, 0.54006, 0.70898, 0.47064, 0.2297, 0.83005, 0.39181, 0.29075, 0.32367)

vect_s = c(-3.2926, -3.0354,-2.0909,-1.9395,-1.3491,-1.1018,-0.92446,-0.79401,-0.31272,1.0671,3.7009)

x0_int_mic1 = c(0.0001,0.0001,    0.0098,    0.0002,    0.0074,    0.000001,         0.000001,    0.0036,    0.000001,   0.000001,    0.0116)
x0_int_mic2 = c(0.2333,    0.4832,    1.7757,    0.0912,    0.0096,    0.0007,         0.000001,   0.0003,    0.0005,    0.0001,    0.0754)
x0_int_mic3 = c(0.0003,    0.0045,    1.3970,    0.0020,    1.0672,    0.0003 ,        0.000001,   0.2005 ,   0.0035 ,   0.0001,    0.6187)

theta = c(matrix(matrix_M,1,121),vect_g,vect_s)

ind_state_sel = c(1:5,9,11)
vect_g_sp = vect_g[ind_state_sel]
vect_s_sp = vect_s[ind_state_sel]
matrix_M_sp = rbind(cbind(matrix_M [1:5,1:5],matrix_M[1:5,9],matrix_M[1:5,11]),
                    c(matrix_M[9,1:5],matrix_M[9,9], matrix_M[9,11]),
                    c(matrix_M[11,1:5], matrix_M[11,9],matrix_M[11,11]) )  


theta_ini_inter = c(t(matrix(matrix_M_sp[1:5,],35,1)), matrix_M_sp[7,])
theta_ini_sp3 = theta_ini_inter[c(1,3:5, 7:13, 16:17, 19, 22:25, 28, 30:38, 40:42)]


# noise level speicfication and data generation for the three subjects

x0_int_mic1_sp = x0_int_mic1[ind_state_sel]
x0_int_mic2_sp = x0_int_mic2[ind_state_sel]
x0_int_mic3_sp = x0_int_mic3[ind_state_sel]
x0_log_tot = c(log(pmax(x0_int_mic1_sp,10^-10)),log(pmax(x0_int_mic2_sp,10^-10)),log(pmax(x0_int_mic3_sp,10^-10)))

theta_sp = c(matrix(matrix_M_sp,1,49),vect_g_sp,vect_s_sp)

out_log_transf_sp1 <- ode(func = Microbiota_sim_log_transf_rest, y = log(pmax(x0_int_mic1_sp,10^-10)), parms =theta_sp, times = Times_obs)
out_log_transf_sp2 <- ode(func = Microbiota_sim_log_transf_rest, y = log(pmax(x0_int_mic2_sp,10^-10)), parms =theta_sp, times = Times_obs)
out_log_transf_sp3 <- ode(func = Microbiota_sim_log_transf_rest, y = log(pmax(x0_int_mic3_sp,10^-10)), parms =theta_sp, times = Times_obs)

mean_state_magn1 = colSums(abs(out_log_transf_sp1[,2:8]))/length(Times_obs)
mean_state_magn2 = colSums(abs(out_log_transf_sp2[,2:8]))/length(Times_obs)
mean_state_magn3 = colSums(abs(out_log_transf_sp3[,2:8]))/length(Times_obs)

 Obs_Y1 =out_log_transf_sp1[,2:8] + sigma_obs*matrix(rnorm(n = nb_obs*dim_syst),nb_obs,dim_syst)*matrix(mean_state_magn1,nb_obs,dim_syst)
 Obs_Y2 =out_log_transf_sp2[,2:8] + sigma_obs*matrix(rnorm(n = nb_obs*dim_syst),nb_obs,dim_syst)*matrix(mean_state_magn2,nb_obs,dim_syst)
 Obs_Y3 =out_log_transf_sp3[,2:8] + sigma_obs*matrix(rnorm(n = nb_obs*dim_syst),nb_obs,dim_syst)*matrix(mean_state_magn3,nb_obs,dim_syst)
 
 Time_list = append(list(Times_obs),list(Times_obs))
 Time_list = append(Time_list,list(Times_obs))
  
  Obs_Y_multi = append(list(t(Obs_Y1)),list(t(Obs_Y2)))
  Obs_Y_multi = append(Obs_Y_multi,list(t(Obs_Y3)))
  
  State_ini_multi = append(list(rbind(t(Obs_Y1),matrix(1,1,nb_obs))),list(rbind(t(Obs_Y2),matrix(1,1,nb_obs))))
  State_ini_multi = append(State_ini_multi,list(rbind(t(Obs_Y3),matrix(1,1,nb_obs))))
  
  ##### Parameter estimation
 
  # perturbation matrix
  mat_B = rbind(diag(dim_control),matrix(0,1,7))
  
  # observation matrix
  mat_C = cbind(diag(dim_syst),matrix(0,7,1))
 
  # trial of tested weighting parameters
  lambda_trial = c(1,2,5)
  
  #estimation with profiling on ic: Results are embedded in object "res_est_err_pred_opt_prof_ic" and "list_res_estimation_prof_ic"
  # -res_est_err_pred_opt_prof_ic$est_param: the estimated value
  # -res_est_err_pred_opt_prof_ic$err_pred: the prediction error 
  # -res_est_err_pred_opt_prof_ic$set_Ricatti_State_opt_value: a list composed of:
  #       -res_est_err_pred_opt_prof_ic$set_Ricatti_State_opt_value$Val_cost_prof: The minimal cost function value 
  #       -res_est_err_pred_opt_prof_ic$set_Ricatti_State_opt_value$control_opt: The optimal control  corresponding to the estimation computed at discretization points for each subject
  #       -res_est_err_pred_opt_prof_ic$set_Ricatti_State_opt_value$State_opt: The optimal trajectory corresponding to the estimation computed at discretization points for each subject
  #       -res_est_err_pred_opt_prof_ic$set_Ricatti_State_opt_value$State_obs: The optimal trajectory corresponding to the estimation computed at observation points for each subject
  #list_res_estimation_prof_ic is a list made of each res_est_err_pred_opt_prof_ic element computed for each lambda value
  param_specified_inter = c(vect_g_sp, vect_s_sp, matrix_M_sp[6,])
  param_specified_sp3 = c(param_specified_inter, theta_ini_inter[c(2,6,14,15,18,20,21,26,27,29,39)])
  Microbiota_matA = function(t,State,param){Microbiota_matA_log_transf_rest_multi(t,State,param,param_specified_sp3)}
  
  err_pred_opt_cur = 10^20
  theta_prof_ini = theta_ini_sp3
  list_res_estimation_prof_ic= list()
  for (lambda  in lambda_trial){
    
    mat_U = lambda*diag(dim_control)
    
    res_est = est_param_iter_oca_nonlinear_profiled_ic_multi_curves(Time_list,Obs_Y_multi,State_ini_multi,param_ini=theta_prof_ini
                                                                            ,mat_U,mat_A = Microbiota_matA ,mat_B,mat_C,mesh_mult=20,x0_known = 1,type_opt =0)
    
    list_res_estimation_prof_ic  = append(list_res_estimation_prof_ic ,list(res_est))
    
    if (res_est$err_pred < err_pred_opt_cur){
      err_pred_opt_cur = res_est$err_pred
      res_est_err_pred_opt_prof_ic = res_est
      theta_prof_ini  = res_est_err_pred_opt_prof_ic$est_param
    }
  }
 