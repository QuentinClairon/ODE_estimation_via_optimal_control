pathnames <- list.files(pattern="[.]R$", path="function_model_and_sims//", full.names=TRUE);
sapply(pathnames, FUN=source);
pathnames <- list.files(pattern="[.]R$", path="function_util_estimation//", full.names=TRUE);
sapply(pathnames, FUN=source);

library('snow')
library('deSolve')
library('optimx')

######## Data generation
# observation interval & true parameter specification

x0 <- c(100,0,0,0,0)
theta <- c(5.93,2.96,2.05,27.5,4)*10^-2

step = 10
beg_obs = 0
end_obs = 100
Times_obs <- seq(beg_obs,end_obs, by = 100/step )

nb_obs <- length(Times_obs)
dim_syst = length(x0)

# noise level specification and data generation
# level of perturbation tested in the article
sig_c = sqrt(0.004)
#sig_c = sqrt(0.006)
#sig_c = sqrt(0.008)

sigma_obs = 5

out_pert_ref = Apinene_stoch_pert_sim (deb=0,h=Times_obs[length(Times_obs)]/5000,fin=Times_obs[length(Times_obs)],Pars = theta,y0= x0,sd_stoch_pert = sig_c)
out_pert = t(out_pert_ref[[2]][,seq(1,length(out_pert_ref[[1]]),5000/step)])

sd_per_step_var = step*colSums(out_pert)*sigma_obs/((end_obs-beg_obs)*100)

Obs_Y =  matrix(data = rep(0, nb_obs*dim_syst), nrow = nb_obs, ncol = dim_syst)
for (col_i in 1:dim_syst) Obs_Y[,col_i] = out_pert[,col_i]+sd_per_step_var[col_i]*rnorm(n = nb_obs)
Obs_Y = t(Obs_Y)


##### Parameter estimation

dim_control = 5

# perturbation matrix
mat_B = diag(dim_control)

# observation matrix
mat_C = diag(dim_syst)


# trial of tested weighting parameters
lambda_trial = c(1,5,10,50,100,500)


#estimation whith profiling on ic: Results are embedded in object "res_est_err_pred_opt_prof_ic" and "list_res_estimation_prof_ic"
# res_est_err_pred_opt_prof_ic is a list composed of:
# -res_est_err_pred_opt_prof_ic$est_param: the estimated value
#- res_est_err_pred_opt_prof_ic$err_pred: the minimal prediction error 
# -res_est_err_pred_opt_prof_ic$set_Ricatti_State_opt_value: a list composed among other things of:
#       -res_est_err_pred_opt_prof_ic$set_Ricatti_State_opt_value$control_opt: The optimal control  corresponding to the estimation
#       -res_est_err_pred_opt_prof_ic$set_Ricatti_State_opt_value$State_opt: The optimal trajectory corresponding to the estimation
#list_res_estimation_prof_ic is a list made of each res_est_err_pred_opt_prof_ic element computed for each lambda value
err_pred_opt_cur = 10^20

theta_prof_ic_ini = theta
list_res_estimation_prof_ic = list()
for (lambda  in lambda_trial){
  
  mat_U = lambda*diag(dim_control)
  
  res_est = est_param_iter_oca_linear_profiled_ic(Times_obs,Obs_Y,param_ini=theta_prof_ic_ini,mat_U,mat_A = Apinene_matA,mat_B,mat_C,mesh_mult=50)
  list_res_estimation_prof_ic = append(list_res_estimation_prof_ic,list(res_est))
  
  if (res_est$err_pred < err_pred_opt_cur){
    err_pred_opt_cur = res_est$err_pred
    res_est_err_pred_opt_prof_ic = res_est
    theta_prof_ic_ini = res_est_err_pred_opt_prof_ic$est_param
  }
  
}


#estimation whithout profiling on ic: Results are embedded in object "res_est_err_pred_opt" and "list_res_estimation"
# res_est_err_pred_opt is a list composed of:
# -res_est_err_pred_opt$est_param: the estimated value
#- res_est_err_pred_opt$err_pred: the minimal prediction error 
# -res_est_err_pred_opt$set_Ricatti_State_opt_value: a list composed among other things of:
#       -res_est_err_pred_opt$set_Ricatti_State_opt_value$control_opt: The optimal control  corresponding to the estimation
#       -res_est_err_pred_opt$set_Ricatti_State_opt_value$State_opt: The optimal trajectory corresponding to the estimation
#list_res_estimation is a list made of each res_est_err_pred_opt element computed for each lambda value
err_pred_opt_cur = 10^20

theta_ini = c(theta,x0)
list_res_estimation = list()
for (lambda  in lambda_trial){
  
  mat_U = lambda*diag(dim_control)
  
  
  res_est = est_param_iter_oca_linear(Times_obs,Obs_Y,param_ini=theta_ini,x0_known = c(),mat_U,mat_A = Apinene_matA,mat_B,mat_C,mesh_mult=50)
  list_res_estimation = append(list_res_estimation,list(res_est))
  
  if (res_est$err_pred < err_pred_opt_cur){
    err_pred_opt_cur = res_est$err_pred
    res_est_err_pred_opt = res_est
    theta_ini = res_est_err_pred_opt$est_param
  }
}
