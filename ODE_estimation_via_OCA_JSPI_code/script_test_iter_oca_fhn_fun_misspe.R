pathnames <- list.files(pattern="[.]R$", path="function_model_and_sims//", full.names=TRUE);
sapply(pathnames, FUN=source);
pathnames <- list.files(pattern="[.]R$", path="function_util_estimation//", full.names=TRUE);
sapply(pathnames, FUN=source);
library('deSolve')
library('optimx')

######## Data generation
# observation interval & true parameter specification

fun_a = function(t){0.2*(1+sin(t/5))}
b = 0.2
c = 3
theta = c(b,c)
theta_to_est = c(b,c)
x0= c(-1,1)

Times_obs <- seq(0, 20, by = 20/50)
Times_integ <- seq(0, 20, by = 20/1000)
nb_obs= length(Times_obs)
dim_syst = 2
dim_obs = 1

# noise level specification and data generation
sigma_obs = 0.03

# level of perturbation tested in the article
sig_r = sqrt(0.1)
#sig_r = sqrt(0.15)

out_pert_ref = fhn_stoch_pert_sim (deb=0,h=Times_obs[length(Times_obs)]/5000,fin=Times_obs[length(Times_obs)],Pars = theta,fun_a,y0= x0,sd_stoch_pert = sig_r)
out_pert = t(out_pert_ref[[2]][,seq(1,length(out_pert_ref[[1]]),5000/50)])

Obs_Y =  matrix(data = rep(0, nb_obs*dim_obs), nrow = nb_obs, ncol = dim_obs )
for (col_i in 1:dim_obs) Obs_Y[,col_i] = out_pert[,col_i]+sigma_obs*rnorm(n = nb_obs)
Obs_Y = t(Obs_Y)
State_ini = rbind(Obs_Y,matrix(1,3,nb_obs))

##### Parameter estimation

dim_control = 3

# perturbation matrix
mat_C = matrix(c(1,0,0,0),1,4)

# observation matrix
mat_B = rbind(c(1,0,0),c(0,1,0),c(0,0,0),c(0,0,1))

# trial of tested weighting matrices
lambda_trial_tot = list()
lambda_trial_tot = append(lambda_trial_tot,list(rbind(c(0.005, 0, 0) ,c(0, 0.005, 0),c(0, 0, 0.005))))
lambda_trial_tot = append(lambda_trial_tot,list(rbind(c(0.01, 0, 0) ,c(0, 0.01, 0),c(0, 0, 0.01))))
lambda_trial_tot = append(lambda_trial_tot,list(rbind(c(0.05, 0, 0) ,c(0, 0.05, 0),c(0, 0, 0.05))))
lambda_trial = lambda_trial_tot

#estimation with profiling on ic: Results are embedded in object "res_est_err_pred_opt_prof_ic" and "list_res_estimation_prof_ic"
# res_est_err_pred_opt_prof_ic is a list composed of:
# -res_est_err_pred_opt_prof_ic$est_param: the estimated value
#- res_est_err_pred_opt_prof_ic$err_pred: the minimal prediction error 
# -res_est_err_pred_opt_prof_ic$set_Ricatti_State_opt_value: a list composed among other things of:
#       -res_est_err_pred_opt_prof_ic$set_Ricatti_State_opt_value$control_opt: The optimal control  corresponding to the estimation
#       -res_est_err_pred_opt_prof_ic$set_Ricatti_State_opt_value$State_opt: The optimal trajectory corresponding to the estimation
#list_res_estimation_prof_ic is a list made of each res_est_err_pred_opt_prof_ic element computed for each lambda value


err_pred_opt_cur = 10^20
list_res_estimation = list()
theta_prof_ini = theta_to_est
for (lambda  in lambda_trial){
  
  mat_U = lambda
  
  res_est= est_param_iter_oca_nonlinear_profiled_ic(Times_obs,Obs_Y,State_ini,
                                                    param_ini=theta_prof_ini,mat_U,mat_A=fhn_fun_a_matA,mat_B,mat_C,mesh_mult=20,x0_known = c())
  
  
  list_res_estimation = append(list_res_estimation,list(res_est))
  
  if (res_est$err_pred < err_pred_opt_cur){
    err_pred_opt_cur = res_est$err_pred
    res_est_err_pred_opt_prof_ic = res_est
    theta_prof_ini = res_est_err_pred_opt_prof_ic$est_param
  }
}
#Plot estimated functional parameter in black and true one in red
plot(Times_integ,res_est_err_pred_opt_prof_ic$set_Ricatti_State_opt_value$State_opt[3,],type='l',col="black",xlab = "Time",ylab = "t ->a(t)",main="Functional estimation with IC profiling")
lines(Times_integ,fun_a(Times_integ),col="red")
legend("topright",c("True a","Est. a"),pch=c(NA,NA),lty=c(1,1),col=c("red", "black"))

#estimation without profiling on ic: Results are embedded in object "res_est_err_pred_opt" and "list_res_estimation"
# res_est_err_pred_opt is a list composed of:
# -res_est_err_pred_opt$est_param: the estimated value
#- res_est_err_pred_opt$err_pred: the minimal prediction error 
# -res_est_err_pred_opt$set_Ricatti_State_opt_value: a list composed among other things of:
#       -res_est_err_pred_opt$set_Ricatti_State_opt_value$control_opt: The optimal control  corresponding to the estimation
#       -res_est_err_pred_opt$set_Ricatti_State_opt_value$State_opt: The optimal trajectory corresponding to the estimation
#list_res_estimation is a list made of each res_est_err_pred_opt element computed for each lambda value

err_pred_opt_cur = 10^20
list_res_estimation = list()
theta_ini = c(theta_to_est,x0,c(0,0))
for (lambda  in lambda_trial){
  
  mat_U = lambda
  
  res_est= est_param_iter_oca_nonlinear(Times_obs,Obs_Y,
                                        param_ini=theta_ini,x0_known = c(),mat_U,mat_A=fhn_fun_a_matA,mat_B,mat_C,mesh_mult=20)
  
  list_res_estimation = append(list_res_estimation,list(res_est))
  
  if (res_est$err_pred < err_pred_opt_cur){
    err_pred_opt_cur = res_est$err_pred
    res_est_err_pred_opt = res_est
    theta_ini = res_est_err_pred_opt$est_param
  }
}

#Plot estimated functional parameter in black and true one in red
plot(Times_integ,res_est_err_pred_opt$set_Ricatti_State_opt_value$State_opt[3,],type='l',col="black",xlab = "Time",ylab = "t ->a(t)",main="Functional estimation with IC estimation")
lines(Times_integ,fun_a(Times_integ),col="red")
legend("topright",c("True a","Est. a"),pch=c(NA,NA),lty=c(1,1),col=c("red", "black"))
