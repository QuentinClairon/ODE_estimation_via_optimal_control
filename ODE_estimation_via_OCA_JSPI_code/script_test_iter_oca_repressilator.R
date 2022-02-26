###### Repressilator Simulation: Well-speicfied #########
pathnames <- list.files(pattern="[.]R$", path="function_model_and_sims//", full.names=TRUE);
sapply(pathnames, FUN=source);
pathnames <- list.files(pattern="[.]R$", path="function_util_estimation//", full.names=TRUE);
sapply(pathnames, FUN=source);

library('snow')
library('deSolve')
library('optimx')


######## Data generation
# observation interval & true parameter speicfication

x0_mrna = c(60,20,6)
x0_pro = c(18,27,1)
x0 = c(x0_mrna,x0_pro)

v_max= c(50, 100, 80)
k_max = c(50,30,40)
k_g = c(1,1,1)
k_mrna = c(5,6,7)
k_p = c(1,2,3)
n = 3

theta = c(v_max,k_max,k_g,k_mrna,k_p,n)
theta_to_est = c(v_max,k_max[1:2],k_g,k_p)

Times_obs = seq(0, 20, by = 20/25)
nb_obs <- length(Times_obs)
dim_syst = length(x0)
dim_obs = 3

# noise level speicfication and data generation

# level of noise tested in the article
sigma_obs = 1
# sigma_obs = 1.5
# sigma_obs = 2


out <- ode(func = Repressilator_sim, y = x0, parms =theta, times = Times_obs)
X_without_noise = out[,2:7]

Obs_Y =  matrix(data = rep(0, nb_obs*dim_obs), nrow = nb_obs, ncol = dim_obs)
for (col_i in 1:dim_obs) Obs_Y[,col_i] =  pmax(X_without_noise[,col_i]+sigma_obs*rnorm(n = nb_obs),matrix(0,nb_obs,1))
 Obs_Y = t(Obs_Y)
 
 State_ini = rbind(Obs_Y,matrix(x0_pro[1],1,length(Times_obs)),matrix(x0_pro[2],1,length(Times_obs)),
                   matrix(x0_pro[3],1,length(Times_obs)), matrix(1,1,length(Times_obs)))
 
 
 ##### Parameter estimation
 
 dim_control = 6

 # perturbation matrix
 mat_B = rbind(diag(dim_control),c(0,0,0,0,0,0))
 
 # observation matrix
 mat_C = rbind(cbind(diag(dim_obs),matrix(0,3,4)))
 
 # trial of tested weighting parameters
 lambda_trial = c(50,100,200)
 
 #estimation with profiling on ic: Results are embedded in object "res_est_err_pred_opt_prof_ic" and "list_res_estimation_prof_ic"
 # res_est_err_pred_opt_prof_ic is a list composed of:
 # -res_est_err_pred_opt_prof_ic$est_param: the estimated value
 #- res_est_err_pred_opt_prof_ic$err_pred: the minimal prediction error 
 # -res_est_err_pred_opt_prof_ic$set_Ricatti_State_opt_value: a list composed among other things of:
 #       -res_est_err_pred_opt_prof_ic$set_Ricatti_State_opt_value$control_opt: The optimal control  corresponding to the estimation
 #       -res_est_err_pred_opt_prof_ic$set_Ricatti_State_opt_value$State_opt: The optimal trajectory corresponding to the estimation
 #list_res_estimation_prof_ic is a list made of each res_est_err_pred_opt_prof_ic element computed for each lambda value
 
 err_pred_opt_cur = 10^20
 theta_ini_prof_ini = theta_to_est
 list_res_estimation_prof_ic  = list()
 for (lambda  in lambda_trial){
   
   mat_U = lambda*diag(dim_control)
   res_est = est_param_iter_oca_nonlinear_profiled_ic(Times_obs,Obs_Y,State_ini,
                                                      param_ini=theta_ini_prof_ini,mat_U,mat_A=Repressilator_matA,mat_B,mat_C,mesh_mult=20,x0_known =  c(x0_pro,1),type_opt =1)
   
   
   list_res_estimation_prof_ic  = append(list_res_estimation_prof_ic ,list(res_est))
   
   if (res_est$err_pred < err_pred_opt_cur){
     err_pred_opt_cur = res_est$err_pred
     res_est_err_pred_opt_prof_ic = res_est
     theta_ini_prof_ini  = res_est_err_pred_opt_prof_ic$est_param
   }
   
 }
 
 #estimation without profiling on ic: Results are embedded in object "res_est_err_pred_opt" and "list_res_estimation"
 # res_est_err_pred_opt is a list composed of:
 # -res_est_err_pred_opt$est_param: the estimated value
 #- res_est_err_pred_opt$err_pred: the minimal prediction error 
 # -res_est_err_pred_opt$set_Ricatti_State_opt_value: a list composed among other things of:
 #       -res_est_err_pred_opt$set_Ricatti_State_opt_value$control_opt: The optimal control  corresponding to the estimation
 #       -res_est_err_pred_opt$set_Ricatti_State_opt_value$State_opt: The optimal trajectory corresponding to the estimation
 #list_res_estimation is a list made of each res_est_err_pred_opt element computed for each lambda value
 
 err_pred_opt_cur = 10^20
 theta_ini = c(theta_to_est,x0_mrna)
 list_res_estimation = list()
 for (lambda  in lambda_trial){
   
   mat_U = lambda*diag(dim_control)
   res_est = est_param_iter_oca_nonlinear(Times_obs,Obs_Y,
                                                param_ini=theta_ini,x0_known =  c(x0_pro,1),mat_U,mat_A=Repressilator_matA,mat_B,mat_C,mesh_mult=20)
   
   list_res_estimation = append(list_res_estimation,list(res_est))
   
   if (res_est$err_pred < err_pred_opt_cur){
     err_pred_opt_cur = res_est$err_pred
     res_est_err_pred_opt = res_est
     theta_ini = res_est_err_pred_opt$est_param
   }
   
 }

