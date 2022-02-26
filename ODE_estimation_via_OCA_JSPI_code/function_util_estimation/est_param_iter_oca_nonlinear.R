est_param_iter_oca_nonlinear <- function(Times_obs,Yn,param_ini,x0_known = c(),mat_U,mat_A,mat_B,mat_C,mesh_mult=1){
  ############################################################
  # DESCRIPTION: compute the parameter estimation via the method presented in article when estimating CI and for nonlinear ODE models
  # 
  # FUNCTION SIGNATURE:  est_param_iter_oca_nonlinear <- function(Times_obs,Yn,param_ini,x0_known = c(),mat_U,mat_A,mat_B,mat_C,mesh_mult=1)
  # 
  # INPUT :
  # Times_obs                 :(vector) observation time points
  # Yn                        :(dim_obs x nb_obs sized matrix) observations matrix
  # param_ini                 :(nb_param sized vector)  initial guess for theta estimator
  # x0_known                   :(vector, default value: empty vector ) list of known initial conditions, for the sake of convenience the model 
  #                              should be written such that the corresponding states variables are written last
  # mat_U                     :(dim_control x dim_control sized matrix) weighting matrix for control magnitude penalization
  # mat_A                     :(function handler) matrix valued function representing the pseudolinear model dX/dt = A(t,X,theta)X of signature: 
  #                                              mat_A<- function(t,State,param) 
  #                                                INPUT: t (real number) current time value
  #                                                     : State (dim_syst sized vector) state_variable value in the pseudo linear representation
  #                                                     : param (vector) current parameter value
  #                                                OUTPUT: (dim_syst square matrix) the value of  A(t,State,param) 
  # mat_B                     :(dim_syst x dim_control sized matrix) perturbation matrix describing how perturbations act on the original ODE
  # mat_C                     :(dim_obs x dim_syst sized matrix) observation matrix describing what is observed among the ODE state variables
  # OPTIONAL INPUT:
  #mesh_mult                  :(integer, default value =1 ) the number of discretisation point added between two observation points
   # OUTPUT:
  # Return a list  res_est  composed of:
  # -res_est$est_param: the estimated value
  #- res_est$err_pred: the prediction error 
  # -res_est$set_Ricatti_State_opt_value: a list composed of the element given for the last iteration of the algorithm described in section 3.3:
  #       -res_est$set_Ricatti_State_opt_value$Val_cost_prof: The minimal cost function value 
  #       -res_est$set_Ricatti_State_opt_value$control_opt: The optimal control  corresponding to the estimation computed at discretization points
  #       -res_est$set_Ricatti_State_opt_value$State_opt: The optimal trajectory corresponding to the estimation computed at discretization points
  #       -res_est$set_Ricatti_State_opt_value$State_obs: The optimal trajectory corresponding to the estimation computed at observation points
  #       -res_est$set_Ricatti_State_opt_value$Ricatti: The Riccati equation solution corresponding to the estimation computed at discretization points
  #       -res_est$set_Ricatti_State_opt_value$Adjoint: The adjoint equation solution corresponding to the estimation computed at discretization points
  
  ############################################################
  
  dim_syst = ncol(mat_C)
  dim_obs = nrow(mat_C)
  dim_control = ncol(mat_U)
  nb_obs = length(Times_obs)
  Yn_fin = Yn[,nb_obs];
  
  nb_time_integ = mesh_mult*(nb_obs-1)+1;
  Times_integ =  matrix(data = rep(0,nb_time_integ ), nrow = 1, ncol = nb_time_integ )
 
  pseudo_Y =  matrix(data = rep(0,dim_obs*nb_time_integ ), nrow = dim_obs, ncol = nb_time_integ )
  weight_integ =  matrix(data = rep(0,nb_time_integ-1 ), nrow = 1, ncol = nb_time_integ-1 )
  
  for (obs_i in 1:(nb_obs-1)){
    
    ecart_i = Times_obs[obs_i+1] - Times_obs[obs_i]  
    refined_interval_i =  seq(Times_obs[obs_i],Times_obs[obs_i+1],by = (ecart_i/mesh_mult))
    refined_interval_i[length( refined_interval_i)] = Times_obs[obs_i+1]
    Times_integ[seq(mesh_mult*(obs_i-1)+1,(mesh_mult*obs_i+1),by =1)] = refined_interval_i
  
    weight_integ[mesh_mult*(obs_i-1)+1] =  1/(refined_interval_i[2]-Times_obs[obs_i])
    pseudo_Y[,mesh_mult*(obs_i-1)+1]= Yn[,obs_i]
    pseudo_Y[,mesh_mult*obs_i+1] = Yn[,obs_i+1]
  }

  
  func_eco_E <- function(param){
    
    nb_param_tot = length(param);
    param = t(param);
    
    if (is.null(x0_known)){
      nb_param_struct = nb_param_tot - dim_syst;
      x0 = param[seq(nb_param_struct+1,nb_param_tot,by =1)];
      param_struct = param[seq(1,nb_param_struct,by =1)]; 
    }
    else{
      nb_x0_known = length(x0_known);
      x0_known = t(x0_known);
      nb_x0_to_est = dim_syst - nb_x0_known ;                
      nb_param_struct =  nb_param_tot -  nb_x0_to_est;
      param_struct = param[seq(1,nb_param_struct,by =1)];
      
      x0 = matrix(data = rep(0,dim_syst), nrow = dim_syst, ncol = 1 );
      if ((nb_param_struct+1) >nb_param_tot){
        x0[seq(1,dim_syst,by =1)] = x0_known;
      }
      else {
        x0_to_est = param[seq(nb_param_struct+1,nb_param_tot,by =1)];
        x0[seq(1,nb_x0_to_est,by =1)] = x0_to_est;
        x0[seq((nb_x0_to_est+1),dim_syst,by =1)] = x0_known;
      }
      
    }
    
    pseudo_State_ini =  matrix(data = rep(0,dim_syst*nb_time_integ ), nrow = dim_syst, ncol = nb_time_integ )
    
    for (obs_i in 1:(nb_obs-1)){
      pseudo_State_ini[,seq(mesh_mult*(obs_i-1)+1,(mesh_mult*obs_i+1),by =1)]= matrix(rep(x0,mesh_mult+1),dim_syst,mesh_mult+1)
    }
    
    
    
    State_i_m1 = pseudo_State_ini;
    val_err_iter_cur = 10^100;
    stop_criteria = 0;
    nb_iter = 0;
    Riccati_cur = list();
    
    while (stop_criteria ==0){
      ensemble_Ricatti_comp =  compute_Riccati_nonlinear(Times_integ,pseudo_Y,State_i_m1 ,param,mat_U,mat_A,mat_B,mat_C,weight_integ)    
      val_err_iter = sum((ensemble_Ricatti_comp$State_opt -  State_i_m1)^2)
      State_i_m1 = ensemble_Ricatti_comp$State_opt
      nb_iter = nb_iter+1
      
      
      if (val_err_iter <10^-4){
        stop_criteria = 1;
      }
      
      if (nb_iter  > 50){
        stop_criteria = 1;
      }
      
    }
    
    return(ensemble_Ricatti_comp)
  }
  
  func_eco <- function(param){
    ensemble_Ricatti_comp = func_eco_E(param)
    return(ensemble_Ricatti_comp$Val_cost_prof)
  }
  
  res_estim = optimr(param_ini, func_eco, gr=NULL, hess=NULL, lower=-Inf, upper=Inf,method="nlminb",control = list(trace=1))
  ensemble_Ricatti_final = func_eco_E(res_estim$par)
  
  est_param = res_estim$par
  State_opt_est =  ensemble_Ricatti_final$State_opt
  
  err_pred = err_pred_estimation_oca_nonlin(Times_integ,pseudo_Y,est_param,State_opt_est,mat_U,mat_A,mat_B,mat_C,weight_integ)
  
  return(list(est_param = est_param,set_Ricatti_State_opt_value = ensemble_Ricatti_final,err_pred = err_pred))
}
