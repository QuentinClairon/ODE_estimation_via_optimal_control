est_param_iter_oca_nonlinear_profiled_ic_multi_curves <- function(Times_obs_list,Yn_list,State_ini_list,param_ini,mat_U,mat_A,mat_B,mat_C,mesh_mult=1,x0_known=c(),type_opt=1){
  ############################################################
  # DESCRIPTION: compute the parameter estimation from observations coming from different subjects via the method presented in article when profiling on ic  and for nonlinear ODE models 
  # 
  # FUNCTION SIGNATURE:  est_param_iter_oca_nonlinear_profiled_ic_multi_curves <- function(Times_obs_list,Yn_list,State_ini_list,param_ini,mat_U,mat_A,mat_B,mat_C,mesh_mult=1,x0_known=c(),type_opt=1)
  # 
  # INPUT :
  # Times_obs_lis             :(list of vector) list of observation time points for each subjects
  # Yn                        :(list of dim_obs x nb_obs sized matrix) list of observations matrix  for each subjects
  # State_ini                 :(list of dim_syst x nb_obs sized matrix) list of  initial guess for the system state  for each subjects
  # param_ini                 :(nb_param sized vector)  initial guess for theta estimator
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
  #x0_known                   :(vector, default value, empty vector ) list of known initial conditions, for the sake of convenience the model 
  #                              should be written such that the corresponding states variables are written last
  #type_opt                   :(integer, default value =1), speicfy the type of optimization algorithm used for estimation (if 1 then "nlminb" or Nelder-Mead algorithm otherwise)
  
  # OUTPUT:
  # Return a list  res_est  composed of:
  # -res_est$est_param: the estimated value
  #- res_est$err_pred: the prediction error 
  # -res_est$set_Ricatti_State_opt_value: a list composed of the element given for the last iteration of the algorithm described in section 3.3:
  #       -res_est$set_Ricatti_State_opt_value$Val_cost_prof: The minimal cost function value 
  #       -res_est$set_Ricatti_State_opt_value$control_opt: The optimal control  corresponding to the estimation computed at discretization points for each subject
  #       -res_est$set_Ricatti_State_opt_value$State_opt: The optimal trajectory corresponding to the estimation computed at discretization points for each subject
  #       -res_est$set_Ricatti_State_opt_value$State_obs: The optimal trajectory corresponding to the estimation computed at observation points for each subject
  #       -res_est$set_Ricatti_State_opt_value$Ricatti: The Riccati equation solution corresponding to the estimation computed at discretization points for each subject
  #       -res_est$set_Ricatti_State_opt_value$Adjoint: The adjoint equation solution corresponding to the estimation computed at discretization points for each subject
  
  ############################################################
  
  dim_syst = ncol(mat_C)
  dim_obs = nrow(mat_C)
  
  dim_control = ncol(mat_U)
  nb_obs = length(Times_obs_list[[1]])
  
  Times_integ_list = list()
  pseudo_Y_list = list()
  weight_integ_list = list()
  pseudo_State_ini_list= list()
  
  for (nbs in 1:(length(Times_obs_list))){
  Yn =  Yn_list[[nbs]] 
  Times_obs = Times_obs_list[[nbs]]
  State_ini = State_ini_list[[nbs]]
  Yn_fin = Yn[,nb_obs];
  
  nb_time_integ = mesh_mult*(nb_obs-1)+1;
  pseudo_State_ini =  matrix(data = rep(0,dim_syst*nb_time_integ ), nrow = dim_syst, ncol = nb_time_integ )
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
    pseudo_State_ini[,seq(mesh_mult*(obs_i-1)+1,(mesh_mult*obs_i+1),by =1)]= matrix(rep(State_ini[,obs_i],mesh_mult+1),dim_syst,mesh_mult+1)
  }
  
  Times_integ_list  = append(Times_integ_list ,list(Times_integ))
  pseudo_Y_list  = append(pseudo_Y_list,list(pseudo_Y))
  weight_integ_list = append( weight_integ_list,list( weight_integ))
  pseudo_State_ini_list = append(pseudo_State_ini_list,list(pseudo_State_ini))
  }
  
  
  func_eco_E <- function(param){
    Val_sum_cost_prof =0
    ensemble_Ricatti_comp_list = list()
    new_pseudo_State_ini_list = list()
    for (nbs in 1:(length(Times_obs_list))){
    val_err_iter_cur = 10^100;
    stop_criteria = 0;
    nb_iter = 0;
    
    Times_integ = Times_integ_list[[nbs]]
    pseudo_Y = pseudo_Y_list[[nbs]]
    State_i_m1 = pseudo_State_ini_list[[nbs]];
    weight_integ=  weight_integ_list[[nbs]]
    
    Riccati_cur = list();
    
    while (stop_criteria ==0){
      ensemble_Ricatti_comp =  compute_Riccati_nonlinear_profiled_ic(Times_integ,pseudo_Y,State_i_m1 ,param,mat_U,mat_A,mat_B,mat_C,weight_integ,x0_known=x0_known)
      
      val_err_iter = sum((ensemble_Ricatti_comp$State_opt -  State_i_m1)^2)
      State_i_m1 = ensemble_Ricatti_comp$State_opt
      nb_iter = nb_iter+1
      
    #  print(nb_iter)
      #  print( val_err_iter)
      if (val_err_iter <10^-4){
        stop_criteria = 1;
      }
      
      if (nb_iter  > 50){
        stop_criteria = 1;
      }
      
    }
    
    new_pseudo_State_ini_list = append(new_pseudo_State_ini_list,list(ensemble_Ricatti_comp$State_opt))
    ensemble_Ricatti_comp_list = append(ensemble_Ricatti_comp_list,list(ensemble_Ricatti_comp))
    Val_sum_cost_prof = Val_sum_cost_prof+ ensemble_Ricatti_comp$Val_cost_prof
    }
    
    pseudo_State_ini_list = new_pseudo_State_ini_list
    return(list(Val_sum_cost_prof=Val_sum_cost_prof, ensemble_Ricatti_comp_list= ensemble_Ricatti_comp_list,State_opt_set =pseudo_State_ini_list))
  }
  
  func_eco <- function(param){
    func_eco_E_val = func_eco_E(param)
    return(func_eco_E_val$Val_sum_cost_prof)
  }
  
 
  if (type_opt==1){
    res_estim = optimr(param_ini, func_eco, gr=NULL, hess=NULL, lower=-Inf, upper=Inf,method="nlminb",control = list(trace=1))
  }else{
    res_estim = optimr(param_ini, func_eco, gr=NULL, hess=NULL, lower=-Inf, upper=Inf,control = list(trace=1)) 
  }
  
  ensemble_Ricatti_final = func_eco_E(res_estim$par)
  
  est_param = res_estim$par
  State_opt_est_set =  ensemble_Ricatti_final$State_opt_set
  err_pred = 0
  for (nbs in 1:(length(Times_obs_list))){
  err_pred_nbs = err_pred_estimation_oca_nonlin(Times_integ_list[[nbs]],pseudo_Y_list[[nbs]],est_param,State_opt_est_set[[nbs]],mat_U,mat_A,mat_B,mat_C,weight_integ_list[[nbs]])
  err_pred = err_pred+err_pred_nbs
  }
  
  return(list(est_param = est_param,set_Ricatti_State_opt_value = ensemble_Ricatti_final,err_pred = err_pred))
}