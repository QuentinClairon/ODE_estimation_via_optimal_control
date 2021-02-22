compute_Riccati_nonlinear_profiled_ic<- function(Times_integ,pseudo_Y,State_i_m1,param_struct,mat_U,mat_A,mat_B,mat_C,weight_integ,x0_known){
  ############################################################
  # DESCRIPTION: compute the Riccati  solution at a given iteration l of the algorithm described in section 3.3 when profiling on CI and for nonlinear ODE models
  # 
  # FUNCTION SIGNATURE:  compute_Riccati_nonlinear_profiled_ci<- function(Times_integ,pseudo_Y,State_i_m1,param_struct,mat_U,mat_A,mat_B,mat_C,weight_integ,x0_known)
  # 
  # INPUT :
  # Times_integ               :(vector) discretization time points
  # pseudo_Y                  :(dim_obs x nb_time_integ sized matrix) extended data as described in section 3.1
  # State_i_m1                :(dim_syst x nb_time_integ sized matrix) computed optimal trajectory at iteration (l-1)
  # param_struct              :(nb_param sized vector)  current assumed value for theta 
  # mat_U                     :(dim_control x dim_control sized matrix) weighting matrix for control magnitude penalization
  # mat_A                     :(function handler) matrix valued function representing the pseudolinear model dX/dt = A(t,X,theta)X of signature: 
  #                                              mat_A<- function(t,State,param) 
  #                                                INPUT: t (real number) current time value
  #                                                     : State (dim_syst sized vector) state_variable value in the pseudo linear representation
  #                                                     : param (vector) current parameter value
  #                                                OUTPUT: (dim_syst square matrix) the value of  A(t,State,param) 
  # mat_B                     :(dim_syst x dim_control sized matrix) perturbation matrix describing how perturbations act on the original ODE
  # mat_C                     :(dim_obs x dim_syst sized matrix) observation matrix describing what is observed among the ODE state variables
  # weight_integ           :(nb_time_integ sized vector) weights w_j as described in section 3.1
  # OPTIONAL INPUT:
  # x0_known                   :(vector, default value, empty vector ) list of known initial conditions, for the sake of convenience the model 
  #                              should be written such that the corresponding states variables are written last
  # OUTPUT:
  # Return a list set_Ricatti_State_opt_value composed of the following elements computed at iteration l:
  #       -set_Ricatti_State_opt_value$Val_cost_prof: The minimal cost function value 
  #       -set_Ricatti_State_opt_value$control_opt: The optimal control computed at discretization points
  #       -set_Ricatti_State_opt_value$State_opt: The optimal trajectory computed at discretization points
  #       -set_Ricatti_State_opt_value$State_obs: The optimal trajectory computed at observation points
  #       -set_Ricatti_State_opt_value$Ricatti: The Riccati equation solution computed at discretization points
  #       -set_Ricatti_State_opt_value$Adjoint: The adjoint equation solution computed at discretization points
  #       -set_Ricatti_State_opt_value$val_RSS: Residual sum of square between the perturbed ODE solutions and observation
  
  ############################################################
  
  dim_syst = ncol(mat_C)
  dim_obs = nrow(mat_C)
  dim_control = ncol(mat_U)
  nb_time_integ = length(Times_integ)
  Yn_fin = pseudo_Y[,nb_time_integ];
  nb_unknown = dim_syst - length(x0_known)
  
  mat_Rn  = t(mat_C)%*%mat_C
  vect_Hn = -t(mat_C)%*%Yn_fin 
  History_R = list()
  History_H = list()
  History_R = list(mat_Rn) 
  History_H = list(vect_Hn) 
  History_G = list()
  mat_Rk  = mat_Rn
  vect_Hk = vect_Hn
  
  val_Sn_prof = t(Yn_fin)%*%Yn_fin
  
  for (nk in 1:(nb_time_integ-1)){
    delta_k = Times_integ[nb_time_integ-nk+1] - Times_integ[nb_time_integ-nk] 
    Yk= pseudo_Y[,nb_time_integ-nk]
    
    mat_Ak = mat_A(Times_integ[nb_time_integ-nk],State_i_m1[,nb_time_integ-nk],param_struct)
    mat_Ak_delta = mat_Ak*delta_k + diag(dim_syst)
    weight_cur = weight_integ[nb_time_integ-nk]
    
#   print(mat_U+delta_k*t(mat_B)%*%mat_Rk%*%mat_B)
    mat_G = solve(mat_U+delta_k*t(mat_B)%*%mat_Rk%*%mat_B)
    History_G = append(list(mat_G),History_G)
    
    val_Sn_prof = val_Sn_prof +delta_k*(weight_cur*t(Yk)%*%Yk - t(vect_Hk)%*%mat_B%*%mat_G%*%t(mat_B)%*%vect_Hk)    
    
   
    mat_Rk_m1_1part = mat_Rk + delta_k*weight_cur*t(mat_C)%*%mat_C + delta_k*(mat_Rk%*%mat_Ak+t(mat_Ak)%*%mat_Rk)+(delta_k^2)*t(mat_Ak)%*%mat_Rk%*%mat_Ak
    mat_Rk_m1_2part = -delta_k*t(mat_Ak_delta)%*%mat_Rk%*%mat_B%*%mat_G%*%t(mat_B)%*%mat_Rk%*%mat_Ak_delta
    mat_Rk_m1 = mat_Rk_m1_1part+mat_Rk_m1_2part
    
    vect_Hk_m1 =  vect_Hk - delta_k*weight_cur*t(mat_C)%*%Yk+delta_k*t(mat_Ak)%*%vect_Hk-  delta_k*t(mat_Ak_delta)%*%mat_Rk%*%mat_B%*%mat_G%*%t(mat_B)%*%vect_Hk 
    
    mat_Rk = mat_Rk_m1
    vect_Hk = vect_Hk_m1
    

    History_R = append(list(mat_Rk),History_R)
    History_H = append(list(vect_Hk),History_H)
   
  }
  
 
  mat_R0  =  History_R[[1]]
  vect_H0 = History_H[[1]]
  
 

  if (is.null(x0_known)){
    est_x0 = -solve(mat_R0)%*%vect_H0
  }else{
   
    
    mat_R0_11 = mat_R0[seq(1,nb_unknown),seq(1,nb_unknown)]
    mat_R0_12 = mat_R0[seq(1,nb_unknown),seq(nb_unknown+1,dim_syst)]
    
    

    if (length(x0_known) ==1){
    est_x0_unknown = -solve(mat_R0_11)%*%(mat_R0_12*x0_known+vect_H0[seq(1,nb_unknown)])
    }else{
    est_x0_unknown = -solve(mat_R0_11)%*%(mat_R0_12%*%x0_known+vect_H0[seq(1,nb_unknown)])
    }
    est_x0 = c(t(est_x0_unknown),t(x0_known))
    
  }
  val_Sn_prof  = val_Sn_prof + t(est_x0)%*%mat_R0%*%est_x0+2*t(vect_H0)%*%est_x0 
  
  
  State_i = matrix(data = rep(0,dim_syst*nb_time_integ ), nrow = dim_syst, ncol = nb_time_integ )
  State_i[,1] =  est_x0
  
  control_i =matrix(data = rep(0,dim_control*(nb_time_integ-1)), nrow = dim_control, ncol = nb_time_integ-1 )
  val_RSS = 0
  for (nkk in 1:(nb_time_integ-1)){
    delta_k = Times_integ[nkk+1] - Times_integ[nkk] 
    mat_Ak = mat_A(Times_integ[nkk],State_i_m1[,nkk],param_struct)
    mat_Ak_delta = mat_Ak*delta_k + diag(dim_syst)
    
    mat_Rk_p1 = History_R[[nkk+1]]
    vect_Hk_p1 = History_H[[nkk+1]]
    mat_G = History_G[[nkk]]
    control_ukk = -mat_G%*%t(mat_B)%*%(mat_Rk_p1%*%mat_Ak_delta%*%State_i[,nkk]+vect_Hk_p1)
    control_i[,nkk] = control_ukk
    State_i[,nkk+1]  = mat_Ak_delta%*%State_i[,nkk]+delta_k*mat_B%*%control_ukk
    
    if (weight_integ[nkk] > 0){
      val_RSS = val_RSS+ sum((mat_C%*%State_i[,nkk] -  pseudo_Y[,nkk])^2);
    }
  }
  
  val_RSS  = val_RSS+ sum((mat_C%*%State_i[,nb_time_integ] - Yn_fin)^2);
  indice_obs  =weight_integ>0.2      
  State_i_obs_time = State_i[,indice_obs]
  
  ensemble_est = list(Val_cost_prof = val_Sn_prof ,control_opt =control_i, Ricatti = History_R,Adjoint =  History_H,State_opt =State_i,Sequence_G = History_G,State_obs = State_i_obs_time,val_RSS=val_RSS)
  
  return(ensemble_est)
}