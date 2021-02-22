err_pred_estimation_oca_nonlin <- function(Times_integ,pseudo_Y,est_param,State_opt,mat_U,mat_A,mat_B,mat_C,weight_integ){
  ############################################################
  # DESCRIPTION: Compute the prediction error according to the formula given in section 5 and for nonlinear ODE models
  # 
  # FUNCTION SIGNATURE:  err_pred_estimation_oca_nonlin <- function(Times_integ,pseudo_Y,est_param,State_opt,mat_U,mat_A,mat_B,mat_C,weight_integ)
  # 
  # INPUT :
  # Times_integ               :(vector) discretization time points
  # pseudo_Y                  :(dim_obs x nb_time_integ sized matrix) extended data as described in section 3.1
  # est_param                 :(nb_param sized vector) estimated parameter value 
  # State_opt                 :(dim_syst x nb_time_integ sized matrix) optimal trajectory corresponding to est_param
  # mat_U                     :(dim_control x dim_control sized matrix) weighting matrix for control magnitude penalization
  # mat_A                     :(function handler) matrix valued function representing the pseudolinear model dX/dt = A(t,X,theta)X of signature: 
  #                                              mat_A<- function(t,State,param) 
  #                                                INPUT: t (real number) current time value
  #                                                     : State (dim_syst sized vector) state_variable value in the pseudo linear representation
  #                                                     : param (vector) current parameter value
  #                                                OUTPUT: (dim_syst square matrix) the value of  A(t,State,param) 
  # mat_B                     :(dim_syst x dim_control sized matrix) perturbation matrix describing how perturbations act on the original ODE
  # mat_C                     :(dim_obs x dim_syst sized matrix) observation matrix describing what is observed among the ODE state variables
  # weight_integ             :(nb_time_integ sized vector) weights w_j as described in section 3.1
  # OUTPUT:
  # Return  error_pred  the prediction error computed according to ERRPRED(lambda) formula
  ############################################################
  
  nb_subject =length(Times_integ)
  dim_syst = ncol(mat_C)
  dim_obs = nrow(mat_C)
  dim_control = ncol(mat_U)
  
  v_field <- function(Time, State, Pars){
    dState = mat_A(Time,State,est_param)%*%State
    
    return(list(c(dState)))
  }
  
  error_pred = 0
  
  Times_obs = Times_integ[weight_integ>0]
  Y_obs = subset(t(pseudo_Y),weight_integ>0)
  Y_obs = t(Y_obs)
  
  
  State_obs = subset(t(State_opt),weight_integ>0)
  State_obs = t(State_obs)
  
  v_field_nbs   =function(Time, State, Pars){v_field(Time, State, Pars)}
  
  middle_interval = ceiling(length(Times_obs)/2)
  
  Times_obs_first_part = Times_obs[1:(middle_interval-1)]
  Times_obs_second_part = Times_obs[middle_interval:length(Times_obs)]
  
  Y_obs_first_part = Y_obs[,1:(middle_interval-1)]
  Y_obs_second_part = Y_obs[,middle_interval:length(Times_obs)]
  
  X0_opt =  State_obs[,1]
  Xmiddle_opt = State_obs[,middle_interval]
  
  
  pred_first_part <- ode(func =v_field, y = X0_opt, parms =c(), times = Times_obs_first_part)
  pred_second_part <- ode(func =v_field, y = Xmiddle_opt, parms =c(), times = Times_obs_second_part)
  
  error_pred = error_pred+sum((Y_obs_first_part - mat_C%*%t(pred_first_part[,2:(dim_syst+1)]))^2)
  error_pred = error_pred+sum((Y_obs_second_part - mat_C%*%t(pred_second_part[,2:(dim_syst+1)]))^2)
  
  
  return(error_pred)
}