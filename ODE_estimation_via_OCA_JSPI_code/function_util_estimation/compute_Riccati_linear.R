compute_Riccati_linear <- function(Times_integ,pseudo_Y,x0,param_struct,mat_U,mat_A,mat_C,mat_B,weight_integ){
  ############################################################
  # DESCRIPTION: compute the solution of the Riccati  when estimating IC and for linear ODE models
  # 
  # FUNCTION SIGNATURE:  compute_Riccati_linear <- function(Times_integ,pseudo_Y,x0,param_struct,mat_U,mat_A,mat_C,mat_B,weight_integ)
  # 
  # INPUT :
  # Times_integ               :(vector) discretization time points
  # pseudo_Y                  :(dim_obs x nb_time_integ sized matrix) extended data as described in section 3.1
  # param_struct              :(nb_param sized vector)  current assumed value for theta 
  # mat_U                     :(dim_control x dim_control sized matrix) weighting matrix for control magnitude penalization
  # mat_A                     :(function handler) matrix valued function representing the linear model dX/dt = A(t,theta)X of signature: 
  #                                               mat_A<- function(t,param) 
  #                                                INPUT: t (real number) current time value
  #                                                     : param (vector) current parameter value
  #                                                OUTPUT: (dim_syst square matrix) the value of  A(t,param) 
  # mat_B                     :(dim_syst x dim_control sized matrix) perturbation matrix describing how perturbations act on the original ODE
  # mat_C                     :(dim_obs x dim_syst sized matrix) observation matrix describing what is observed among the ODE state variables
  # weight_integ           :(nb_time_integ sized vector) weights w_j as described in section 3.1

  # OUTPUT:
  # Return a list set_Ricatti_State_opt_value composed of:
  #       -set_Ricatti_State_opt_value$Val_cost_prof: The minimal cost function value 
  #       -set_Ricatti_State_opt_value$control_opt: The optimal control   computed at discretization points
  #       -set_Ricatti_State_opt_value$State_opt: The optimal trajectory computed at discretization points
  #       -set_Ricatti_State_opt_value$State_obs: The optimal trajectory computed at observation points
  #       -set_Ricatti_State_opt_value$Ricatti: The Riccati equation solution computed at discretization points
  #       -set_Ricatti_State_opt_value$Adjoint: The adjoint equation solution computed at discretization points
  #       -set_Ricatti_State_opt_value$val_RSS: Residual sum of square between the perturbed ODE solutions and observation
  ############################################################
dim_syst = ncol(mat_C)
dim_obs = nrow(mat_C)
dim_control = ncol(mat_U)

nb_time_integ = length(Times_integ);
Yn_fin =pseudo_Y[,nb_time_integ];

mat_Rn  = t(mat_C)%*%mat_C
vect_Hn = -t(mat_C)%*%Yn_fin 
History_R = list()
History_H = list()
History_R = list(mat_Rn) 
History_H = list(vect_Hn) 
mat_Rk  = mat_Rn
vect_Hk = vect_Hn

val_Sn_prof = t(Yn_fin)%*%Yn_fin

for (nk in 1:(nb_time_integ-1)){
  delta_k = Times_integ[nb_time_integ-nk+1] - Times_integ[nb_time_integ-nk] 
  Yk= pseudo_Y[,nb_time_integ-nk]
  
  mat_Ak = mat_A(Times_integ[nb_time_integ-nk],param_struct)
  mat_Ak_delta = mat_Ak*delta_k + diag(dim_syst)
  weight_cur = weight_integ[nb_time_integ-nk]
  mat_G = solve(mat_U+delta_k*t(mat_B)%*%mat_Rk%*%mat_B)
   
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
est_x0 = x0
val_Sn_prof  = val_Sn_prof + t(est_x0)%*%mat_R0%*%est_x0+2*t(vect_H0)%*%est_x0 


State_i = matrix(data = rep(0,dim_syst*nb_time_integ ), nrow = dim_syst, ncol = nb_time_integ )
State_i[,1] =  est_x0

control_i =matrix(data = rep(0,dim_control*(nb_time_integ-1)), nrow = dim_control, ncol = nb_time_integ-1 )

for (nkk in 1:(nb_time_integ-1)){
  delta_k = Times_integ[nkk+1] - Times_integ[nkk] 
  mat_Ak = mat_A(Times_integ[nkk],param_struct)
  mat_Ak_delta = mat_Ak*delta_k + diag(dim_syst)
  
  mat_Rk_p1 = History_R[[nkk+1]]
  vect_Hk_p1 = History_H[[nkk+1]]
  mat_G = solve(mat_U+delta_k*t(mat_B)%*%mat_Rk_p1%*%mat_B)
  control_ukk = -mat_G%*%t(mat_B)%*%(mat_Rk_p1%*%mat_Ak_delta%*%State_i[,nkk]+vect_Hk_p1)
  control_i[,nkk] = control_ukk
  State_i[,nkk+1]  = mat_Ak_delta%*%State_i[,nkk]+delta_k*mat_B%*%control_ukk

}
indice_obs  =weight_integ>0      
State_i_obs_time = State_i[,indice_obs]

ensemble_est = list(Val_cost_prof = val_Sn_prof ,control_opt =control_i, Ricatti = History_R,Adjoint =  History_H,State_opt =State_i,State_obs = State_i_obs_time)

return(ensemble_est)
}