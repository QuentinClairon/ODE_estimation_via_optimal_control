Microbiota_matA_log_transf_rest_multi <- function(t,State,param,Pars_set)
{
 
  theta =param

  
  mat_AF = matrix(0,8,8)
  
  vect_g_rec = Pars_set[1:7]
  vect_s_rec = Pars_set[8:14]
  M_sp_x9 = Pars_set[15:21]
  
  th_21 =  Pars_set[22]
  th_12 =  Pars_set[23]
  th_43 =  Pars_set[24]
  th_53 =  Pars_set[25]
  th_34 =  Pars_set[26]
  th_54 =  Pars_set[27]
  th_15 =  Pars_set[28]
  th_16 =  Pars_set[29]
  th_26 =  Pars_set[30]
  th_46 =  Pars_set[31]
  th_114 = Pars_set[32]
  
  M_x1_col = c(theta[1],th_21,theta[2:4])
  M_x2_col = c(th_12,  theta[5:8])
  M_x3_col = c(theta[9:11], th_43, th_53)
  M_x4_col = c(theta[12:13], th_34, theta[14], th_54)
  M_x5_col = c(th_15,theta[15:18])
  M_x6_col = c(th_16, th_26, theta[19], th_46, theta[20])
  M_x7_col = theta[21:25]
  
  M_sp_x1_x5 = cbind(M_x1_col, M_x2_col, M_x3_col, M_x4_col, M_x5_col, M_x6_col, M_x7_col)
  M_sp_x11 = c(theta[26:28], th_114, theta[29:31])
  
  Matrix_M_rec = rbind(M_sp_x1_x5,M_sp_x9,M_sp_x11)
  vect_1 = matrix(1,7,1)
  
  Mat_exp = diag((exp(State[1:7])-1)/State[1:7])
  
  pert = 0
  if(t==0){
    pert=1
  }
 
  mat_AF[1:7,1:7] = Matrix_M_rec%*%Mat_exp
  mat_AF[8,1:7] = vect_g_rec +vect_s_rec*pert +Matrix_M_rec%*%vect_1 
 
  return( mat_AF)
}
