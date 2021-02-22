Microbiota_sim_log_transf_rest <- function(Time, State, Pars)
{
  Matrix_M_rec = matrix(Pars[1:49], 7, 7)
  vect_g_rec = Pars[50:56]
  vect_s_rec = Pars[57:63]
  
  pert = 0
  if(Time==0){
    pert=1
  }
  
  dlogx = (vect_g_rec+Matrix_M_rec%*%exp(State) +vect_s_rec*pert)
  
  return(list(c(dlogx)))
  
}