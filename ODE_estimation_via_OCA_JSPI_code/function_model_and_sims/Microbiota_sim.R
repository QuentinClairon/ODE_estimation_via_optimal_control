Microbiota_sim <- function(Time, State, Pars)
{
  Matrix_M_rec = matrix(Pars[1:121], 11, 11)
  vect_g = Pars[122:132]
  vect_s = Pars[133:143]
  
  pert = 0
  if(Time==0){
    pert=1
  }
  
  dx = (vect_g+Matrix_M_rec%*%State +vect_s*pert)*State
  
  return(list(c(dx)))
  
}