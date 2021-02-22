Apinene_matA<- function(t,param)
{
  
  res_t = matrix (rep(0, 25), 5, 5)
  #dx1 = -(Pars[1]+Pars[2])*State[1]
  res_t[1,1] =-(param[1]+param[2]) 
  
  #dx2 = Pars[1]*State[1]
  res_t[2,1] =param[1]
  
  #dx3 = Pars[2]*State[1]- (Pars[3]+Pars[4])*State[3]+Pars[5]*State[5]
  res_t[3,1] =  param[2]
  res_t[3,3] =  - (param[3]+param[4])
  res_t[3,5] =   param[5]
  
  #dx4 = Pars[3]*State[3]
  res_t[4,3] = param[3]
  
  #dx5 = Pars[4]*State[3]-Pars[5]*State[5]
  res_t[5,3] = param[4]
  res_t[5,5] = -param[5]
  
  return(res_t)
  
}