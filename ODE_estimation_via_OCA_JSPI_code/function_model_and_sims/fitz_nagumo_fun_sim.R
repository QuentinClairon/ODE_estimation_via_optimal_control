fitz_nagumo_fun_sim <- function(Time, State, Pars,fun_a)
{
  b= Pars[1]
  c= Pars[2]
  
  V = State[1]
  R = State[2]
  dV = c*(V-V^3/3+R)
  dR = -(V-fun_a(Time)+b*R)/c
  
  return(list(c(dV,dR)))
  
}