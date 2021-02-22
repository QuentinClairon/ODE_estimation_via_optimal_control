Repressilator_sim <- function(Time, State, Pars)
{
  v_max= Pars[1:3]
  k_max = Pars[4:6]
  k_g = Pars[7:9]
  k_mrna = Pars[10:12]
  k_p = Pars[13:15]
  n = Pars[16]
  
  r1 = State[1]
  r2 = State[2]
  r3 = State[3]
  
  p1 = State[4]
  p2 = State[5]
  p3 = State[6]
  
  dr1 = v_max[1]*k_max[1]^n/(p2^n+k_max[1]^n) - k_g[1]*r1
  dr2 = v_max[2]*k_max[2]^n/(p3^n+k_max[2]^n) - k_g[2]*r2
  dr3 = v_max[3]*k_max[3]^n/(p1^n+k_max[3]^n) - k_g[3]*r3
  
  dp1 = k_mrna[1]*r1 - k_p[1]*p1
  dp2 = k_mrna[2]*r2 - k_p[2]*p2
  dp3 = k_mrna[3]*r3 - k_p[3]*p3
  
  return(list(c(dr1,dr2,dr3,dp1,dp2,dp3)))
  
}