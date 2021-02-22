rep_stoch_pert_sim <- function(deb,h,fin,Pars,y0,sd_stoch_pert){
  
  Times_obs = seq(deb,fin, by=h) 
  nb_obs = length(Times_obs)
  
  dim_syst =  6
  
  v_max= Pars[1:3]
  k_max = Pars[4:6]
  k_g = Pars[7:9]
  k_mrna = Pars[10:12]
  k_p = Pars[13:15]
  n = Pars[16]
  
  Res_sim = matrix(rep(0,dim_syst*nb_obs),dim_syst,nb_obs)
  Res_sim[,1] = y0
  for (i in 1 : (nb_obs-1)){
    delta_i = Times_obs[i+1]-Times_obs[i]
    State = Res_sim[,i]
    
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
    
    part_det_i = matrix(0,dim_syst,1)
    part_det_i[1] =  dr1
    part_det_i[2] =  dr2
    part_det_i[3] =  dr3
    
    part_det_i[4] =  dp1
    part_det_i[5] =  dp2
    part_det_i[6] =  dp3
    
    Res_sim[,i+1] = State +delta_i*part_det_i+ delta_i*sd_stoch_pert*rnorm(6,0)*State
  }
  return(list(Times_obs,Res_sim))
}