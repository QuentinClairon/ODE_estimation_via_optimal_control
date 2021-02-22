fhn_stoch_pert_sim <- function(deb,h,fin,Pars,fun_a,y0,sd_stoch_pert){
  
  Times_obs = seq(deb,fin, by=h) 
  nb_obs = length(Times_obs)
  
  dim_syst =  2
  
  b= Pars[1]
  c= Pars[2]
  
  Res_sim = matrix(rep(0,dim_syst*nb_obs),dim_syst,nb_obs)
  Res_sim[,1] = y0
  for (i in 1 : (nb_obs-1)){
    delta_i = Times_obs[i+1]-Times_obs[i]
    State = Res_sim[,i]

    V = State[1]
    R = State[2]
    dV = c*(V-V^3/3+R)
    dR = -(V-fun_a(Times_obs[i])+b*R)/c
    
    part_det_i = matrix(c(0,0),dim_syst,1)
    part_det_i[1] =  dV
    part_det_i[2] =  dR

    Res_sim[,i+1] = State +delta_i*part_det_i+ delta_i*sd_stoch_pert*rnorm(2,0)*State
  }
  return(list(Times_obs,Res_sim))
}