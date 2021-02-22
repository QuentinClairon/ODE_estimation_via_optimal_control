Repressilator_matA<- function(t,State,param)
{
 
  k_max_31_set = 40
  k_mrna_set = c(5,6,7)
  #k_p_set = c(1,2,3)
  n_set = 3
  
  v_max= param[1:3]
  k_max = c(param[4:5],k_max_31_set)
  k_g = param[6:8]
  k_mrna =  k_mrna_set
  k_p = param[9:11]
  sigm_shape = n_set
  
  res_t = matrix (rep(0,7^2), 7, 7)
  
  
  #dy(1) = v1_max*km1^sigm_shape./(y(5)^sigm_shape+km1^sigm_shape) - k_deg_m1*y(1); 
  res_t[1,1] =  -k_g[1]
  res_t[1,7] =  v_max[1]*k_max[1]^sigm_shape/(State[5]^sigm_shape+k_max[1]^sigm_shape)
  
  #dy(2) = v2_max*km2^sigm_shape./(y(6)^sigm_shape+km2^sigm_shape) - k_deg_m2*y(2);
  res_t[2,2] =  -k_g[2]
  res_t[2,7] =  v_max[2]*k_max[2]^sigm_shape/(State[6]^sigm_shape+k_max[2]^sigm_shape)
  
  #dy(3) = v3_max*km3^sigm_shape./(y(4)^sigm_shape+km3^sigm_shape) - k_deg_m3*y(3); 
  res_t[3,3] =  -k_g[3]
  res_t[3,7] =  v_max[3]*k_max[3]^sigm_shape/(State[4]^sigm_shape+k_max[3]^sigm_shape)
  
  
  #dy(4) = k_transc1*y(1) - k_deg_p1*y(4);
  res_t[4,1] =k_mrna[1]
  res_t[4,4] =-k_p[1]
  
  #dy(5) = k_transc2*y(2) - k_deg_p2*y(5);
  res_t[5,2] =k_mrna[2]
  res_t[5,5] =-k_p[2]
  
  #%dy(6) = k_transc3*y(3) - k_deg_p3*y(6);
  res_t[6,3] =k_mrna[3]
  res_t[6,6] =-k_p[3]
  
  #print(res_t)
  return(res_t)
  
}