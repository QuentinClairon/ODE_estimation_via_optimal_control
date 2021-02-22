fhn_fun_a_matA<- function(t,State,param){
  
b = param[1]
c = param[2]

res_t = matrix(0,4,4)
res_t[1,1] = c*(1-(State[1]^2/3))
res_t[1,2]  = c;

res_t[2,1]= -1/c;
res_t[2,2]= -b/c;
res_t[2,3] = 1/c;

res_t[3,4]  = 1;

return(res_t)

}