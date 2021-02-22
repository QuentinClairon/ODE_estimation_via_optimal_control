The code contained in this folder is here to reproduce the simulated section results given by the optimal control based methods presented in the article
 "A regularization method for the parameter estimation problem in ordinary differential equations via discrete optimal control theory".

For the code to work you simply need to set this folder as the working directory. It is composed of several scripts and two folders. 

Script:
The scripts simulate data set according to the models and experimental designs presented in "Section: 5 Experiments"
and then proceed to estimation by using the optimal control based methods with the same hyperparameters as in the article:
1) Scripts "script_test_iter_oca_apinene" and "script_test_iter_oca_apinene_misspe" reproduce experimental designs and optimal control based estimation procedures of section 5.1
2) Scripts "script_test_iter_oca_repressilator" and "script_test_iter_oca_repressilator_misspe" reproduce experimental designs and optimal control based estimation procedures of section 5.2
3) Scripts "script_test_iter_oca_fhn_fun" and "script_test_iter_oca_fhn_fun_misspe" reproduce experimental designs and optimal control based estimation procedures of section 5.3
4) Scripts "script_test_iter_oca_microbiota" and "script_test_iter_oca_microbiota_misspe" reproduce experimental designs and optimal control based estimation procedures of section 6.1

Folders:
1) The folder "function_model_and_sims" contains all the needed model specific functions to generate the observations and estimate parameters from the previous scripts.

2) The folder "function_util_estimation" contains the generic functions for estimating parameters via our optimal control based approaches for any ODE models. 
They are commented in order to be reusable for other practitioner models 