% % function [eta_f, tau, gamma, sigma_omega, lambda] = unpack_parm(theta0)
function [eta_f, tau, gamma, sigma_omega, alpha_a, kappa, alpha_d, sigma_d, alpha_p] = unpack_parm(theta0)

%--------------------------------------------------------------------------
% GLOBALS
%--------------------------------------------------------------------------

global N_app J K_A K_D K_M K_R K_K K_Z K_C K_C_set
global A_DATA DATA C_set A_ij A_i A_j I_ij id_j 
global e N_draws_r N_draws_app e_draws
global P_rate L_rate T_rate

%--------------------------------------------------------------------------
% EXTRACT PARAMETERS
%--------------------------------------------------------------------------

% Supply
eta_f =         theta0(1                                :K_K);
tau =           theta0(K_K+1                            :K_K+J*K_K);
gamma =         theta0(K_K+J*K_K+1                      :K_K+J*K_K+K_M*K_K);
sigma_omega =   theta0(K_K+J*K_K+K_M*K_K+1              :K_K+J*K_K+K_M*K_K+K_K);

% Demand               
alpha_a =       theta0(K_K+J*K_K+K_M*K_K+K_K+1          :K_K+J*K_K+K_M*K_K+K_K+K_A);
kappa =         theta0(K_K+J*K_K+K_M*K_K+K_K+K_A+1      :K_K+J*K_K+K_M*K_K+K_K+K_A+K_Z);
alpha_d =       theta0(K_K+J*K_K+K_M*K_K+K_K+K_A+K_Z+1  :K_K+J*K_K+K_M*K_K+K_K+K_A+K_Z+K_D);
sigma_d =       theta0(K_K+J*K_K+K_M*K_K+K_K+K_A+K_Z+K_D+1  :K_K+J*K_K+K_M*K_K+K_K+K_A+K_Z+K_D+1);
alpha_p =       theta0(K_K+J*K_K+K_M*K_K+K_K+K_A+K_Z+K_D+1+1  :K_K+J*K_K+K_M*K_K+K_K+K_A+K_Z+K_D+1+1);

end