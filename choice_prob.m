function [cp, cp_choice] = choice_prob(theta0)


%--------------------------------------------------------------------------
% GLOBALS
%--------------------------------------------------------------------------

global N_app J K_A K_D K_M K_R K_K K_Z K_C K_C_set
global A_DATA DATA C_set A_ij A_i A_j I_ij id_j 
global N_draws_r_cost N_draws_r_app N_draws_app draws_e_r_cost draws_e_r_app
global P_rate L_rate T_rate
global app_prob_app_fit N_app_app
global d_knitro

%--------------------------------------------------------------------------
% Unpack parameters
%--------------------------------------------------------------------------

[eta_f, tau, gamma, sigma_omega, alpha_a, kappa, alpha_d, sigma_d, alpha_p] = unpack_parm(theta0);
[eta_f_m, tau_m, M_gamma_m, sigma_omega_m] = costtype_parm(theta0);

%--------------------------------------------------------------------------
% Unpack data
%--------------------------------------------------------------------------

[X_A_m, X_D_m, Z_A_m, apply_m, f_m, risk_k_m, L_m, T_m, P_m, b_m, p_cap_m, share_repay_m, M_m, cost_k_m] = unpack_data(A_DATA);

%--------------------------------------------------------------------------
% Construct choice probabilities
%--------------------------------------------------------------------------

% Assign delta to banks
delta_m = tau_m + M_gamma_m;

% Choice probability
num = exp(delta_m ./ sigma_omega_m);

den = kron(accumarray(A_ij(:,1), num),ones(J,1));

cp = num ./ den;

cp_choice = cp(A_ij(:,2) == b_m(:,1));

end
