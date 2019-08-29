function [loc_max_full, loc_max_second] = loc_max_omega(theta0)


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

[X_A, X_D, Z_A, apply, f, risk_k, L, T, P, b, case_l, p_cap, share_repay, cost_k, M_m, b_m] = unpack_data_ind(A_DATA);
[X_A_m, X_D_m, Z_A_m, apply_m, f_m, risk_k_m, L_m, T_m, P_m, b_m, p_cap_m, share_repay_m, M_m, cost_k_m] = unpack_data(A_DATA);

%--------------------------------------------------------------------------
% FULL CHOICE SET
%--------------------------------------------------------------------------

% Assign delta to banks
delta_m = tau_m + M_gamma_m;

% Exp max
loc_m = delta_m; %- sigma_omega_m * 0.5772156649015328606065120;
num = exp(loc_m ./ sigma_omega_m);

den = kron(accumarray(A_ij(:,1), num),ones(J,1));
loc_max_full = sigma_omega_m .* (log(den)); %- sigma_omega * 0.5772156649015328606065120  ; % + sigma_omega * 0.5772156649015328606065120

%--------------------------------------------------------------------------
% CHOICE SET EXCLUDING TOP CHOICE
%--------------------------------------------------------------------------

% Assign delta to banks
delta_m_second = tau_m + M_gamma_m;

% Exp max
loc_m_second = delta_m_second; %- sigma_omega_m * 0.5772156649015328606065120;
num_second = exp(loc_m_second ./ sigma_omega_m);

num_second(A_ij(:,2) == b_m) = 0;
den_second = kron(accumarray(A_ij(:,1), num_second),ones(J,1));
loc_max_second = sigma_omega_m .* (log(den_second)); % - sigma_omega * 0.5772156649015328606065120  ; % + sigma_omega * 0.5772156649015328606065120  

end