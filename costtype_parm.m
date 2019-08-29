% % function [eta_f_m, tau_m, R_gamma_m, sigma_omega_m, X_C_lambda_m] = risktype_parm(theta0)
function [eta_f_m, tau_m, M_gamma_m, sigma_omega_m] = costtype_parm(theta0)

%--------------------------------------------------------------------------
% GLOBALS
%--------------------------------------------------------------------------

global N_app J K_A K_D K_M K_R K_K K_Z K_C K_C_set
global A_DATA DATA C_set A_ij A_i A_j I_ij id_j 
global N_draws_r_cost N_draws_r_app N_draws_app draws_e_r_cost draws_e_r_app
global P_rate L_rate T_rate
global approved_app_draws P_app_draws
global d_knitro

%--------------------------------------------------------------------------
% EXTRACT PARAMETERS
%--------------------------------------------------------------------------

% Unpack parameters
[eta_f, tau, gamma, sigma_omega, alpha_a, kappa, alpha_d, sigma_d] = unpack_parm(theta0);

% Unpack data
[X_A, X_D, Z_A, apply, f, risk_k, L, T, P, b, case_l, p_cap, share_repay, cost_k, M_m, b_m] = unpack_data_ind(A_DATA);
[X_A_m, X_D_m, Z_A_m, apply_m, f_m, risk_k_m, L_m, T_m, P_m, b_m, p_cap_m, share_repay_m, M_m, cost_k_m] = unpack_data(A_DATA);

% Construct interacted parameters
J_dum = dummies(A_ij(:,2));

eta_f_aux = reshape(eta_f,[1,K_K]);
eta_f_aux_m = repmat(eta_f_aux,N_app*J,1);
idx = sub2ind(size(eta_f_aux_m), 1:size(eta_f_aux_m, 1), cost_k_m');
eta_f_m = eta_f_aux_m(idx');

tau_aux = reshape(tau,[J,K_K]);
tau_aux_m = J_dum * tau_aux;
idx = sub2ind(size(tau_aux_m), 1:size(tau_aux_m, 1), cost_k_m');
tau_m = tau_aux_m(idx');

gamma_aux = reshape(gamma,[K_M,K_K]);
M_gamma_aux_m = M_m * gamma_aux;
idx = sub2ind(size(M_gamma_aux_m), 1:size(M_gamma_aux_m, 1), cost_k_m');
M_gamma_m = M_gamma_aux_m(idx');

sigma_omega_aux = reshape(sigma_omega,[1,K_K]);
sigma_omega_aux_m = repmat(sigma_omega_aux,N_app*J,1);
idx = sub2ind(size(sigma_omega_aux_m), 1:size(sigma_omega_aux_m, 1), cost_k_m');
sigma_omega_m = sigma_omega_aux_m(idx');

% % lambda_aux = reshape(lambda,[K_C,K_K]);
% % X_C_lambda_aux_m = X_C_m * lambda_aux;
% % idx = sub2ind(size(X_C_lambda_aux_m), 1:size(X_C_lambda_aux_m, 1), cost_k_m');
% % X_C_lambda_m = X_C_lambda_aux_m(idx');


end