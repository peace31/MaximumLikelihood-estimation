function [pdf_full, pdf_second, cdf_full, cdf_second, cdf_choice,g_pdf_full,g_pdf_full_sigma,g_pdf_second,g_pdf_second_sigma,g_cdf_full,g_cdf_full_sigma,g_cdf_second,g_cdf_second_sigma,g_cdf_choice,g_cdf_choice_tao,g_cdf_choice_gamma,g_cdf_choice_sigma,g_cdf_full1,g_cdf_full1_sigma] = order_stat(theta0)

%--------------------------------------------------------------------------
% GLOBALS
%--------------------------------------------------------------------------

global N_app J K_A K_D K_M K_R K_K K_Z K_C K_C_set          % Dimensions
global A_DATA DATA C_set A_ij A_i A_j I_ij id_j zeta        % Data and indices
global N_draws_repay N_draws_app                            % Number of draws
global P_app_draws prob_app_draws draws_omega_app           % Application draws
global N_draws_r_cost draws_e_repay share_repay_draws dr    % Supply draws
global P_rate L_rate T_rate                                 % Payments to interest rates
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
% Probabilities of order statistics
%--------------------------------------------------------------------------

% Location parameters
[loc_max_full, loc_max_second] = loc_max_omega(theta0);

% Assign delta to banks
delta_m = tau_m + M_gamma_m;

% Borrower-specific cost
%r_m = repmat(X_R_m * beta_r,1,N_draws_r_cost) + kron(sigma_r_m .* draws_e_r_cost,ones(J,1));
r_m = zeros(N_app*J,N_draws_r_cost);

% PV factor 
pv_factor_repay_draws = mean((1 ./ repmat(dr,N_app*J,1)) .* (1 - exp(-repmat(dr,N_app*J,1) .* kron(share_repay_draws, ones(J,1)) .* T_m * 12)),2);
pv_factor_full_draws  = (1 ./ repmat(dr,N_app*J,1)) .* (1 - exp(-repmat(dr,N_app*J,1) .*                                            T_m * 12));
pv_factor_ratio_draws_inv = pv_factor_repay_draws ./ pv_factor_full_draws;

% Construct PDF for full choice set
pdf_full_m =  (1 ./ sigma_omega_m) ...
            .* exp(     -((repmat((-(pv_factor_ratio_draws_inv .* P_m - eta_f_m .* f_m) ./ L_m),1,N_draws_r_cost) + r_m - repmat(loc_max_full,1,N_draws_r_cost)) ./ (sigma_omega_m))) ...
            .* exp(-exp(-((repmat((-(pv_factor_ratio_draws_inv .* P_m - eta_f_m .* f_m) ./ L_m),1,N_draws_r_cost) + r_m - repmat(loc_max_full,1,N_draws_r_cost)) ./ (sigma_omega_m))));                
cc=-(repmat((-(pv_factor_ratio_draws_inv .* P_m - eta_f_m .* f_m) ./ L_m),1,N_draws_r_cost) + r_m - repmat(loc_max_full,1,N_draws_r_cost));
g_full_c=-f_m./sigma_omega_m./L_m.*(1-exp(-((repmat((-(pv_factor_ratio_draws_inv .* P_m - eta_f_m .* f_m) ./ L_m),1,N_draws_r_cost) + r_m - repmat(loc_max_full,1,N_draws_r_cost)) ./ (sigma_omega_m))));
g_full_sigma_c=-1./sigma_omega_m.^2.*(sigma_omega_m+cc-cc.*exp(-((repmat((-(pv_factor_ratio_draws_inv .* P_m - eta_f_m .* f_m) ./ L_m),1,N_draws_r_cost) + r_m - repmat(loc_max_full,1,N_draws_r_cost)) ./ (sigma_omega_m))));
g_pdf_full_m=g_full_c.*pdf_full_m;
g_full_sigma_m=g_full_sigma_c.*pdf_full_m;
g_pdf_full=g_pdf_full_m(A_ij(:,2) == b_m(:,1),:);
pdf_full = pdf_full_m(A_ij(:,2) == b_m(:,1),:);
g_pdf_full_sigma=g_full_sigma_m(A_ij(:,2) == b_m(:,1),:);

% Construct PDF excluding top choice
pdf_second_m = (1 ./ sigma_omega_m) ...
            .* exp(     -((repmat((-(pv_factor_ratio_draws_inv .* P_m - eta_f_m .* f_m) ./ L_m),1,N_draws_r_cost) + r_m - repmat(loc_max_second,1,N_draws_r_cost)) ./ (sigma_omega_m))) ...
            .* exp(-exp(-((repmat((-(pv_factor_ratio_draws_inv .* P_m - eta_f_m .* f_m) ./ L_m),1,N_draws_r_cost) + r_m - repmat(loc_max_second,1,N_draws_r_cost)) ./ (sigma_omega_m))));
cc1=-(repmat((-(pv_factor_ratio_draws_inv .* P_m - eta_f_m .* f_m) ./ L_m),1,N_draws_r_cost) + r_m - repmat(loc_max_second,1,N_draws_r_cost));
g_second_c=-f_m./sigma_omega_m./L_m.*(1-exp(-((repmat((-(pv_factor_ratio_draws_inv .* P_m - eta_f_m .* f_m) ./ L_m),1,N_draws_r_cost) + r_m - repmat(loc_max_second,1,N_draws_r_cost)) ./ (sigma_omega_m))));
g_second_simga_c=-1./sigma_omega_m.^2.*(sigma_omega_m+cc1-cc1.*exp(-((repmat((-(pv_factor_ratio_draws_inv .* P_m - eta_f_m .* f_m) ./ L_m),1,N_draws_r_cost) + r_m - repmat(loc_max_second,1,N_draws_r_cost)) ./ (sigma_omega_m))));
g_pdf_second_m=g_second_c.*pdf_second_m;
g_pdf_second_sigma_m=g_second_simga_c.*pdf_second_m;
g_pdf_second=g_pdf_second_m(A_ij(:,2) == b_m(:,1),:);
g_pdf_second_sigma=g_pdf_second_sigma_m(A_ij(:,2) == b_m(:,1),:);
pdf_second = pdf_second_m(A_ij(:,2) == b_m(:,1),:);

% Construct CDF for full choice set
cdf_full_m = exp(-exp(-((repmat((-(pv_factor_ratio_draws_inv .* P_m - eta_f_m .* f_m) ./ L_m),1,N_draws_r_cost) + r_m - repmat(loc_max_full,1,N_draws_r_cost)) ./ (sigma_omega_m))));
cc2=-(repmat((-(pv_factor_ratio_draws_inv .* P_m - eta_f_m .* f_m) ./ L_m),1,N_draws_r_cost) + r_m - repmat(loc_max_full,1,N_draws_r_cost));
g_cdf_full_m=f_m./sigma_omega_m./L_m.*cdf_full_m.*exp(-((repmat((-(pv_factor_ratio_draws_inv .* P_m - eta_f_m .* f_m) ./ L_m),1,N_draws_r_cost) + r_m - repmat(loc_max_full,1,N_draws_r_cost)) ./ (sigma_omega_m)));
g_cdf_full_sigma_m=1./sigma_omega_m.^2.*cc2.*cdf_full_m.*exp(-((repmat((-(pv_factor_ratio_draws_inv .* P_m - eta_f_m .* f_m) ./ L_m),1,N_draws_r_cost) + r_m - repmat(loc_max_full,1,N_draws_r_cost)) ./ (sigma_omega_m)));
g_cdf_full1_m=f_m./sigma_omega_m./L_m.*exp(-((repmat((-(pv_factor_ratio_draws_inv .* P_m - eta_f_m .* f_m) ./ L_m),1,N_draws_r_cost) + r_m - repmat(loc_max_full,1,N_draws_r_cost)) ./ (sigma_omega_m)));
g_cdf_full1_sigma_m=1./sigma_omega_m.^2.*cc2.*exp(-((repmat((-(pv_factor_ratio_draws_inv .* P_m - eta_f_m .* f_m) ./ L_m),1,N_draws_r_cost) + r_m - repmat(loc_max_full,1,N_draws_r_cost)) ./ (sigma_omega_m)));

cdf_full = cdf_full_m(A_ij(:,2) == b_m(:,1),:);
g_cdf_full=g_cdf_full_m(A_ij(:,2) == b_m(:,1),:);
g_cdf_full_sigma=g_cdf_full_sigma_m(A_ij(:,2) == b_m(:,1),:);
g_cdf_full1=g_cdf_full1_m(A_ij(:,2) == b_m(:,1),:);
g_cdf_full1_sigma=g_cdf_full1_sigma_m(A_ij(:,2) == b_m(:,1),:);

% Construct CDF excluding top choice
cdf_second_m = exp(-exp(-(repmat((-(pv_factor_ratio_draws_inv .* P_m - eta_f_m .* f_m) ./ L_m),1,N_draws_r_cost) + r_m - repmat(loc_max_second,1,N_draws_r_cost)) ./ (sigma_omega_m)));
cc3=-(repmat((-(pv_factor_ratio_draws_inv .* P_m - eta_f_m .* f_m) ./ L_m),1,N_draws_r_cost) + r_m - repmat(loc_max_second,1,N_draws_r_cost));
g_cdf_second_m=f_m./sigma_omega_m./L_m.*cdf_second_m.*exp(-(repmat((-(pv_factor_ratio_draws_inv .* P_m - eta_f_m .* f_m) ./ L_m),1,N_draws_r_cost) + r_m - repmat(loc_max_second,1,N_draws_r_cost)) ./ (sigma_omega_m));
g_cdf_second_sigma_m=1./sigma_omega_m.^2.*cc3.*cdf_second_m.*exp(-(repmat((-(pv_factor_ratio_draws_inv .* P_m - eta_f_m .* f_m) ./ L_m),1,N_draws_r_cost) + r_m - repmat(loc_max_second,1,N_draws_r_cost)) ./ (sigma_omega_m));
cdf_second = cdf_second_m(A_ij(:,2) == b_m(:,1),:);
g_cdf_second=g_cdf_second_m(A_ij(:,2) == b_m(:,1),:);
g_cdf_second_sigma=g_cdf_second_sigma_m(A_ij(:,2) == b_m(:,1),:);

% Construct CDF conditional on top choice
cdf_choice_m = exp(-exp(-(repmat((-(pv_factor_ratio_draws_inv .* P_m - eta_f_m .* f_m) ./ L_m),1,N_draws_r_cost) + r_m - repmat(delta_m,1,N_draws_r_cost)) ./ (sigma_omega_m)));
cc4=-(repmat((-(pv_factor_ratio_draws_inv .* P_m - eta_f_m .* f_m) ./ L_m),1,N_draws_r_cost) + r_m - repmat(delta_m,1,N_draws_r_cost));
g_cdf_choice_m=f_m./sigma_omega_m./L_m.*cdf_second_m.*exp(-(repmat((-(pv_factor_ratio_draws_inv .* P_m - eta_f_m .* f_m) ./ L_m),1,N_draws_r_cost) + r_m - repmat(delta_m,1,N_draws_r_cost)) ./ (sigma_omega_m));
g_cdf_choice_sigma_m=1./sigma_omega_m.^2.*cc4.*cdf_second_m.*exp(-(repmat((-(pv_factor_ratio_draws_inv .* P_m - eta_f_m .* f_m) ./ L_m),1,N_draws_r_cost) + r_m - repmat(delta_m,1,N_draws_r_cost)) ./ (sigma_omega_m));
g_cdf_choice_m_tao=-1./sigma_omega_m.*cdf_second_m.*exp(-(repmat((-(pv_factor_ratio_draws_inv .* P_m - eta_f_m .* f_m) ./ L_m),1,N_draws_r_cost) + r_m - repmat(delta_m,1,N_draws_r_cost)) ./ (sigma_omega_m));
g_cdf_choice_m_gamma=-dr./sigma_omega_m.*cdf_second_m.*exp(-(repmat((-(pv_factor_ratio_draws_inv .* P_m - eta_f_m .* f_m) ./ L_m),1,N_draws_r_cost) + r_m - repmat(delta_m,1,N_draws_r_cost)) ./ (sigma_omega_m));
cdf_choice = cdf_choice_m(A_ij(:,2) == b_m(:,1),:);
g_cdf_choice=g_cdf_choice_m(A_ij(:,2) == b_m(:,1),:);
g_cdf_choice_tao=g_cdf_choice_m_tao(A_ij(:,2) == b_m(:,1),:);
g_cdf_choice_gamma=g_cdf_choice_m_gamma(A_ij(:,2) == b_m(:,1),:);
g_cdf_choice_sigma=g_cdf_choice_sigma_m(A_ij(:,2) == b_m(:,1),:);
end

