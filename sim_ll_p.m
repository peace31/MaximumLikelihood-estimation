function [LL,G,gceq] = sim_ll_p(theta0)

%--------------------------------------------------------------------------
% GLOBALS
%--------------------------------------------------------------------------
ceq=[];
gceq=[];
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
% Construct Likelihood
%--------------------------------------------------------------------------

% Inputs
[cp, cp_choice] = choice_prob(theta0);
[pdf_full, pdf_second, cdf_full, cdf_second, cdf_choice,g_pdf_full,g_pdf_full_sigma,g_pdf_second,g_pdf_second_sigma,g_cdf_full,g_cdf_full_sigma,g_cdf_second,g_cdf_second_sigma,g_cdf_choice,g_cdf_choice_tao,g_cdf_choice_gamma,g_cdf_choice_sigma,g_cdf_full1,g_cdf_full1_sigma] = order_stat(theta0);

% Case 1: Unconstrained
L_case_1 = 1e-99 + pdf_second + (repmat(cp_choice,1,N_draws_r_cost) - 1) .* pdf_full;
g_eta_1=(g_pdf_second+(repmat(cp_choice,1,N_draws_r_cost) - 1) .*g_pdf_full)./L_case_1;
g_sigma_1=(g_pdf_second_sigma+(repmat(cp_choice,1,N_draws_r_cost) - 1) .*g_pdf_full_sigma)./L_case_1;
% Case 2: Constrained
const1=(cdf_second + (repmat(cp_choice,1,N_draws_r_cost) - 1) .* cdf_full);
const2=(1 - cdf_choice);
L_case_2 = 1e-99 + const1 .* const2;
g_eta_2=(g_cdf_second + (repmat(cp_choice,1,N_draws_r_cost) - 1) .* g_cdf_full)./(1e-99 + const1)-g_cdf_choice./const2;
g_sigma_2=(g_cdf_second_sigma + (repmat(cp_choice,1,N_draws_r_cost) - 1) .* g_cdf_full_sigma) ./(1e-99 + const1)-g_cdf_choice_sigma./const2;
g_tao_2=- g_cdf_choice_tao./const2;
g_gmama_2=-g_cdf_choice_gamma./const2;
% Case 3: Rejected
L_case_3 = 1e-99 + cdf_full;
g_eta_3=g_cdf_full1;
g_sigma_3=g_cdf_full1_sigma;
% Indicators
I_case_1 = case_l == 1;
I_case_2 = case_l == 2;
I_case_3 = case_l == 3;

% Assign likelihoods to cases
L = zeros (N_app, N_draws_r_cost);
gL_eta = zeros (N_app, N_draws_r_cost);
gL_tao = zeros (N_app, N_draws_r_cost);
gL_gamma = zeros (N_app, N_draws_r_cost);
gL_sigma = zeros (N_app, N_draws_r_cost);
L(I_case_1,:) = L_case_1(I_case_1,:);
gL_eta(I_case_1,:) = g_eta_1(I_case_1,:);
gL_sigma(I_case_1,:) = g_sigma_1(I_case_1,:);
L(I_case_2,:) = L_case_2(I_case_2,:);
gL_eta(I_case_2,:) = g_eta_2(I_case_2,:);
gL_tao(I_case_2,:) = g_tao_2(I_case_2,:);
gL_gamma(I_case_2,:) = g_gmama_2(I_case_2,:);
gL_sigma(I_case_2,:) = g_sigma_2(I_case_2,:);
L(I_case_3,:) = L_case_3(I_case_3,:);
gL_eta(I_case_3,:) = g_eta_3(I_case_3,:);
gL_sigma(I_case_3,:) = g_sigma_3(I_case_3,:);


% Add across draws
LI = (1 / N_draws_r_cost) * sum(L,2);
gLI_eta = (1 / N_draws_r_cost) * sum(gL_eta,2);
gLI_tao = (1 / N_draws_r_cost) * sum(gL_tao,2);
gLI_gamma = (1 / N_draws_r_cost) * sum(gL_gamma,2);
gLI_sigma = (1 / N_draws_r_cost) * sum(gL_sigma,2);

% Add across choice sets
% % L_C_set = L .* P_C_set;
% % LI = sum(L_C_set,2);

% Log likelihood
log_LI = log(LI);

% Fix problems
log_LI(log_LI==-Inf) = -999;
log_LI(isnan(log_LI))= -999;
gLI_eta(gLI_eta==-Inf) = 0;
gLI_eta(isnan(gLI_eta))= 0;
gLI_tao(gLI_tao==-Inf) = 0;
gLI_tao(isnan(gLI_tao))= 0;
gLI_gamma(gLI_gamma==-Inf) = 0;
gLI_gamma(isnan(gLI_gamma))= 0;
gLI_sigma(gLI_sigma==-Inf) = 0;
gLI_sigma(isnan(gLI_sigma))= 0;

% Output
N_apply = sum(I_case_1 + I_case_2 + I_case_3);
LL = -(1 / N_apply) * sum(log_LI,1);
gLL_eta = -(1 / N_apply) * sum(gLI_eta,1);
gLL_tao = -(1 / N_apply) * sum(gLI_tao,1);
gLL_gama = -(1 / N_apply) * sum(gLI_gamma,1);
gLL_sigma = -(1 / N_apply) * sum(gLI_sigma,1);
%LL = - sum(log(LI));
G=zeros(1,length(theta0));
G(1                                :K_K)=gLL_eta;
G(K_K+1                            :K_K+J*K_K)=gLL_tao;
G(K_K+J*K_K+1                      :K_K+J*K_K+K_M*K_K)=gLL_gama;
G(K_K+J*K_K+K_M*K_K+1              :K_K+J*K_K+K_M*K_K+K_K)=gLL_sigma;
G=G';
%--------------------------------------------------------------------------
% Output parameter and likelihood function values
%--------------------------------------------------------------------------

fid = fopen('likefull_p.txt','wt');
fprintf(fid,'LL = %6.3f\n',LL);
fprintf(fid,'pr(unconstrained|apply) = %6.3f\n',mean(LI(I_case_1 == 1, :)));
fprintf(fid,'pr(constrained|apply) = %6.3f\n',  mean(LI(I_case_2 == 1, :)));
fprintf(fid,'pr(reject|apply) = %6.3f\n',       mean(LI(I_case_3 == 1, :)));
fprintf(fid,'%6.3f\n',theta0);
fclose(fid);

end


