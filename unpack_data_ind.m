%function [X_A, X_D, Z_A, apply, f, risk_k, L, T, P, b, case_l, p_cap, share_repay, pv_factor_full, pv_factor_repay, M_m, b_m] = unpack_data_ind(DATA_IN)
function [X_A, X_D, Z_A, apply, f, risk_k, L, T, P, b, case_l, p_cap, share_repay, cost_k, M_m, b_m] = unpack_data_ind(DATA_IN)

%--------------------------------------------------------------------------
% GLOBALS
%--------------------------------------------------------------------------

global N_app J K_A K_D K_M K_R K_K K_Z K_C K_C_set      % Dimensions
global A_DATA DATA C_set A_ij A_i A_j I_ij id_j zeta    % Data and indices
global N_draws_repay N_draws_app                        % Number of draws
global P_app_draws prob_app_draws draws_omega_app       % Application draws
global N_draws_r_cost share_repay_hat draws_e_repay dr  % Supply draws
global P_rate L_rate T_rate                             % Payments to interest rates
global d_knitro

%--------------------------------------------------------------------------
% EXTRACT DATA
%--------------------------------------------------------------------------

% DATA_IN = A_DATA
% VARIABLES: [I_i, A_id_i, X_A, X_D, X_R, Z_A, apply, f, df_cprob, risk_k, L, P, b, d_cap, p_cap, approved, case_l, default, M, d_choice]

DATA_IND = DATA_IN(DATA_IN(:,end) == 1,:);

X_A =         DATA_IND(:,2                     :2+K_A-1);
X_D =         DATA_IND(:,2+K_A                 :2+K_A+K_D-1);
Z_A =         DATA_IND(:,2+K_A+K_D-1+1         :2+K_A+K_D-1+K_Z);
apply =       DATA_IND(:,2+K_A+K_D-1+K_Z+1     :2+K_A+K_D-1+K_Z+1);
f =           DATA_IND(:,2+K_A+K_D-1+K_Z+2     :2+K_A+K_D-1+K_Z+2);
risk_k =      DATA_IND(:,2+K_A+K_D-1+K_Z+3     :2+K_A+K_D-1+K_Z+3);
L =           DATA_IND(:,2+K_A+K_D-1+K_Z+4     :2+K_A+K_D-1+K_Z+4);
T =           DATA_IND(:,2+K_A+K_D-1+K_Z+5     :2+K_A+K_D-1+K_Z+5);
P =           DATA_IND(:,2+K_A+K_D-1+K_Z+6     :2+K_A+K_D-1+K_Z+6);
b =           DATA_IND(:,2+K_A+K_D-1+K_Z+7     :2+K_A+K_D-1+K_Z+7);
p_cap =       DATA_IND(:,2+K_A+K_D-1+K_Z+9     :2+K_A+K_D-1+K_Z+9);
case_l =      DATA_IND(:,2+K_A+K_D-1+K_Z+11    :2+K_A+K_D-1+K_Z+11);
share_repay = DATA_IND(:,2+K_A+K_D-1+K_Z+12    :2+K_A+K_D-1+K_Z+12);
cost_k =      DATA_IND(:,2+K_A+K_D-1+K_Z+13    :2+K_A+K_D-1+K_Z+13);

b_m =           DATA_IN(:,2+K_A+K_D-1+K_Z+7    :2+K_A+K_D-1+K_Z+7);
M_m =           DATA_IN(:,2+K_A+K_D-1+K_Z+14   :2+K_A+K_D-1+K_Z+13+K_M);

end
