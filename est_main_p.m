%--------------------------------------------------------------------------
% PROJECT: PCB
% OBJECTIVE: Simple simulation for auction model
% DATE STARTED: Jul/15/2017
% DATE LAST MODIFIED: Jun/22/2018
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% PREAMBLE
%--------------------------------------------------------------------------

clear all
close all
%cd('~/---/Code')
%savepath();

% Diary
diary('output_est_main')
diary on
datetime('today')

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
% INPUTS
%--------------------------------------------------------------------------

% Dimensions
K_A = 4;
K_Z = 4;
K_M = 1;
K_D = 4;

% Observations
N_DATA_p = 10000;

% Draws
N_draws_r_cost = 1;
N_draws_repay = 100;

% Discount rate
dr = .05 / 12;

%--------------------------------------------------------------------------
% ESTIMATION OF SUPPLY SIDE
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% IMPORT DATA
%--------------------------------------------------------------------------

A_DATA_RAW = dlmread('A_DATA_RAW.csv');
A_DATA_RAW = A_DATA_RAW(A_DATA_RAW(:,1) <= N_DATA_p ,:);

% Import variables
id_app_id       = A_DATA_RAW(:,1);
X_A             = A_DATA_RAW(:,2              :K_A+1);
X_D             = A_DATA_RAW(:,1+K_A+1        :1+K_A+K_D);
Z_A             = A_DATA_RAW(:,1+K_A+K_D+1    :1+K_A+K_D+K_Z);
apply           = A_DATA_RAW(:,1+K_A+K_D+K_Z+1);
f               = A_DATA_RAW(:,1+K_A+K_D+K_Z+1+1);
risk_k          = A_DATA_RAW(:,1+K_A+K_D+K_Z+1+1+1);
L               = A_DATA_RAW(:,1+K_A+K_D+K_Z+1+1+1+1);
T               = A_DATA_RAW(:,1+K_A+K_D+K_Z+1+1+1+1+1);
P               = A_DATA_RAW(:,1+K_A+K_D+K_Z+1+1+1+1+1+1);
b               = A_DATA_RAW(:,1+K_A+K_D+K_Z+1+1+1+1+1+1+1);
d_cap           = A_DATA_RAW(:,1+K_A+K_D+K_Z+1+1+1+1+1+1+1+1);
p_cap           = A_DATA_RAW(:,1+K_A+K_D+K_Z+1+1+1+1+1+1+1+1+1);
approved        = A_DATA_RAW(:,1+K_A+K_D+K_Z+1+1+1+1+1+1+1+1+1+1);
case_l          = A_DATA_RAW(:,1+K_A+K_D+K_Z+1+1+1+1+1+1+1+1+1+1+1);
share_repay     = A_DATA_RAW(:,1+K_A+K_D+K_Z+1+1+1+1+1+1+1+1+1+1+1+1);
cost_k          = A_DATA_RAW(:,1+K_A+K_D+K_Z+1+1+1+1+1+1+1+1+1+1+1+1+1);
M               = A_DATA_RAW(:,1+K_A+K_D+K_Z+1+1+1+1+1+1+1+1+1+1+1+1+1+1:1+K_A+K_D+K_Z+1+1+1+1+1+1+1+1+1+1+1+1+1+K_M);
d_choice        = A_DATA_RAW(:,1+K_A+K_D+K_Z+1+1+1+1+1+1+1+1+1+1+1+1+1+K_M+1);

%--------------------------------------------------------------------------
% PREPARE FOR ESTIMATION
%--------------------------------------------------------------------------

% Borrowers and banks
N_app = length(unique(id_app_id));  % Number of borrowers
J = length(unique(b));              % Number of banks
K_K = length(unique(cost_k));       % Number of borrower cost types

% Indices
A_j = repmat((1:J)',N_app,1);
A_i = kron((1:N_app)', ones(J,1));
A_ij = [A_i,A_j];

% Lists
A_id_i = unique(A_i);
A_id_j = unique(A_j);

% Halton draws
q_app_repay = qrandstream('halton',N_draws_repay,'Skip',1e3,'Leap',1e2);
draws_e_repay_u = qrand(q_app_repay,N_app);
draws_e_repay = norminv(draws_e_repay_u);

% Data for estimation
A_DATA = [A_i, X_A, X_D, Z_A, apply, f, risk_k, L, T, P, b, d_cap, p_cap, approved, case_l, share_repay, cost_k, M, d_choice];

% Summary of data
display 'applications'
tabulate(apply)
display 'approvals|apply'
tabulate(approved(apply==1))
display 'market shares|apply'
tabulate(b(apply==1))
display 'case|apply'
tabulate(case_l(apply==1))

%--------------------------------------------------------------------------
% ESTIMATION
%--------------------------------------------------------------------------

% Unpack data
[X_A, X_D, Z_A, apply, f, risk_k, L, T, P, b, case_l, p_cap, share_repay, cost_k, M_m, b_m] = unpack_data_ind(A_DATA);
[X_A_m, X_D_m, Z_A_m, apply_m, f_m, risk_k_m, L_m, T_m, P_m, b_m, p_cap_m, share_repay_m, M_m, cost_k_m] = unpack_data(A_DATA);

% Unpack parameters
theta_true = dlmread('theta_true.csv');
[eta_f0, tau0, gamma0, sigma_omega0, alpha_a0, kappa0, alpha_d0, sigma_d0] = unpack_parm(theta_true);

% Generate share_repay using default estimates
share_repay_draws = exp(X_D * alpha_d0 + sigma_d0 * draws_e_repay);
share_repay_draws(share_repay_draws > 1) = 1;
share_repay_draws(share_repay_draws < 0) = 0;

% Initial guess
theta0 = theta_true
%theta0 = (1 + .3 * (rand(length(theta_true),1)-.5)) .* theta_true; 
theta_L = [-10 * ones(K_K,1); -10 * ones(J*K_K,1); -10 * ones(K_M*K_K,1);  0 * ones(K_K,1); -10 * ones(K_A,1); -10 * ones(K_Z,1); -10 * ones(K_D,1); 0 * ones(1,1); -10 * ones(1,1)];
theta_H = [+10 * ones(K_K,1); +10 * ones(J*K_K,1); +10 * ones(K_M*K_K,1); 10 * ones(K_K,1); +10 * ones(K_A,1); +10 * ones(K_Z,1);  10 * ones(K_D,1); 10 * ones(1,1); +10 * ones(1,1)];

% Estimation
options = optimset('Display','iter-detailed','TolFun',1E-10,'TolCon',1E-24,'TolX',1E-24,...
                    'GradObj', 'off', 'DerivativeCheck', 'on',...
                    'Algorithm', 'interior-point', 'MaxFunEvals', 10000);  
options1 = optimoptions('fmincon','Display','iter-detailed','CheckGradients', true,'SpecifyObjectiveGradient',true,'ConstraintTolerance',1E-24,'FunctionTolerance',1E-10,'MaxFunctionEvaluations',10000,'Algorithm', 'interior-point');
[x_p, fval_p, exitflag_p, output_p, lambda_p, grad_p, hessian_p] = ...
                    fmincon(@sim_ll_p,theta0,[],[],[],[],theta_L,theta_H,[],options1);
                
% Compile results            
N_apply = sum(apply);
se_p = sqrt((1 / N_apply) * diag(inv(hessian_p)));
x_p(K_K+J*K_K+K_M*K_K+K_K+1:end) = theta0(K_K+J*K_K+K_M*K_K+K_K+1:end);        
se_p(K_K+J*K_K+K_M*K_K+K_K+1:end,1) = 0;            
RES_p = [x_p, se_p];

save('res_app_p_upwork.mat', 'RES_p', 'grad_p', 'exitflag_p')
csvwrite('res_app_p_upwork.csv',RES_p)


% Diary
diary off

