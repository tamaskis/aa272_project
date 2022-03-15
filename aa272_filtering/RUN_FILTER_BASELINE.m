%% RUN_FILTER_BASELINE.m 
% Astrodynamics Toolbox
%
% Runs the baseline filter for the AA 272 project.
%
% Author: Tamas Kis
% Last Update: 2022-03-14



%% SCRIPT SETUP

% clears Workspace and Command Window, closes all figures
clear; clc; close all;

% adds path to entire directory
addpath(genpath('..'));



%% SELECT SIMULATION DATA / SATELLITE

% -----------
% File paths.
% -----------

% path to folder storing data
data_path = 'simdata';

% file name for chief and deputy simulation data
file_name_chief = 'chief_simdata';
file_name_deputy = 'deputy_simdata';

% file name for GPS data
file_name_gps = 'gps';

% --------------------
% Satellite selection.
% --------------------

chief = chief_parameters;
deputy = deputy_parameters;



%% FILTER SETTINGS

% -------------------------------------
% Number of GPS satellites to consider.
% -------------------------------------

l = 8;

% -------
% Models.
% -------

% density model ('Exponential', 'Harris-Priester', 'Jacchia-Bowman 2008',
% 'Jacchia-Roberts', 'NRLMSISE-00', or 'NRLMSISE-00 MATLAB')
models.density = 'NaN';

% maximum degree/order for gravity model
models.grav_N = 20;

% Earth orientation model ('IAU2006/2000' or 'Simple Rotation')
models.orientation = 'IAU2006/2000';

% Greenwich mean sidereal time model ('Linear', 'Vallado', or 
% 'IAU2006/2000') --> only used if "orientation" set to 'Simple Rotation'
models.gmst = 'NaN';

% --------------
% Perturbations.
% --------------

perturb.drag = false;
perturb.relativity = false;
perturb.srp = false;
perturb.moon = false;
perturb.sun = false;
perturb.empirical = false;

% --------------
% Process noise.
% --------------

% position [m] and velocity [m/s] noise standard deviations
sigma_r = 0.5;
sigma_v = 0.005;

% clock bias [m] and clock bias drift rate [m/s] noise standard deviations
sigma_b = 10;
sigma_b_dot = 0.3;

% process noise covariance [m^2][m^2/s^2]
Q = construct_Q(sigma_r,sigma_v,sigma_b,sigma_b_dot);

% ------------------
% Measurement noise.
% ------------------

% pseudorange [m] and pseudorange rate [m/s] noise standard deviations
sigma_rho = 5;
sigma_rho_dot = 0.4;

% measurement noise covariance [m^2][m^2/s^2]
R = construct_R_baseline(sigma_rho,sigma_rho_dot,l);

% -------------------------------------------
% Initial conditions and integrator settings.
% -------------------------------------------

% estimated initial conditions
x0_est = [chief.ECI;
          0;
          0.01;
          deputy.ECI;
          0;
          0.01];
P0 = Q*100;

% sample initial condition from prior distribution
x0 = mvnrnd(x0_est,P0).';

% integrator ('RK1_Euler', 'RK2', 'RK2_heun', 'RK2_ralston', 'RK3',
% 'RK3_heun', 'RK3_ralston', 'SSPRK3', 'RK4', 'RK4_ralston', 'RK4_38')
integrator = 'RK4';



%% FILTERING

% ----------
% Load data.
% ----------

% load simulation data
chief_simdata = load_simdata(strcat(data_path,'/',file_name_chief,'.mat'));
deputy_simdata = load_simdata(strcat(data_path,'/',file_name_deputy,...
    '.mat'));

% load GPS data
gps = load_simdata(strcat(data_path,'/',file_name_gps,'.mat'));

% --------------------------------------------
% Initialize propagator parameters for filter.
% --------------------------------------------

% initial UTC [y,mo,d,h,m,s]
UTC_start = chief_simdata.cal_UTC(1,:);

% simulation duration [h]
duration = (chief_simdata.t(end)-chief_simdata.t(1))/3600;

% time step [s]
dt = chief_simdata.t(2)-chief_simdata.t(1);

% simulation start time [s]
t0 = chief_simdata.t(1);

% initialize propagator parameters
prop = initialize_propagator(models,perturb,UTC_start,duration,dt,...
    integrator,t0);

% -----------------------
% Continuous-time system.
% -----------------------

% continuous nonlinear dynamics equation
f = @(x,u,t) continuous_dynamics(x,t,prop,chief,deputy);

% continuous nonlinear measurement equation
h = @(x,t) measurement_model_baseline(x,t,false,false,prop,gps,l);

% continuous Jacobians
A = @(x,u,t) dynamics_jacobian(x,t,prop,chief,deputy);
C = @(x,t) measurement_jacobian_baseline(x,t,prop,gps,l);

% OPTION FOR NUMERICAL LINEARIZATION
%F = @(x,u,k) f2stm_num(f,x,u,k2t(k,dt),dt);
%H = @(x,k) hd2H_num(hd,x,k);

% ---------------
% Discretization.
% ---------------

% discrete nonlinear dynamics equation
fd = f2fd_num(f,dt,t0,integrator);

% discrete nonlinear measurement equation
hd = h2hd_num(h,dt,t0);

% discrete Jacobians
F = @(x,u,k) Af2stm(A,f,x,u,k2t(k,dt),dt);
H = C2H_num(C,dt);

% ----------------------
% Generate ground truth.
% ----------------------

% true states
x_true = true_state(chief_simdata,deputy_simdata);

% true measurements
y = true_measurement_baseline(x_true,chief_simdata.t,prop,gps,l);

% -----------------
% State estimation.
% -----------------

% runs extended Kalman filter
[x,P,tsol,rank_Ob,z_pre,z_post] = EKF_sim(fd,hd,F,H,Q,R,[],y,x0,P0,true);



%% TRANSFORM ABSOLUTE STATE ESTIMATES TO RELATIVE STATE ESTIMATES

% -------------
% Ground truth.
% -------------

% ground truth deputy ECI state [m][m/s]
Xd_true = x_true(9:14,:);

% ground truth chief ECI state [m][m/s]
Xc_true = x_true(1:6,:);

% ground truth deputy relative ECI state resolved in chief's RTN frame 
% [m][m/s]
dX_true = abs2rel_filter_results(Xd_true,Xc_true);

% -----------------
% Filter estimates.
% -----------------

% estimated deputy and chief ECI states
Xd = x(9:14,:);
Xc = x(1:6,:);

% estimated deputy and chief ECI state error covariances
Pd = P(9:14,9:14,:);
Pc = P(1:6,1:6,:);

% estimated deputy relative ECI state [m][m/s] and associated error 
% covariance [m^2][m^2/s^2] resolved in chief's RTN frame
[dX,dP] = abs2rel_filter_results(Xd,Xc,Pd,Pc);



%% SAMPLE RESIDUAL PLOT

% time vector for plotting
t = chief_simdata.t;

% residual plot for chief pseudorange
figure;
hold on;
plot(t(5:end)/3600,z_pre(1,5:end),'k.');
plot(t(5:end)/3600,z_post(1,5:end),'r.');
hold off;
grid on;
xlabel('time $[\mathrm{h}]$','interpreter','latex','fontsize',18);
ylabel('$\rho^{(k)}$ measurement residuals $[\mathrm{m}]$',...
    'interpreter','latex','fontsize',18);
legend('pre-fit residuals','post-fit residuals','interpreter','latex',...
    'fontsize',14,'location','best');



%% CHIEF CLOCK BIAS PLOT

% 3-sigma bounds for chief clock bias error [m]
[bc_lower,bc_upper] = covariance_bounds(x(7,:),P(7,7,:),3);

% plots
opts.name = '$b_{c}$';
opts.tunits = 'h';
opts.xunits = 'm';
opts.shaded = true;
opts.color = [0,0.4470,0.7410];
opts.error = true;
opts.M = 3;
plot_filter_results(t'/3600,x(7,:),bc_lower,bc_upper,x_true(7,:),opts);



%% PLOTS

% sigma (covariance) bound
M = 3;

% plots results
plot_absolute_results(t,Xc,Pc,Xc_true,M);
plot_absolute_results(t,Xd,Pd,Xd_true,M);
plot_relative_results(t,dX,dP,dX_true,M);