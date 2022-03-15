%% RUN_PROPAGATOR.m
% Astrodynamics Toolbox
%
% Runs an orbit propagator to simulate an orbit.
%
% Author: Tamas Kis
% Last Update: 2022-03-10



%% SCRIPT SETUP

% clears Workspace and Command Window, closes all figures
clear; clc; close all;

% adds path to entire directory
addpath(genpath('..'));



%% PROPAGATOR SETTINGS

% ------------
% Data saving.
% ------------

% set file name for simulation data
file_name = 'deputy_simdata';

% produce GPS data (true or false)?
produce_GPS_data = false;

% --------------------
% Satellite selection.
% --------------------

sat = deputy_parameters;

% -----
% Time.
% -----

% initial UTC [y,mo,d,h,m,s]
UTC_start = [2017,1,1,0,0,0];

% simulation duration [h]
duration = 2;

% --------------------------------
% Integator (ODE solver) settings.
% --------------------------------

% time step [s]
time_step = 10;

% integrator ('RK4' or 'ABM8')
integrator = 'ABM8';

% simulation start time [s]
sim_start = 0;

% -------
% Models.
% -------

% density model ('Exponential', 'Harris-Priester', 'Jacchia-Bowman 2008',
% 'Jacchia-Roberts', 'NRLMSISE-00', or 'NRLMSISE-00 MATLAB')
models.density = 'NRLMSISE-00';

% maximum degree/order for gravity model
models.grav_N = 120;

% Earth orientation model ('IAU2006/2000' or 'Simple Rotation')
models.orientation = 'IAU2006/2000';

% Greenwich mean sidereal time model ('Linear', 'Vallado', or 
% 'IAU2006/2000') --> only used if "orientation" set to 'Simple Rotation'
models.gmst = 'N/A';

% --------------
% Perturbations.
% --------------

perturb.drag = true;
perturb.relativity = true;
perturb.srp = true;
perturb.moon = true;
perturb.sun = true;
perturb.empirical = false;



%% SIMULATION

% initialize propagator
prop = initialize_propagator(models,perturb,UTC_start,duration,...
    time_step,integrator,sim_start);

% run simulation
simdata = simulate_orbit(sat,prop);



%% CREATE GPS CONSTELLATION OBJECT FOR A TWO-SPACECRAFT SWARM

if produce_GPS_data
    gps = GPS_Constellation(simdata,2);
end



%% SAVING DATA

% creates "simdata" directory if needed
simdata_file_path = 'simdata/';
if ~exist(simdata_file_path,'dir')
    mkdir(simdata_file_path)
end

% saves simulation data
save(strcat(simdata_file_path,file_name,'.mat'),'simdata');

% saves GPS constellation object
if produce_GPS_data
    save(strcat(simdata_file_path,'gps','.mat'),'gps');
end