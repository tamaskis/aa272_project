%==========================================================================
%
% continuous_dynamics  State vector derivative from current state.
%
%   dxdt = continuous_dynamics(x,t,prop,chief,deputy)
%
% Author: Tamas Kis
% Last Update: 2022-03-10
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   x       - (16×1 double) state vector
%               -->   1-3. rc_eci - chief position resolved in ECI frame 
%                                   [m]
%               -->   4-6. vc_eci - chief inertial velocity resolved in ECI
%                                   frame [m/s]
%               -->     7. bc     - chief clock bias [m]
%               -->     8. bc_dot - chief clock bias drift rate [m/s]
%               -->  9-11. rd_eci - deputy position resolved in ECI frame
%                                   [m]
%               --> 11-14. vd_eci - deputy inertial velocity resolved in
%                                   ECI frame [m/s]
%               -->    15. bd     - deputy clock bias [m]
%               -->    16. bd_dot - deputy clock bias drift rate [m/s]
%   t       - (1×1 double) simulation time [s]
%   prop    - (1×1 struct) propagator settings (see 
%                 "initialize_propagator" function for full definition)
%   chief   - (1×1 struct) chief satellite parameters
%       • m  	- (1×1 double) mass [kg]
%    	• AD   	- (1×1 double) atmospheric drag reference area [m^2]
%    	• CD  	- (1×1 double) drag coefficient [-]
%    	• Asrp  - (1×1 double) solar radiation pressure reference area 
%                 [m^2]
%     	• CR    - (1×1 double) coefficient of reflectivity [-]
%   deputy  - (1×1 struct) deputy satellite parameters
%       • m  	- (1×1 double) mass [kg]
%    	• AD   	- (1×1 double) atmospheric drag reference area [m^2]
%    	• CD  	- (1×1 double) drag coefficient [-]
%    	• Asrp  - (1×1 double) solar radiation pressure reference area 
%                 [m^2]
%     	• CR    - (1×1 double) coefficient of reflectivity [-]
%
% -------
% OUTPUT:
% -------
%   dxdt    - (16×1 double) state vector derivative
%
%==========================================================================
function dxdt = continuous_dynamics(x,t,prop,chief,deputy)
    
    % ----------------------------------------------------
    % Extracts required state variables from state vector.
    % ----------------------------------------------------
    
    % chief ECI state
    Xc = x(1:6);
    
    % chief clock bias drift rate [m/s]
    bc_dot = x(8);
    
    % deputy ECI state
    Xd = x(9:14);
    
    % deputy clock bias drift rate [m/s]
    bd_dot = x(16);
    
    % ----------------------
    % ECI state derivatives.
    % ----------------------

    % chief ECI state derivative
    dXcdt = newton_propagator(Xc,t,prop,chief);

    % deputy ECI state derivative
    dXddt = newton_propagator(Xd,t,prop,deputy);

    % ----------------------------------
    % Assembles state vector derivative.
    % ----------------------------------

    dxdt = [dXcdt;
            bc_dot;
            0;
            dXddt;
            bd_dot;
            0];
    
end