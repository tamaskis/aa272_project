%==========================================================================
%
% dynamics_jacobian  Continuous dynamics Jacobian from current state.
%
%   A = dynamics_jacobian(x,t,prop,chief,deputy)
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
%   A       - (16×16 double) continuous dynamics Jacobian
%
%==========================================================================
function A = dynamics_jacobian(x,t,prop,chief,deputy)
    
    % -------------------------------------------------
    % Extract required variables from input parameters.
    % -------------------------------------------------
    
    % chief ECI state
    Xc = x(1:6);
    
    % deputy ECI state
    Xd = x(9:14);

    % --------------------------------------------------------
    % Continuous dynamics Jacobians for variational equations.
    % --------------------------------------------------------

    % dynamics Jacobian for chief variational equations
    Ac_var = A_var_eqns(Xc,t,prop,chief);

    % dynamics Jacobian for deputy variational equations
    Ad_var = A_var_eqns(Xd,t,prop,deputy);

    % ---------------------------------------
    % Assembles continuous dynamics Jacobian.
    % ---------------------------------------

    % initializes dynamics Jacobian as identity matrix
    A = zeros(16);

    % edits blocks/elements of A
    A(1:6,1:6) = Ac_var;
    A(7,8) = 1;
    A(9:14,9:14) = Ad_var;
    A(15,16) = 1;
    
end