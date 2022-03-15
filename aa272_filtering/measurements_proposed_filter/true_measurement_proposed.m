%==========================================================================
%
% true_measurement_proposed  Ground truth measurements.
%
%   y = true_measurement_proposed(x_true,t,prop,gps)
%   y = true_measurement_proposed(x_true,t,prop,gps,l)
%
% Author: Tamas Kis
% Last Update: 2022-03-14
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   x_true  - (16×K double) ground truth state trajectory
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
%   t       - (K×1 double) time vector
%   prop    - (1×1 struct) propagator settings (see "initialize_propagator"
%             function for full definition)
%   gps     - (1×1 GPS_Constellation) GPS constellation object
%   l       - (1×1 double) (OPTIONAL) number of GNSS satellites used for
%              measurements (defaults to 4)
%
% -------
% OUTPUT:
% -------
%   y       - (4l×K double) ground truth measurement trajectory
%               -->       1-l. "l" GRAPHIC Hatch measurements for chief 
%                              satellite
%               -->  (l+1)-2l. "l" pseudorange rates for chief satellite
%               --> (2l+1)-3l. "l" SDCP Hatch measurements between chief
%                              and deputy satellites
%               --> (3l+1)-4l. "l" pseudorange rates for deputy satellite
%
%==========================================================================
function y = true_measurement_proposed(x_true,t,prop,gps,l)
    
    % defaults "l" to 4 if not input
    if (nargin < 5) || isempty(l)
        l = 4;
    end
    
    % length of time vector
    K = length(t);
    
    % preallocates array to store measurements
    y = zeros(4*l,K);

    % preallocates arrays to store auxiliary measurements
    rho_c = zeros(l,K);
    rho_d = zeros(l,K);
    SVID = zeros(l,K);
    
    % generates ground-truth measurements at each sample time --> no Hatch
    % filtering on these measurements
    for i = 1:K
        [y(:,i),rho_c(:,i),rho_d(:,i),SVID(:,i)] =...
            measurement_model_proposed(x_true(:,i),t(i),true,true,prop,...
            gps,l);
    end
    
    % ---------------------
    % "Hatch" measurements.
    % ---------------------

    % extract GRAPHIC and SDCP measurements
    rho_gr_c = y(1:l,:);
    rho_sdcp = y((2*l+1):(3*l),:);
    
    % run Hatch filter using GRAPHIC measurements in place of carrier phase
    % measurements
    rho_gr_hatch = batch_hatch(rho_c,rho_gr_c,SVID);

    % run Hatch filter using single-difference pseudorange and carrier
    % phase measurements
    rho_sdcp_hatch = batch_hatch(rho_d-rho_c,rho_sdcp,SVID);

    % replace 1st and 3rd block rows of measurement vector time history
    % with "Hatch" measurements
    y(1:l,:) = rho_gr_hatch;
    y((2*l+1):(3*l),:) = rho_sdcp_hatch;

end