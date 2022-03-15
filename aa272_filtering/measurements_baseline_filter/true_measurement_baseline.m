%==========================================================================
%
% true_measurement_baseline  Ground truth measurements.
%
%   y = true_measurement_baseline(x_true,t,prop,gps)
%   y = true_measurement_baseline(x_true,t,prop,gps,l)
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
%               -->       1-l. "l" pseudoranges for chief satellite
%               -->  (l+1)-2l. "l" pseudorange rates for chief satellite
%               --> (2l+1)-3l. "l" pseudoranges for deputy satellite
%               --> (3l+1)-4l. "l" pseudorange rates for deputy satellite
%
%==========================================================================
function y = true_measurement_baseline(x_true,t,prop,gps,l)
    
    % defaults "l" to 4 if not input
    if (nargin < 5) || isempty(l)
        l = 4;
    end
    
    % length of time vector
    K = length(t);
    
    % preallocates array to store measurements
    y = zeros(4*l,K);
    
    % generates ground-truth measurements at each sample time
    for i = 1:K
        y(:,i) = measurement_model_baseline(x_true(:,i),t(i),true,true,...
            prop,gps,l);
    end
    
end