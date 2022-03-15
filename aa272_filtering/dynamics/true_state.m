%==========================================================================
%
% true_state  Ground truth state.
%
%   x_true = true_state(chief_simdata,deputy_simdata)
%
% Author: Tamas Kis
% Last Update: 2022-03-07
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   chief_simdata   - (1×1 struct) chief satellite simulation data
%   deputy_simdata  - (1×1 struct) deputy satellite simulation data
%
% -------
% OUTPUT:
% -------
%   x_true          - (14×N double) ground truth state
%                       -->   1-3. rc_eci - chief position resolved in ECI 
%                                           frame [m]
%                       -->   4-6. vc_eci - chief inertial velocity 
%                                           resolved in ECI frame [m/s]
%                       -->     7. bc     - chief clock bias [m]
%                       -->     8. bc_dot - chief clock bias drift rate 
%                                           [m/s]
%                       -->  9-11. rd_eci - deputy position resolved in ECI
%                                           frame [m]
%                       --> 11-14. vd_eci - deputy inertial velocity 
%                                           resolved in ECI frame [m/s]
%                       -->    15. bd     - deputy clock bias [m]
%                       -->    16. bd_dot - deputy clock bias drift rate 
%                                           [m/s]
%
% -----
% NOTE:
% -----
%   --> N = number of iterations (i.e. length of time vector)
%
%==========================================================================
function x_true = true_state(chief_simdata,deputy_simdata)
    N = size(chief_simdata.r_eci,2);
    x_true = [chief_simdata.r_eci;
              chief_simdata.v_eci;
              200*ones(1,N);
              0.0001*ones(1,N);
              deputy_simdata.r_eci;
              deputy_simdata.v_eci;
              50*ones(1,N);
              0.00005*ones(1,N)];
end