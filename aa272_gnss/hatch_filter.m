%==========================================================================
%
% hatch_filter  Hatch filter for smoothing GNSS pseudorange measurements.
%
%   srho = hatch_filter(rho,srho_prev,phi,phi_prev,M)
%
% Author: Matthew Hunter
% Last Update: 2022-03-07
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   rho         - (1×1 double) pseudorange at current sample time [m]
%   srho_prev   - (1×1 double) smoothed pseudorange at previous sample time
%                 [m]
%   phi         - (1×1 double) carrier phase measurement at current sample 
%                 time [m]
%   phi_prev    - (1×1 double) carrier phase measurement at previous sample 
%                 time [m]
%   M           - (1×1 double) current iteration number of hatch filter [-]
%
% -------
% OUTPUT:
% -------
%   srho        - (1×1 double) smoothed pseudorange at current sample time 
%                 [m]
%
%==========================================================================
function srho = hatch_filter(rho,srho_prev,phi,phi_prev,M)
    srho = (rho+(M-1)*(srho_prev+phi-phi_prev))/M;
end