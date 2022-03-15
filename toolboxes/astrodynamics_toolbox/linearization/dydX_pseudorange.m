%==========================================================================
%
% dydX_pseudorange  Partial derivatives of pseudorange measurement with 
% respect to position and inertial velocity resolved in the ECI frame.
%
%   [drhodr,drhodv] = dydX_pseudorange(r_rcv_ecef,r_sat_ecef,R_eci2ecef,...
%       w_eci)
%
% Author: Tamas Kis
% Last Update: 2022-03-10
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   r_rcv_ecef  - (3×1 double) receiver position resolved in ECEF frame [m]
%   r_sat_ecef  - (3×n double) GNSS satellite positions resolved in ECEF
%                 frame [m]
%   R_eci2ecef  - (3×3 double) rotation matrix (ECI --> ECEF)
%   w_eci       - (3×1 double) Earth angular velocity resolved in ECI frame
%                 [rad/s]
%
% -------
% OUTPUT:
% -------
%   drhodr      - (n×3 double) partial derivative of pseudoranges w.r.t.
%                 position resolved in ECI frame [-]
%   drhodv      - (n×3 double) partial derivative of pseudoranges w.r.t.
%                 inertial velocity resolved in ECI frame [s]
%
%==========================================================================
function [drhodr,drhodv] = dydX_pseudorange(r_rcv_ecef,r_sat_ecef,...
    R_eci2ecef,w_eci)
    
    % partial derivative of ECEF position w.r.t. ECI position [-]
    drdr = dECEFdECI(R_eci2ecef,w_eci);

    % determines number of GNSS satellites
    n = size(r_sat_ecef,2);

    % partial derivative of pseudoranges w.r.t. position resolved in ECEF
    % frame [-]
    drhodr_ecef = zeros(n,3);
    for k = 1:n
        drhodr_ecef(k,:) = -((r_sat_ecef(:,k)-r_rcv_ecef)/...
            inorm(r_sat_ecef(:,k)-r_rcv_ecef)).';
    end

    % partial derivative of pseudorange w.r.t. position resolved in ECI 
    % frame [s]
    drhodr = drhodr_ecef*drdr;

    % partial derivative of pseudorange w.r.t. inertial velocity resolved 
    % in ECI frame [s]
    drhodv = zeros(n,3);
    
end