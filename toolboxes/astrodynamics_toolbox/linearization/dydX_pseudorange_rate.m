%==========================================================================
%
% dydX_pseudorange_rate  Partial derivatives of pseudorange rate 
% measurement with respect to position and inertial velocity resolved in
% the ECI frame.
%
%   [drhodotdr,drhodotdv] = dydX_pseudorange_rate(r_rcv_eci,r_sat_eci,...
%       v_rcv_eci,v_sat_eci)
%
% Author: Tamas Kis
% Last Update: 2022-03-11
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   r_rcv_eci   - (3×1 double) receiver position resolved in ECI frame [m]
%   r_sat_eci   - (3×n double) GNSS satellite positions resolved in ECI 
%                 frame [m]
%   v_rcv_eci   - (3×1 double) receiver inertial velocity resolved in ECI
%                 frame [m/s]
%   v_sat_eci   - (3×n double) GNSS satellite inertial velocities resolved 
%                 in ECI frame [m/s]
%
% -------
% OUTPUT:
% -------
%   drhodotdr   - (n×3 double) partial derivative of pseudorange rates
%                 w.r.t. position resolved in ECI frame [1/s]
%   drhodotdv   - (n×3 double) partial derivative of pseudorange rates
%                 w.r.t. inertial velocity resolved in ECI frame [-]
%
%==========================================================================
function [drhodotdr,drhodotdv] = dydX_pseudorange_rate(r_rcv_eci,...
    r_sat_eci,v_rcv_eci,v_sat_eci)
    
    % determines number of GNSS satellites
    n = size(r_sat_eci,2);

    % preallocates arrays
    drhodotdr = zeros(n,3);
    drhodotdv = zeros(n,3);

    % partial derivatives for measurements w.r.t. each GNSS satellite
    % resolved in the ECEF frame
    for k = 1:n
        
        % relative position of GNSS satellite w.r.t. receiver resolved in
        % ECI frame [m]
        r_rel = r_sat_eci(:,k)-r_rcv_eci;
        
        % magnitude of relative position vector [m]
        r_rel_mag = inorm(r_rel);
    
        % relative velocity of GNSS satellite w.r.t. receiver [m/s]
        v_rel = v_sat_eci(:,k)-v_rcv_eci;
        
        % partial derivative of pseudorange rate w.r.t. position resolved 
        % in ECI frame [1/s]
        drhodotdr(k,:) = ((v_rel.')*((r_rel*r_rel.')-eye(3,3)*...
            r_rel_mag^2))/(r_rel_mag^3);
        
        % partial derivative of pseudorange w.r.t. inertial velocity 
        % resolved in ECI frame[-]
        drhodotdv(k,:) = -(r_rel/r_rel_mag).';

    end
    
end