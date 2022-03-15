%==========================================================================
%
% pseudorange_rate  GNSS pseudorange rate measurement.
%
%   rho_dot = pseudorange_rate(r_rcv_ecef,r_sat_ecef,v_rcv_ecef,...
%       v_sat_ecef,b_dot_rcv,b_dot_sat)
%   rho_dot = pseudorange_rate(r_rcv_ecef,r_sat_ecef,v_rcv_ecef,...
%       v_sat_ecef,b_dot_rcv,b_dot_sat,noise)
%
% Author: Tamas Kis
% Last Update: 2022-03-09
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   r_rcv_ecef  - (3×1 double) receiver position resolved in ECEF frame [m]
%   r_sat_ecef  - (3×1 double) GPS satellite position resolved in ECECF 
%                 frame [m]
%   v_rcv_ecef  - (3×1 double) receiver ECEF velocity resolved in ECEF 
%                 frame [m/s]
%   v_sat_ecef  - (3×1 double) GPS satellite ECEF velocity resolved in ECEF
%                 frame [m/s]
%   b_dot_rcv   - (1×1 double) receiver clock bias drift rate [m/s]
%   b_dot_sat   - (1×1 double) satellite clock bias drift rate [m/s]
%   noise       - (1×1 logical) "true" if random noise should be included,
%                 "false" otherwise (defaults to false)
%
% -------
% OUTPUT:
% -------
%   rho_dot     - (1×1 double) pseudorange rate [m/s]
%
%==========================================================================
function rho_dot = pseudorange_rate(r_rcv_ecef,r_sat_ecef,v_rcv_ecef,...
    v_sat_ecef,b_dot_rcv,b_dot_sat,noise)
    
    % defaults "noise" to "false" if not input
    if (nargin < 7) || isempty(noise)
        noise = false;
    end

    % receiver-to-satellite line-of-site unit vector
    LOS = (r_sat_ecef-r_rcv_ecef)/inorm(r_sat_ecef-r_rcv_ecef);
        
    % pseudorange rate measurement [m/s]
    rho_dot = idot((v_sat_ecef-v_rcv_ecef),LOS)+(b_dot_rcv-b_dot_sat);

    % adds noise (0 mean, 0.5 m/s std. dev.) to pseudorange rate [m/s]
    if noise
        rho_dot = rho_dot+normrnd(0,0.5);
    end
    
end