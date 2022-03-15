%==========================================================================
%
% graphic  GNSS group and phase ionospheric correction (GRAPHIC)
% measurement.
%
%   rho_gr = graphic(r_rcv_ecef,r_sat_ecef,b_rcv,b_sat,MJD_GPS,N)
%   rho_gr = graphic(r_rcv_ecef,r_sat_ecef,b_rcv,b_sat,MJD_GPS,N,
%       iono_delay,noise)
%   [rho_gr,rho,rho_phi] = graphic(__)
%
% Author: Tamas Kis
% Last Update: 2022-03-14
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   r_rcv_ecef  - (3×1 double) receiver position resolved in ECEF frame [m]
%   r_sat_ecef  - (3×1 double) GPS satellite position resolved in ECEF 
%                 frame [m]
%   b_rcv       - (1×1 double) receiver clock bias [m]
%   b_sat       - (1×1 double) GPS satellite clock bias [m]
%   MJD_GPS     - (1×1 double) GPS time [MJD]
%   N           - (1×1 double) integer ambiguity [-]
%   iono_delay  - (1×1 logical) (OPTIONAL) "true" if ionospheric path delay
%                 should be included, "false" otherwise (defaults to 
%                 "false")
%   noise       - (1×1 logical) (OPTIONAL) "true" if random noise should be
%                 included, "false" otherwise (defaults to false)
%
% -------
% OUTPUT:
% -------
%   rho_gr      - (1×1 double) GRAPHIC measurement [m]
%   rho         - (1×1 double) pseudorange measurement [m]
%   rho_phi     - (1×1 double) carrier phase measurement [m]
%
%==========================================================================
function [rho_gr,rho,rho_phi] = graphic(r_rcv_ecef,r_sat_ecef,b_rcv,...
    b_sat,MJD_GPS,N,iono_delay,noise)
    
    % pseudorange measurement [m]
    rho = pseudorange(r_rcv_ecef,r_sat_ecef,b_rcv,b_sat,MJD_GPS,...
        iono_delay,noise);
    
    % carrier phase measurement [m]
    rho_phi = carrier_phase(r_rcv_ecef,r_sat_ecef,b_rcv,b_sat,MJD_GPS,N,...
        iono_delay,noise);
    
    % GRAPHIC measurement [m]
    rho_gr = (rho+rho_phi)/2;
    
end