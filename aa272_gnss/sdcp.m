%==========================================================================
%
% sdcp  GNSS single-difference carrier phase (SDCP) measurement.
%
%   rho_sdcp = sdcp(rd_ecef,rc_ecef,r_sat_ecef,bd,bc,b_sat,MJD_GPS,Nd,Nc)
%   rho_sdcp = sdcp(rd_ecef,rc_ecef,r_sat_ecef,bd,bc,b_sat,MJD_GPS,Nd,...
%       Nc,iono_delay,noise)
%   [rho_sdcp,rho_phi_d,rho_phi_c] = sdcp(__)
%
% Author: Matthew Hunter
% Last Update: 2022-03-14
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   rd_ecef     - (3×1 double) deputy position resolved in ECEF frame [m]
%   rc_ecef     - (3×1 double) chief position resolved in ECEF frame [m]
%   r_sat_ecef  - (3×1 double) GPS satellite position resolved in ECEF 
%                 frame [m]
%   bd          - (1×1 double) deputy clock bias [m]
%   bc          - (1×1 double) chief clock bias [m]
%   b_sat       - (1×1 double) GPS satellite clock bias [m]
%   MJD_GPS     - (1×1 double) GPS time [MJD]
%   Nd          - (1×1 double) integer ambiguity between deputy and GPS 
%                 satellite [-]
%   Nc          - (1×1 double) integer ambiguity between chief and GPS 
%                 satellite [-]
%   iono_delay  - (1×1 logical) (OPTIONAL) "true" if ionospheric path delay
%                 should be included, "false" otherwise (defaults to 
%                 "false")
%   noise       - (1×1 logical) (OPTIONAL) "true" if random noise should be
%                 included, "false" otherwise (defaults to false)
%
% -------
% OUTPUT:
% -------
%   rho_sdcp    - (1×1 double) SDCP measurement [m]
%   rho_phi_d   - (1×1 double) deputy carrier phase measurement [m]
%   rho_phi_c   - (1×1 double) chief carrier phase measurement [m]
%
%==========================================================================
function [rho_sdcp,rho_phi_d,rho_phi_c] = sdcp(rd_ecef,rc_ecef,...
    r_sat_ecef,bd,bc,b_sat,MJD_GPS,Nd,Nc,iono_delay,noise)
    
    % deputy carrier phase measurement [m]
    rho_phi_c = carrier_phase(rc_ecef,r_sat_ecef,bc,b_sat,MJD_GPS,Nc,...
        iono_delay,noise);

    % chief carrier phase measurement [m]
    rho_phi_d = carrier_phase(rd_ecef,r_sat_ecef,bd,b_sat,MJD_GPS,Nd,...
        iono_delay,noise);

    % single difference carrier phase (SDCP) measurement [m]
    rho_sdcp = rho_phi_d-rho_phi_c;
    
end