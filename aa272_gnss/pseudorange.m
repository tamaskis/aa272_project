%==========================================================================
%
% pseudorange  GNSS pseudorange measurement.
%
%   rho = pseudorange(r_rcv_ecef,r_sat_ecef,b_rcv,b_sat,MJD_GPS,iono_delay)
%   rho = pseudorange(r_rcv_ecef,r_sat_ecef,b_rcv,b_sat,MJD_GPS,...
%       iono_delay,noise)
%
% Author: Josh Geiser
% Last Update: 2022-03-09
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
%   iono_delay  - (1×1 logical) (OPTIONAL) "true" if ionospheric path delay
%                 should be included, "false" otherwise (defaults to 
%                 "false")
%   noise       - (1×1 logical) (OPTIONAL) "true" if random noise should be
%                 included, "false" otherwise (defaults to false)
%
% -------
% OUTPUT:
% -------
%   rho         - (1×1 double) pseudorange [m]
%
%==========================================================================
function rho = pseudorange(r_rcv_ecef,r_sat_ecef,b_rcv,b_sat,MJD_GPS,...
    iono_delay,noise)
    
    % defaults "iono_delay" to "false" if not input
    if (nargin < 6) || isempty(iono_delay)
        iono_delay = false;
    end

    % defaults "noise" to "false" if not input
    if (nargin < 7) || isempty(noise)
        noise = false;
    end

    % GPS weeks [wk] and seconds [s]
    [~,GPS_s] = gps2wks(MJD_GPS);
    
    % Klobuchar coefficients
    [Kalpha,Kbeta] = klobuchar_coeffs;
    
    % ionospheric path delay [m]
    if iono_delay
        I = ionosphere(GPS_s,r_rcv_ecef,r_sat_ecef,Kalpha,Kbeta);
    else
        I = 0;
    end
    
    % pseudorange measurement [m]
    rho = inorm(r_sat_ecef-r_rcv_ecef)+(b_rcv-b_sat)+I;

    % adds noise (0 mean, 0.5 m std. dev.) to pseudorange measurement [m]
    if noise
        rho = rho+normrnd(0,0.5);
    end
    
end