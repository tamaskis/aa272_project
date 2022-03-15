%==========================================================================
%
% carrier_phase  GNSS carrier phase measurement.
%
%   rho_phi = carrier_phase(r_rcv_ecef,r_sat_ecef,b_rcv,b_sat,...
%       MJD_GPS,N,iono_delay,noise)
%
% Author: Matthew Hunter
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
%   rho_phi     - (1×1 double) carrier phase measurement [m]
%
%==========================================================================
function rho_phi = carrier_phase(r_rcv_ecef,r_sat_ecef,b_rcv,b_sat,...
    MJD_GPS,N,iono_delay,noise)

    % defaults "iono_delay" to "false" if not input
    if (nargin < 7) || isempty(iono_delay)
        iono_delay = false;
    end

    % defaults "noise" to "false" if not input
    if (nargin < 8) || isempty(noise)
        noise = false;
    end

    % GPS weeks [wk] and seconds [s]
    [~,GPS_s] = gps2wks(MJD_GPS);
    
    % L1 signal wavelength [m]
    lambda = C_LIGHT/F_L1;
    
    % Klobuchar coefficients
    [Kalpha,Kbeta] = klobuchar_coeffs;

    % ionospheric path delay [m]
    if iono_delay
        I = ionosphere(GPS_s,r_rcv_ecef,r_sat_ecef,Kalpha,Kbeta);
    else
        I = 0;
    end
    
    % carrier phase measurement [m]
    rho_phi = inorm(r_sat_ecef-r_rcv_ecef)+(b_rcv-b_sat)-I+lambda*N;

    % adds noise (0 mean, 0.005 m std. dev.) to carrier phase meas. [m]
    if noise
        rho_phi = rho_phi+normrnd(0,0.005);
    end
    
end