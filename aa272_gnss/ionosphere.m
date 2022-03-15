%==========================================================================
%
% ionosphere  Ionospheric path delay.
%
%   I = ionosphere(GPS_s,r_rcv,r_sat,Kalpha,Kbeta)
%
% Author: Matthew Hunter (adapted from function by Vince Giralo)
% Last Update: 2022-03-07
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   GPS_s   - (1×1 double) seconds of GPS week [s]
%   r_rcv   - (3×1 double) receiver position resolved in ECEF frame [m]
%   r_sat   - (3×1 double) GPS satellite pos. resolved in ECEF frame [m]
%   Kalpha  - (1×4 double) Klobuchar alpha parameters
%   Kbeta   - (1×4 double) Klobuchar beta parameters
%
% -------
% OUTPUT:
% -------
%   I       - (1×1 double) ionospheric path delay [m]
%
%==========================================================================
function I = ionosphere(GPS_s,r_rcv,r_sat,Kalpha,Kbeta)
    
    % -----------------
    % Input processing.
    % -----------------
    
    % azimuth and elevation of GNSS satellite with respect to receiver [°]
    rho_enu = ecef2enu(r_sat,r_rcv);
    [Az,El] = enu2aer(rho_enu);

    % converts azimuth and elevation to radians
    Az = deg2rad(Az);
    El = deg2rad(El);
    
    % geodetic latitude and longitude of receiver [°]
    [lat,lon] = ecef2geod(r_rcv);

    % converts geodetic latitude and longitude to radians
    lat = deg2rad(lat);
    lon = deg2rad(lon);
    
    % converts geodetic latitude, geodetic longitude, and elevation to 
    % semi-circles
    lat = lat/pi;
    lon = lon/pi;
    El = El/pi;

    % ---------------------
    % GPS Ionosphere model.
    % ---------------------

    A1 = 5e-9;
    A3 = 50400;
    PSI = 0.0137/(El+0.11)-0.022;
    PHIl = lon+PSI*cos(Az);
    if PHIl > 0.416
        PHIl = 0.416;
    elseif PHIl < -0.416
        PHIl = -0.416;
    end
    LAMl = lat+(PSI*sin(Az))/cos(PHIl*pi);
    PHIm = PHIl+0.064*cos((LAMl-1.617)*pi);
    t = imod(4.32e4*LAMl+GPS_s,86400);
    if t > 86400
        t = t-86400;
    elseif t < 0
        t = t + 86400;
    end
    F = 1+16*(0.53-El)^3;
    A2 = Kalpha(1)*(PHIm)^0+Kalpha(2)*(PHIm)^1+Kalpha(3)*(PHIm)^2+...
        Kalpha(4)*(PHIm)^3;
    if A2 < 0
        A2 = 0;
    end
    A4 = Kbeta(1)*(PHIm)^0+Kbeta(2)*(PHIm)^1+Kbeta(3)*(PHIm)^2+Kbeta(4)*...
        (PHIm)^3;
    if A4 >= 172800
        A4 = 172800;
    elseif A4 < 72000
        A4 = 72000;
    end
    x = 2*pi*(t-A3)/A4;
    if abs(x) > 1.57
        T = F*A1;
    else
        T = F*(A1 + A2*(1-x^2/2 + x^4/24));
    end
    I = T*C_LIGHT;
end