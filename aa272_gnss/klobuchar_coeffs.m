%==========================================================================
%
% klobuchar_coeffs  Klobuchar coefficients used for calculating the
% ionospheric path delay.
%
%   [Kalpha,Kbeta] = klobuchar_coeffs
%
% Author: Matthew Hunter
% Last Update: 2022-03-06
%
%--------------------------------------------------------------------------
%
% -------
% OUTPUT:
% -------
%   Kalpha  - (1×4 double) Klobuchar alpha parameters
%   Kbeta   - (1×4 double) Klobuchar beta parameters
%
%==========================================================================
function [Kalpha,Kbeta] = klobuchar_coeffs
    Kalpha = [0.2142e-7,0.7451e-8,-0.1192e-6,0];
    Kbeta = [0.1229e6,0,-0.2621e6,0.1966e6];
end