%==========================================================================
%
% construct_R_baseline  Constructs the measurement noise covariance matrix
% for the AA 272 project (baseline filter).
%
%   R = construct_R_baseline(sigma_rho,sigma_rho_dot,l)
%
% Author: Tamas Kis
% Last Update: 2022-03-14
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   sigma_rho       - (1×1 double) pseudorange noise standard deviation [m]
%   sigma_rho_dot   - (1×1 double) pseudorange rate noise standard
%                     deviation [m/s]
%   l               - (1×1 double) (OPTIONAL) number of GNSS satellites 
%                     used for measurements (defaults to 4)
%
% -------
% OUTPUT:
% -------
%   R               - ((4l)×(4l) double) measurement noise covariance 
%                     [m^2][m^2/s^2]
%
%==========================================================================
function R = construct_R_baseline(sigma_rho,sigma_rho_dot,l)

    % defaults "l" to 4 if not input
    if (nargin < 3) || isempty(l)
        l = 4;
    end
    
    % pseudorange measurement noise covariance [m^2]
    R_rho = diag(ones(1,l)*sigma_rho^2);
    
    % pseudorange rate measurement noise covariance [m^2/s^2]
    R_rho_dot = diag(ones(1,l)*sigma_rho_dot^2);
    
    % indices to assist construction of measurement noise covariance
    l1 = 1;
    l2 = l;
    l3 = l+1;
    l4 = 2*l;
    l5 = 2*l+1;
    l6 = 3*l;
    l7 = 3*l+1;
    l8 = 4*l;

    % assembles measurement noise covariance
    R = zeros(4*l,4*l);
    R(l1:l2,l1:l2) = R_rho;
    R(l3:l4,l3:l4) = R_rho_dot;
    R(l5:l6,l5:l6) = R_rho;
    R(l7:l8,l7:l8) = R_rho_dot;
    
end