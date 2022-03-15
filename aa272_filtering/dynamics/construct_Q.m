%==========================================================================
%
% construct_Q  Constructs the process noise covariance matrix for the
% AA 272 project.
%
%   Q = construct_Q(sigma_r,sigma_v,sigma_b,sigma_b_dot)
%
% Author: Tamas Kis
% Last Update: 2022-03-14
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   sigma_r     - (1×1 double) pseudorange noise standard deviation [m]
%   sigma_v     - (1×1 double) pseudorange rate error standard deviation 
%                 [m/s]
%   sigma_b     - (1×1 double) clock bias noise standard deviation [m]
%   sigma_b_dot - (1×1 double) clock bias drift rate noise standard 
%                 deviation [m/s]
%
% -------
% OUTPUT:
% -------
%   Q           - (16×16 double) process noise covariance [m^2][m^2/s^2]
%
%==========================================================================
function Q = construct_Q(sigma_r,sigma_v,sigma_b,sigma_b_dot)
    
    % position noise covariance [m^2]
    Q_r = diag(ones(1,3)*sigma_r^2);

    % velocity noise covariance [m^2/s^2]
    Q_v = diag(ones(1,3)*sigma_v^2);

    % ECI state noise covariance [m^2][m^2/s^2]
    Q_rv = zeros(6,6);
    Q_rv(1:3,1:3) = Q_r;
    Q_rv(4:6,4:6) = Q_v;

    % clock bias noise covariance [m^2]
    Q_b = sigma_b^2;

    % clock bias drift rate noise covariance [m^2/s^2]
    Q_b_dot = sigma_b_dot^2;
    
    % assembles process covariance
    Q = zeros(16,16);
    Q(1:6,1:6) = Q_rv;
    Q(7,7) = Q_b;
    Q(8,8) = Q_b_dot;
    Q(9:14,9:14) = Q_rv;
    Q(15,15) = Q_b;
    Q(16,16) = Q_b_dot;
    
end