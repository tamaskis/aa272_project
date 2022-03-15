%==========================================================================
%
% covariance_bounds  Lower and upper sigma bounds for a time series of 
% means and covariances.
%
%   [mu_lower,mu_upper] = covariance_bounds(mu,Sigma,M)
%
% Author: Tamas Kis
% Last Update: 2021-12-09
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   mu          - (n×N double) time history of mean (μ)
%   Sigma       - (n×n×N double) time history of covariance (Σ)
%   M           - (OPTIONAL) (1×1 double) number of standard deviations for
%                 bounds (defaults to 1)
%                   --> for example, M = 2 returns ±2σ
%
% -------
% OUTPUT:
% -------
%   mu_lower    - (n×N double) lower Mσ covariance bound
%   mu_upper    - (n×N double) upper Mσ covariance bound
%
% -----
% NOTE:
% -----
%   --> "mu", "mu_lower", and "mu_upper" are n×N matrices of n×1 vectors.
%   --> "Sigma" is an n×n×N array of n×n matrices.
%
%==========================================================================
function [mu_lower,mu_upper] = covariance_bounds(mu,Sigma,M)
    
    % defaults M to 1 if not specified
    if nargin == 2, M = 1; end

    % preallocates matrices
    mu_lower = zeros(size(mu));
    mu_upper = zeros(size(mu));
    
    % finds lower and upper N-sigma bounds
    for i = 1:size(mu,2)
        mu_lower(:,i) = mu(:,i)-M*sqrt(diag(Sigma(:,:,i)));
        mu_upper(:,i) = mu(:,i)+M*sqrt(diag(Sigma(:,:,i)));
    end
    
end