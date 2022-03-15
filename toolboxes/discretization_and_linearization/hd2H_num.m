%==========================================================================
%
% hd2H_num  Discrete measurement Jacobian from discrete measurement
% equation (numerical approximation).
%
%   H = hd2H_num(hd,xk)
%   H = hd2H_num(hd,xe)
%   H = hd2H_num(__,k)
%
% See also fd2F_num, fd2G_num, dlinsys_num.
%
% Author: Tamas Kis
% Last Update: 2022-02-25
%
% REFERENCES:
%   [1] TODO
%
%--------------------------------------------------------------------------
%
% -------------------------------------------------------------
% INPUT (linearization about nominal or estimated state/input):
% -------------------------------------------------------------
%   hd      - (1×1 function_handle) discrete measurement equation,
%             yₖ = hd(xₖ,k) (hd : ℝⁿ×ℤ → ℝᵖ)
%   xk      - (n×1 double) state vector at kth sample
%   k       - (1×1 double) (OPTIONAL) sample number
%
% ---------------------------------------------------
% INPUT (linearization about equilbrium state/input):
% ---------------------------------------------------
%   hd      - (1×1 function_handle) discrete measurement equation,
%             yₖ = hd(xₖ,k) (hd : ℝⁿ×ℤ → ℝᵖ)
%   xe      - (n×1 double) equilibrium state vector
%   k       - (1×1 double) (OPTIONAL) sample number
%
% -------
% OUTPUT:
% -------
%   H       - (p×n sym) discrete measurement Jacobian at kth sample
%
%==========================================================================
function H = hd2H_num(hd,x,k)
    
    % defaults sample number to empty vector if not specified
    if (nargin < 3)
        k = [];
    end
    
    % evaluates discrete measurement Jacobian
    H = ijacobian(@(x)hd(x,k),x);
    
end