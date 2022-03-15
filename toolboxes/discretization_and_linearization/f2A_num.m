%==========================================================================
%
% f2A_num  Continuous dynamics Jacobian from continuous dynamics equation
% (numerical approximation).
%
%   A = f2A_num(f,x)
%   A = f2A_num(f,x,u)
%   A = f2A_num(f,xe)
%   A = f2A_num(f,xe,ue)
%   A = f2A_num(__,t)
%
% See also f2B_num, h2C_num, clinsys_num.
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
%   f       - (1×1 function_handle) continuous dynamics equation,
%             dx/dt = f(x,u,t) (f : ℝⁿ×ℝᵐ×ℝ → ℝⁿ)
%   x       - (n×1 double) state vector
%   u       - (m×1 double) (OPTIONAL) control input
%   t       - (1×1 double) (OPTIONAL) time
%
% ---------------------------------------------------
% INPUT (linearization about equilbrium state/input):
% ---------------------------------------------------
%   f       - (1×1 function_handle) continuous dynamics equation,
%             dx/dt = f(x,u,t) (f : ℝⁿ×ℝᵐ×ℝ → ℝⁿ)
%   xe      - (n×1 double) equilibrium state vector
%   ue      - (m×1 double) (OPTIONAL) equilibrium control input
%   t       - (1×1 double) (OPTIONAL) time
%
% -------
% OUTPUT:
% -------
%   A       - (n×n double) continuous dynamics Jacobian
%
%==========================================================================
function A = f2A_num(f,x,u,t)
    
    % defaults control input to empty vector if not specified
    if (nargin < 3)
        u = [];
    end
    
    % defaults time to empty vector if not specified
    if (nargin < 4)
        t = [];
    end
    
    % evaluates continuous dynamics Jacobian
    A = ijacobian(@(x)f(x,u,t),x);
    
end