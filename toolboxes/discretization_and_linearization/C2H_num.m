%==========================================================================
%
% C2H_num  Discretization of continuous measurement Jacobian.
%
%   H = C2H_num(C,x,t,dt)
%   H = C2H_num(C,x,t,dt,t0)
%
% Author: Tamas Kis
% Last Update: 2022-03-10
%
% REFERENCES:
%   [1] TODO
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   C       - (1×1 function_handle) continuous measurement Jacobian, C(x,t)
%             (C : ℝⁿ×ℝ → ℝᵖˣⁿ)
%   dt      - (1×1 double) time step
%   t0      - (1×1 double) (OPTIONAL) initial time (defaults to 0)
%
% -------
% OUTPUT:
% -------
%   hd      - (1×1 function_handle) discrete nonlinear measurement 
%             equation, yₖ = hd(xₖ,k) (hd : ℝⁿ×ℤ → ℝᵖ)
%
%==========================================================================
function H = C2H_num(C,dt,t0)
    
    % defaults initial time to 0 if not specified
    if (nargin < 3) || isempty(t0)
        t0 = 0;
    end
    
    % discretization
    H = @(xk,k) C(xk,k2t(k,dt,t0));
    
end