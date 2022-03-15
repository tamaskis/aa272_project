%==========================================================================
%
% h2hd_num  Discretization of continuous nonlinear measurement equation
% (numerical approximation).
%
%   hd = h2hd_num(h,dt)
%   hd = h2hd_num(h,dt,t0)
%
% See also f2fd_num.
%
% Author: Tamas Kis
% Last Update: 2022-02-25
%
% REFERENCES:
%   [1] TODO
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   h       - (1×1 function_handle) continuous nonlinear measurement 
%             equation, y = h(x,t) (h : ℝⁿ×ℝ → ℝᵖ)
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
function hd = h2hd_num(h,dt,t0)
    
    % defaults initial time to 0 if not specified
    if (nargin < 3) || isempty(t0)
        t0 = 0;
    end
    
    % discretization
    hd = @(xk,k) h(xk,k2t(k,dt,t0));
    
end